/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Model/ModelOptimSillsVario.hpp"

#include "Variogram/Vario.hpp"
#include "Model/Model.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"
#include "Model/ModelOptimVario.hpp"
#include "Model/Constraints.hpp"

#define IJDIR(ijvar, ipadir)    ((ijvar)*_npadir + (ipadir))
#define _WT(ijvar, ipadir)      _wt[IJDIR(ijvar, ipadir)]
#define _GG(ijvar, ipadir)      _gg[IJDIR(ijvar, ipadir)]
#define _WT2(ijvar, ipadir)     _wt2[IJDIR(ijvar, ipadir)]
#define _GG2(ijvar, ipadir)     _gg2[IJDIR(ijvar, ipadir)]
#define TAB(ijvar, ipadir)      tabin[IJDIR(ijvar, ipadir)]
#define DD(idim, ijvar, ipadir) _dd[idim][IJDIR(ijvar, ipadir)]

#define CORRECT(idir, k)                                                       \
  (!isZero(vario->getHhByIndex(idir, k)) &&                                    \
   !FFFF(vario->getHhByIndex(idir, k)) &&                                      \
   !isZero(vario->getSwByIndex(idir, k)) &&                                    \
   !FFFF(vario->getSwByIndex(idir, k)) && !FFFF(vario->getGgByIndex(idir, k)))
#define INCORRECT(idir, k)                                                     \
  (isZero(vario->getHhByIndex(idir, k)) ||                                     \
   FFFF(vario->getHhByIndex(idir, k)) ||                                       \
   isZero(vario->getSwByIndex(idir, k)) ||                                     \
   FFFF(vario->getSwByIndex(idir, k)) || FFFF(vario->getGgByIndex(idir, k)))

ModelOptimSillsVario::ModelOptimSillsVario(Model* model,
                                           Constraints* constraints,
                                           const Option_AutoFit& mauto,
                                           const Option_VarioFit& optvar)
  : AModelOptimSills(model, constraints, mauto, optvar)
  , _vario()
  , _wmode(2)
{
}

ModelOptimSillsVario::ModelOptimSillsVario(const ModelOptimSillsVario& m)
  : AModelOptimSills(m)
  , _vario(m._vario)
  , _wmode(m._wmode)
{
}

ModelOptimSillsVario& ModelOptimSillsVario::operator=(const ModelOptimSillsVario& m)
{
  if (this != &m)
  {
    AModelOptimSills::operator=(m);
    _vario = m._vario;
    _wmode = m._wmode;
  }
  return (*this);
}

ModelOptimSillsVario::~ModelOptimSillsVario()
{
}

int ModelOptimSillsVario::loadEnvironment(Vario* vario, int wmode, bool verbose)
{
  _vario = vario;
  _wmode = wmode;
  _modelPart._verbose = verbose;

  // Get internal dimension
  if (_getDimensions()) return 1;

  // Allocate internal arrays
  _allocateInternalArrays(true);

  // Initialize Model-free quantities
  _wt = vario->computeWeightsFromVario(wmode);
  _compressArray(_wt, _wtc);
  _computeGg();
  _compressArray(_gg, _ggc);

  // Initialize the array of sills
  _resetSill(_ncova, _sill);

  return 0;
}

  /****************************************************************************/
/*!
 **  General Routine for fitting a model using an experimental variogram
 **
 ** \return  Error return code
 **
 ** \param[in]  vario       Experimental variogram
 ** \param[in]  wmode       Weighting mode
 ** \param[in]  verbose     Verbose flag
   **
   *****************************************************************************/
  int ModelOptimSillsVario::fit(Vario * vario, int wmode, bool verbose)
  {
    // Define the environment
    if (loadEnvironment(vario, wmode, verbose)) return 1;

    // Initialize Model-dependent quantities
    updateFromModel();

    // Perform the sill fitting
    return fitPerform();
  }

  /****************************************************************************/
  /*!
   **  Calculate the main dimensions
   **
   *****************************************************************************/
  int ModelOptimSillsVario::_getDimensions()
  {
    _ndim        = _modelPart._model->getDimensionNumber();
    _nvar        = _modelPart._model->getNVar();
    _ncova       = _modelPart._model->getCovaNumber();
    Vario* vario = _vario;

    int nbexp  = 0;
    int npadir = 0;

    // Possibly update the distance for first lag
    // if equal to 0 but corresponds to lots of pairs attached
    // This patch is not performed for asymetrical case as the h=0 is only
    // conventional.
    for (int idir = 0; idir < vario->getDirectionNumber(); idir++)
    {
      for (int ivar = 0; ivar < _nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++)
        {
          int iad0   = vario->getCenter(ivar, jvar, idir);
          double sw0 = vario->getSwByIndex(idir, iad0);
          double hh0 = vario->getHhByIndex(idir, iad0);
          // The test on the number of pairs avoids hacking in the case
          // of a conventional construction where the number of pairs
          // for the first lag is arbitrarily set to 1.
          if (isZero(hh0) && sw0 > 1.)
          {
            int iad    = vario->getNext(ivar, jvar, idir);
            double sw1 = vario->getSwByIndex(idir, iad);
            double hh1 = vario->getHhByIndex(idir, iad);

            if (!vario->getFlagAsym())
            {
              hh0 = hh1 * sw0 / sw1;
              vario->setHhByIndex(idir, iad0, hh0);
            }
          }
        }
    }

    /* Calculate the total number of lags */

    for (int idir = 0; idir < vario->getDirectionNumber(); idir++)
    {
      npadir += vario->getLagTotalNumber(idir);
      for (int ipas = 0; ipas < vario->getLagNumber(idir); ipas++)
        for (int ivar = 0; ivar < _nvar; ivar++)
          for (int jvar = 0; jvar <= ivar; jvar++)
          {
            int i = vario->getDirAddress(idir, ivar, jvar, ipas, false, 1);
            if (CORRECT(idir, i)) nbexp++;
          }
    }

    if (nbexp <= 0)
    {
      messerr("No active experimental variogram");
      return (1);
    }

    _nbexp  = nbexp;
    _npadir = npadir;
    return (0);
  }

  /****************************************************************************/
  /*!
   **  Fill the array of pointers on the experimental conditions
   **
   *****************************************************************************/
  void ModelOptimSillsVario::_computeGg()
  {
    Vario* vario = _vario;

    int ipadir = 0;
    for (int idir = 0, ndir = vario->getDirectionNumber(); idir < ndir; idir++)
    {
      for (int ipas = 0, npas = vario->getLagNumber(idir); ipas < npas; ipas++, ipadir++)
      {
        int ijvar = 0;
        for (int ivar = ijvar = 0; ivar < _nvar; ivar++)
          for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
          {

            // Calculate the variogram value
            double dist        = 0.;
            _GG(ijvar, ipadir) = TEST;
            if (vario->getFlagAsym())
            {
              int iad = vario->getDirAddress(idir, ivar, jvar, ipas, false, 1);
              int jad = vario->getDirAddress(idir, ivar, jvar, ipas, false, -1);
              double c00 = vario->getC00(idir, ivar, jvar);
              double n1  = vario->getSwByIndex(idir, iad);
              double n2  = vario->getSwByIndex(idir, jad);
              if (n1 + n2 > 0)
              {
                double g1 = vario->getGgByIndex(idir, iad);
                double g2 = vario->getGgByIndex(idir, jad);
                if (CORRECT(idir, iad) && CORRECT(idir, jad))
                {
                  _GG(ijvar, ipadir) = c00 - (n1 * g1 + n2 * g2) / (n1 + n2);
                  dist               = (ABS(vario->getHhByIndex(idir, iad)) +
                          ABS(vario->getHhByIndex(idir, jad))) / 2.;
                }
              }
            }
            else
            {
              int iad = vario->getDirAddress(idir, ivar, jvar, ipas, false, 1);
              if (CORRECT(idir, iad))
              {
                _GG(ijvar, ipadir) = vario->getGgByIndex(idir, iad);
                dist               = ABS(vario->getHhByIndex(idir, iad));
              }
            }

            // Store the distances
            int i = vario->getDirAddress(idir, ivar, jvar, ipas, false, 1);
            for (int idim = 0; idim < _ndim; idim++)
            {
              if (INCORRECT(idir, i)) continue;
              DD(idim, ijvar, ipadir) = dist * vario->getCodir(idir, idim);
            }
          }
      }
    }
  }

  /*****************************************************************************/
  /*!
   **  Calculates the values of a generic covariance model corresponding
   **  to the lags of an experimental variogram
   **
   ** \param[in]  vario   Vario structure
   ** \param[in]  model   Model structure
   ** \param[in]  npadir  Total number of lags
   **
   ** \param[out] dd      Array of distances (optional)
   ** \param[out] ge      Array of generic covariance values (optional)
   **
   *****************************************************************************/
  void ModelOptimSillsVario::updateFromModel()
  {
    Model* model = _modelPart._model;
    Vario* vario = _vario;

    int norder = 0;
    if (vario->getCalcul() == ECalcVario::GENERAL1) norder = 1;
    if (vario->getCalcul() == ECalcVario::GENERAL2) norder = 2;
    if (vario->getCalcul() == ECalcVario::GENERAL3) norder = 3;
    VectorDouble d1(_ndim);
    CovCalcMode mode = CovCalcMode(ECalcMember::LHS);
    mode.setAsVario(true);
    mode.setUnitary(true);
    mode.setOrderVario(norder);

    /* Loop on the basic structures */

    for (int icov = 0; icov < model->getCovaNumber(); icov++)
    {
      ACov* cova = model->getCova(icov);
      for (int idim = 0; idim < _ndim; idim++) d1[idim] = 0.;

      /* Loop on the experiments */

      int ipadir = 0;
      for (int idir = 0, ndir = vario->getDirectionNumber(); idir < ndir;
           idir++)
      {
        for (int ipas = 0, npas = vario->getLagNumber(idir); ipas < npas;
             ipas++, ipadir++)
        {
          int ijvar = 0;
          for (int ivar = 0; ivar < _nvar; ivar++)
            for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
            {
              int shift = ijvar * vario->getLagTotalNumber(idir);
              if (!_ge.empty()) _ge[icov].setValue(ijvar, ipadir, 0.);

              double dist = 0.;
              if (vario->getFlagAsym())
              {
                int iad = shift + vario->getLagNumber(idir) + ipas + 1;
                int jad = shift + vario->getLagNumber(idir) - ipas - 1;
                if (INCORRECT(idir, iad) || INCORRECT(idir, jad)) continue;
                dist = (ABS(vario->getHhByIndex(idir, iad)) +
                        ABS(vario->getHhByIndex(idir, jad))) / 2.;
              }
              else
              {
                int iad = shift + ipas;
                if (INCORRECT(idir, iad)) continue;
                dist = ABS(vario->getHhByIndex(idir, iad));
              }
              for (int idim = 0; idim < _ndim; idim++)
                d1[idim] = dist * vario->getCodir(idir, idim);

              if (!_ge.empty())
                _ge[icov].setValue(
                  ijvar, ipadir, cova->evalIvarIpas(1., d1, ivar, jvar, &mode));

              if (!_dd.empty())
                for (int idim = 0; idim < _ndim; idim++)
                  DD(idim, ijvar, ipadir) = d1[idim];
            }
        }
      }
    }
  }

  /****************************************************************************/
  /*!
   **  Prepare the array for Goulard's algorithm
   **  in the case of Variogram calculation
   **
   *****************************************************************************/
  void ModelOptimSillsVario::_prepareGoulard()
  {
    Model* model = _modelPart._model;
    VectorDouble tab(_nvar * _nvar);
    VectorDouble d0(_ndim);
    CovCalcMode mode(ECalcMember::LHS);
    mode.setAsVario(true);
    mode.setUnitary(true);
    // mode.setOrderVario(STRMOD->norder);

    /* Loop on the basic structures */

    for (int icov = 0, ncov = _ncova; icov < ncov; icov++)
    {
      mode.setActiveCovListFromOne(icov);

      /* Loop on the experiments */

      for (int ipadir = 0; ipadir < _npadir; ipadir++)
      {
        int ijvar = 0;
        for (int ivar = 0; ivar < _nvar; ivar++)
          for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
          {
            int flag_test = 0;
            for (int idim = 0; idim < _ndim && flag_test == 0; idim++)
            {
              d0[idim] = DD(idim, ijvar, ipadir);
              if (FFFF(d0[idim])) flag_test = 1;
            }
            if (flag_test)
            {
              _ge[icov].setValue(ijvar, ipadir, TEST);
            }
            else
            {
              _ge[icov].setValue(
                ijvar, ipadir, model->evalIvarIpas(1., d0, ivar, jvar, &mode));
            }
          }
      }
    }
  }

  /****************************************************************************/
  /*!
   **  Compress the weights for the experimental variograms
   **
   ** \param[in]  tabin     Uncompressed array
   **
   ** \param[out] tabout    Compressed array
   **
   *****************************************************************************/
  void ModelOptimSillsVario::_compressArray(const VectorDouble& tabin,
                                            VectorDouble& tabout)
  {
    Vario* vario = _vario;

    int ecr    = 0;
    int ipadir = 0;
    for (int idir = 0, ndir = vario->getDirectionNumber(); idir < ndir; idir++)
      for (int ipas = 0, npas = vario->getLagNumber(idir); ipas < npas;
           ipas++, ipadir++)
      {
        int ijvar = 0;
        for (int ivar = 0; ivar < _nvar; ivar++)
          for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
          {
            double tabval = TAB(ijvar, ipadir);
            if (!FFFF(tabval)) tabout[ecr++] = tabval;
          }
      }
  }
