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
#include "Model/ModelFitSillsVario.hpp"

#include "Variogram/Vario.hpp"
#include "Covariances/CovBase.hpp"
#include "Model/Model.hpp"
#include "Model/ModelOptimVario.hpp"
#include "Model/Constraints.hpp"

#define IJDIR(ijvar, ipadir)    ((ijvar) * _npadir + (ipadir))
#define _WT(ijvar, ipadir)      _wt[IJDIR(ijvar, ipadir)]
#define _GG(ijvar, ipadir)      _gg[IJDIR(ijvar, ipadir)]
#define _WT2(ijvar, ipadir)     _wt2[IJDIR(ijvar, ipadir)]
#define _GG2(ijvar, ipadir)     _gg2[IJDIR(ijvar, ipadir)]
#define TAB(ijvar, ipadir)      tabin[IJDIR(ijvar, ipadir)]
#define DD(idim, ijvar, ipadir) _dd[idim][IJDIR(ijvar, ipadir)]

#define CORRECT(idir, k)                    \
  (!isZero(_vario->getHhByIndex(idir, k)) && \
   !FFFF(_vario->getHhByIndex(idir, k)) &&   \
   !isZero(_vario->getSwByIndex(idir, k)) && \
   !FFFF(_vario->getSwByIndex(idir, k)) && !FFFF(_vario->getGgByIndex(idir, k)))
#define INCORRECT(idir, k)                 \
  (isZero(_vario->getHhByIndex(idir, k)) || \
   FFFF(_vario->getHhByIndex(idir, k)) ||   \
   isZero(_vario->getSwByIndex(idir, k)) || \
   FFFF(_vario->getSwByIndex(idir, k)) || FFFF(_vario->getGgByIndex(idir, k)))

ModelFitSillsVario::ModelFitSillsVario(Vario* vario,
                                       ModelCovList* model,
                                       const Constraints* constraints,
                                       const ModelOptimParam& mop)
  : AModelFitSills(model, constraints, mop)
  , _vario(vario)
{
  (void)_prepare();
}

ModelFitSillsVario::ModelFitSillsVario(const ModelFitSillsVario& m)
  : AModelFitSills(m)
  , _vario(m._vario)
{
  (void)_prepare();
}

ModelFitSillsVario& ModelFitSillsVario::operator=(const ModelFitSillsVario& m)
{
  if (this != &m)
  {
    AModelFitSills::operator=(m);
    _vario = m._vario;
    (void)_prepare();
  }
  return (*this);
}

ModelFitSillsVario::~ModelFitSillsVario()
{
}

ModelFitSillsVario* ModelFitSillsVario::createForOptim(Vario* vario,
                                                       ModelGeneric* model,
                                                       const Constraints* constraints,
                                                       const ModelOptimParam& mop)
{
  ModelCovList* mcv = dynamic_cast<ModelCovList*>(model);
  if (mcv == nullptr)
  {
    messerr("The argument 'model' should be a 'ModelCovList'");
    return nullptr;
  }
  ModelFitSillsVario* optim = new ModelFitSillsVario(vario, mcv, constraints, mop);

  return optim;
}

int ModelFitSillsVario::_prepare()
{
  // Get internal dimension
  if (_getDimensions()) return 1;

  // Allocate internal arrays
  _allocateInternalArrays(true);

  // Initialize Model-free quantities
  int wmode = _mop.getWmode();
  _wt       = _vario->computeWeightsFromVario(wmode);
  _compressArray(_wt, _wtc);
  _computeGg();
  _compressArray(_gg, _ggc);

  // Initialize the array of sills
  _resetInitialSill(_sill);

  int norder = 0;
  if (_vario->getCalcul() == ECalcVario::GENERAL1) norder = 1;
  if (_vario->getCalcul() == ECalcVario::GENERAL2) norder = 2;
  if (_vario->getCalcul() == ECalcVario::GENERAL3) norder = 3;
  _calcmode = CovCalcMode(ECalcMember::LHS);
  _calcmode.setAsVario(true);
  _calcmode.setUnitary(true);
  _calcmode.setOrderVario(norder);
  return 0;
}

/****************************************************************************/
/*!
 **  General Routine for fitting a model using an experimental variogram
 **
 ** \return  Error return code
 **
 *****************************************************************************/
int ModelFitSillsVario::fitSills()
{
  // Initialize Model-dependent quantities
  _updateFromModel();

  // In this iterative manner of Fitting Sills, the verbose flag is switched OFF
  // in order to avoid intermediate printouts
  setVerbose(false);
  int status = _fitSills();

  return status;
}

/****************************************************************************/
/*!
 **  Calculate the main dimensions
 **
 *****************************************************************************/
int ModelFitSillsVario::_getDimensions()
{
  _ndim  = _model->getNDim();
  _nvar  = _model->getNVar();
  _ncova = _model->getNCov();
  _nvs2 = _nvar * (_nvar + 1) / 2;

  int nbexp  = 0;
  int npadir = 0;

  // Possibly update the distance for first lag
  // if equal to 0 but corresponds to lots of pairs attached
  // This patch is not performed for asymetrical case as the h=0 is only
  // conventional.
  for (int idir = 0; idir < _vario->getNDir(); idir++)
  {
    for (int ivar = 0; ivar < _nvar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++)
      {
        int iad0   = _vario->getCenter(ivar, jvar, idir);
        double sw0 = _vario->getSwByIndex(idir, iad0);
        double hh0 = _vario->getHhByIndex(idir, iad0);
        // The test on the number of pairs avoids hacking in the case
        // of a conventional construction where the number of pairs
        // for the first lag is arbitrarily set to 1.
        if (isZero(hh0) && sw0 > 1.)
        {
          int iad    = _vario->getNext(ivar, jvar, idir);
          double sw1 = _vario->getSwByIndex(idir, iad);
          double hh1 = _vario->getHhByIndex(idir, iad);

          if (!_vario->getFlagAsym())
          {
            hh0 = hh1 * sw0 / sw1;
            _vario->setHhByIndex(idir, iad0, hh0);
          }
        }
      }
  }

  /* Calculate the total number of lags */

  for (int idir = 0; idir < _vario->getNDir(); idir++)
  {
    npadir += _vario->getNLagTotal(idir);
    for (int ilag = 0; ilag < _vario->getNLag(idir); ilag++)
      for (int ivar = 0; ivar < _nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++)
        {
          int i = _vario->getDirAddress(idir, ivar, jvar, ilag, false, 1);
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
void ModelFitSillsVario::_computeGg()
{
  int ipadir = 0;
  for (int idir = 0, ndir = _vario->getNDir(); idir < ndir; idir++)
  {
    for (int ilag = 0, nlag = _vario->getNLag(idir); ilag < nlag; ilag++, ipadir++)
    {
      int ijvar = 0;
      for (int ivar = ijvar = 0; ivar < _nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
        {

          // Calculate the variogram value
          double dist        = 0.;
          _GG(ijvar, ipadir) = TEST;
          if (_vario->getFlagAsym())
          {
            int iad    = _vario->getDirAddress(idir, ivar, jvar, ilag, false, 1);
            int jad    = _vario->getDirAddress(idir, ivar, jvar, ilag, false, -1);
            double c00 = _vario->getC00(idir, ivar, jvar);
            double n1  = _vario->getSwByIndex(idir, iad);
            double n2  = _vario->getSwByIndex(idir, jad);
            if (n1 + n2 > 0)
            {
              double g1 = _vario->getGgByIndex(idir, iad);
              double g2 = _vario->getGgByIndex(idir, jad);
              if (CORRECT(idir, iad) && CORRECT(idir, jad))
              {
                _GG(ijvar, ipadir) = c00 - (n1 * g1 + n2 * g2) / (n1 + n2);
                dist               = (ABS(_vario->getHhByIndex(idir, iad)) +
                        ABS(_vario->getHhByIndex(idir, jad))) / 2.;
              }
            }
          }
          else
          {
            int iad = _vario->getDirAddress(idir, ivar, jvar, ilag, false, 1);
            if (CORRECT(idir, iad))
            {
              _GG(ijvar, ipadir) = _vario->getGgByIndex(idir, iad);
              dist               = ABS(_vario->getHhByIndex(idir, iad));
            }
          }

          // Store the distances
          int i = _vario->getDirAddress(idir, ivar, jvar, ilag, false, 1);
          for (int idim = 0; idim < _ndim; idim++)
          {
            if (INCORRECT(idir, i)) continue;
            DD(idim, ijvar, ipadir) = dist * _vario->getCodir(idir, idim);
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
 *****************************************************************************/
void ModelFitSillsVario::_updateFromModel()
{
  VectorDouble d1(_ndim);

  /* Loop on the basic structures */

  for (int icov = 0; icov < _model->getNCov(); icov++)
  {
    const CovBase* cova = _model->getCovList()->getCov(icov);
    d1.fill(0.);

    /* Loop on the experiments */

    int ipadir = 0;
    for (int idir = 0, ndir = _vario->getNDir(); idir < ndir; idir++)
    {
      for (int ilag = 0, nlag = _vario->getNLag(idir); ilag < nlag; ilag++, ipadir++)
      {
        int ijvar = 0;
        for (int ivar = 0; ivar < _nvar; ivar++)
          for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
          {
            int shift = ijvar * _vario->getNLagTotal(idir);
            if (!_ge.empty()) _ge[icov].setValue(ijvar, ipadir, 0.);

            double dist = 0.;
            if (_vario->getFlagAsym())
            {
              int iad = shift + _vario->getNLag(idir) + ilag + 1;
              int jad = shift + _vario->getNLag(idir) - ilag - 1;
              if (INCORRECT(idir, iad) || INCORRECT(idir, jad)) continue;
              dist = (ABS(_vario->getHhByIndex(idir, iad)) +
                      ABS(_vario->getHhByIndex(idir, jad))) / 2.;
            }
            else
            {
              int iad = shift + ilag;
              if (INCORRECT(idir, iad)) continue;
              dist = ABS(_vario->getHhByIndex(idir, iad));
            }
            for (int idim = 0; idim < _ndim; idim++)
              d1[idim] = dist * _vario->getCodir(idir, idim);

            if (!_ge.empty())
              _ge[icov].setValue(
                ijvar, ipadir, cova->evalIvarIpas(1., d1, ivar, jvar, &_calcmode));

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
void ModelFitSillsVario::_prepareGoulard()
{
  VectorDouble tab(_nvar * _nvar);
  VectorDouble d0(_ndim);
  CovCalcMode mode(ECalcMember::RHS);
  mode.setAsVario(true);
  mode.setUnitary(true);

  /* Loop on the basic structures */

  const CovList* cova = _model->getCovList();
  for (int icov = 0, ncov = _ncova; icov < ncov; icov++)
  {
    cova->setActiveCovListFromOne(icov);

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
              ijvar, ipadir, _model->evalIvarIpas(1., d0, ivar, jvar, &mode));
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
void ModelFitSillsVario::_compressArray(const VectorDouble& tabin,
                                        VectorDouble& tabout)
{
  int ecr    = 0;
  int ipadir = 0;
  for (int idir = 0, ndir = _vario->getNDir(); idir < ndir; idir++)
    for (int ilag = 0, nlag = _vario->getNLag(idir); ilag < nlag;
         ilag++, ipadir++)
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
