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
#include "Model/ModelOptimSills.hpp"

#include "Model/Model.hpp"
#include "Variogram/Vario.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"
#include "Model/Constraints.hpp"
#include "Basic/MathFunc.hpp"

#define IJDIR(ijvar, ipadir) ((ijvar) * _npadir + (ipadir))
#define WT(ijvar, ipadir)       wt[IJDIR(ijvar, ipadir)]
#define GG(ijvar, ipadir)       gg[IJDIR(ijvar, ipadir)]
#define _WT(ijvar, ipadir)       _wt[IJDIR(ijvar, ipadir)]
#define _GG(ijvar, ipadir)       _gg[IJDIR(ijvar, ipadir)]
#define _WT2(ijvar, ipadir) _wt2[IJDIR(ijvar, ipadir)]
#define _GG2(ijvar, ipadir)      _gg2[IJDIR(ijvar, ipadir)]
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

ModelOptimSills::ModelOptimSills()
  : ModelOptimVario()
{
}

ModelOptimSills::ModelOptimSills(const ModelOptimSills& m)
  : ModelOptimVario(m)
{
}

ModelOptimSills& ModelOptimSills::operator=(const ModelOptimSills& m)
{
  if (this != &m)
  {
    ModelOptimVario::operator=(m);
  }
  return (*this);
}

ModelOptimSills::~ModelOptimSills()
{
}

/****************************************************************************/
/*!
 **  General Routine for fitting a model using an experimental variogram
 **
 ** \return  Error return code
 **
 ** \param[in]  vario       Experimental variogram
 ** \param[in]  model       Model to be fitted
 ** \param[in]  constraints Constraints structure
 ** \param[in]  mauto       Option_AutoFit structure
 ** \param[in]  optvar      Option_VarioFit structure
 **
 *****************************************************************************/
int ModelOptimSills::fit(Vario* vario,
                         Model* model,
                         Constraints* constraints,
                         const Option_AutoFit* mauto,
                         const Option_VarioFit* optvar)
{
  // Assign to local variables
  _modelPart._model   = model;
  _varioPart._vario   = vario;
  _constraints = constraints;
  _mauto       = mauto;
  _optvar      = optvar;

  // Get internal dimension
  if (_getDimensions()) return 1;

  // Allocate internal arrays
  _allocateInternalArrays();

  // Load the arrays
  _computeWt();
  _compressArray(_wt, _wtc);
  _computeGg();
  _compressArray(_gg, _ggc);
  _computeGe();

  // Initialize the array of sills
  _resetSill(_ncova, _sill);

  /* Dispatch */

  int status = 0;
  double crit = 0.;
  if (!optvar->getFlagIntrinsic())
  {

    /* No intrinsic hypothesis */

    if (FFFF(constraints->getConstantSillValue()))
    {
      /* Without constraint on the sill */

      status = _goulardWithoutConstraint(_mauto, _nvar, _ncova, _npadir,
                                         _wt, _gg, _ge, _sill, &crit);
    }
    else
    {

      /* With constraint on the sill */

      status = _goulardWithConstraints();
    }

    // Store the sills into the Model
    _storeSillsInModel();
  }
  else
  {
    status = _sillFittingIntrinsic();
  }

  return (status);
}

void ModelOptimSills::_storeSillsInModel() const
{
  for (int icov = 0; icov < _ncova; icov++)
  {
    int ijvar = 0;
    for (int ivar = ijvar = 0; ivar < _nvar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
        _setSill(icov, ivar, jvar, _sill[icov].getValue(ivar, jvar));
  }
}

void ModelOptimSills::_resetSill(int ncova, std::vector<MatrixSquareSymmetric>& sill) const
{
  for (int icova = 0; icova < ncova; icova++)
  {
    for (int ivar = 0; ivar < _nvar; ivar++)
      for (int jvar = 0; jvar < _nvar; jvar++)
        sill[icova].setValue(ivar, jvar, ivar == jvar);
  }
}

/****************************************************************************/
/*!
 **  Calculate the main dimensions
 **
 *****************************************************************************/
int ModelOptimSills::_getDimensions()
{
  _ndim  = _modelPart._model->getDimensionNumber();
  _nvar  = _modelPart._model->getVariableNumber();
  _ncova = _modelPart._model->getCovaNumber();
  Vario* vario = _varioPart._vario;

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

  _nbexp = nbexp;
  _npadir = npadir;
  return (0);
}

/****************************************************************************/
/*!
 **  Manage memory for variogram fitting
 **
 ** \return  Error returned code
 **
 ** \param[in]  flag_exp  1 for experimental variogram
 **
 *****************************************************************************/
void ModelOptimSills::_allocateInternalArrays(bool flag_exp)
{
  int nvs2 = _nvar * (_nvar + 1) / 2;

  _wt.fill(TEST, _npadir * nvs2);
  _gg.fill(TEST, _npadir * nvs2);
  _ge.clear();
  for (int icova = 0; icova < _ncova; icova++)
    _ge.push_back(MatrixRectangular(nvs2, _npadir));
  _sill.clear();
  for (int icova = 0; icova < _ncova; icova++)
    _sill.push_back(MatrixSquareSymmetric(_nvar));

  if (flag_exp)
  {
    _wtc.fill(TEST, _nbexp);
    _ggc.fill(TEST, _nbexp);
    _dd.clear();
    for (int idim = 0; idim < _ndim; idim++)
      _dd.push_back(VectorDouble(_npadir * nvs2, TEST));
  }

  if (_optvar->getFlagIntrinsic())
  {
    _alphau.clear();
    for (int icova = 0; icova < _ncova; icova++)
      _alphau.push_back(MatrixSquareSymmetric(1));
    _ge1.clear();
    _ge1.push_back(MatrixRectangular(nvs2, _npadir));
    _ge2.clear();
    for (int icova = 0; icova < _ncova; icova++)
      _ge2.push_back(MatrixRectangular(nvs2, _npadir));
    _wt2.fill(TEST, nvs2 * _npadir);
    _gg2.fill(TEST, nvs2 * _npadir);
  }
}

/****************************************************************************/
/*!
 **  Fill the array of pointers on the experimental conditions
 **
 *****************************************************************************/
void ModelOptimSills::_computeGg()
{
  Vario* vario = _varioPart._vario;

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
          double dist       = 0.;
          _GG(ijvar, ipadir) = TEST;
          if (vario->getFlagAsym())
          {
            int iad = vario->getDirAddress(idir, ivar, jvar, ipas, false, 1);
            int jad = vario->getDirAddress(idir, ivar, jvar, ipas, false, -1);
            double c00 = _getC00(idir, ivar, jvar);
            double n1  = vario->getSwByIndex(idir, iad);
            double n2  = vario->getSwByIndex(idir, jad);
            if (n1 + n2 > 0)
            {
              double g1 = vario->getGgByIndex(idir, iad);
              double g2 = vario->getGgByIndex(idir, jad);
              if (CORRECT(idir, iad) && CORRECT(idir, jad))
              {
                _GG(ijvar, ipadir) = c00 - (n1 * g1 + n2 * g2) / (n1 + n2);
                dist              = (ABS(vario->getHhByIndex(idir, iad)) +
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
              dist              = ABS(vario->getHhByIndex(idir, iad));
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
void ModelOptimSills::_computeGe()
{
  Vario* vario = _varioPart._vario;
  Model* model = _modelPart._model;

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
    for (int idir = 0, ndir = vario->getDirectionNumber(); idir < ndir; idir++)
    {
      for (int ipas = 0, npas = vario->getLagNumber(idir); ipas < npas; ipas++, ipadir++)
      {
        int ijvar = 0;
        for (int ivar = 0; ivar < _nvar; ivar++)
          for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
          {
            int shift = ijvar * vario->getLagTotalNumber(idir);
            if (! _ge.empty()) _ge[icov].setValue(ijvar, ipadir, 0.);

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
              _ge[icov].setValue(ijvar, ipadir,
                                cova->evalIvarIpas(1., d1, ivar, jvar, &mode));

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
 **  Compress the weights for the experimental variograms
 **
 ** \param[in]  tabin     Uncompressed array
 **
 ** \param[out] tabout    Compressed array
 **
 *****************************************************************************/
void ModelOptimSills::_compressArray(const VectorDouble& tabin, VectorDouble& tabout)
{
  Vario* vario = _varioPart._vario;

  int ecr    = 0;
  int ipadir = 0;
  for (int idir = 0, ndir = vario->getDirectionNumber(); idir < ndir; idir++)
    for (int ipas = 0, npas = vario->getLagNumber(idir); ipas < npas; ipas++, ipadir++)
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

/****************************************************************************/
/*!
 **  Prepare the array for Goulard's algorithm
 **  in the case of Variogram calculation
 **
 *****************************************************************************/
void ModelOptimSills::_preparGoulardVario()
{
  Model* model = _modelPart._model;
  VectorDouble tab(_nvar * _nvar);
  VectorDouble d0(_ndim);
  CovCalcMode mode(ECalcMember::LHS);
  mode.setAsVario(true);
  mode.setUnitary(true);
  //mode.setOrderVario(STRMOD->norder);

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
            _ge[icov].setValue(ijvar, ipadir,
                              model->evalIvarIpas(1., d0, ivar, jvar, &mode));
          }
        }
    }
  }
}

int ModelOptimSills::_goulardWithConstraints()
{
  double crit;
  VectorDouble consSill = _constraints->getConstantSills();

  /* Core allocation */
  std::vector<MatrixSquareSymmetric> matcor;
  matcor.reserve(_ncova);
  for (int icova = 0; icova < _ncova; icova++)
    matcor.push_back((MatrixSquareSymmetric(_nvar)));

  /* Initialize the Goulard system */

  _initializeGoulard();

  /* Update according to the eigen values */

  bool flag_positive = true;
  for (int icov = 0; icov < _ncova; icov++)
    if (!_makeDefinitePositive(icov))
      flag_positive = false;

  if (!flag_positive)
  {
    for (int icov = 0; icov < _ncova; icov++)
      for (int ivar = 0; ivar < _nvar; ivar++)
      {
        _sill[icov].fill(0.);
        if (!FFFF(consSill[ivar]))
          _sill[icov].setValue(ivar, ivar, consSill[ivar] / _ncova);
        else
          _sill[icov].setValue(ivar, ivar, 1.);
      }

    /* Perform the optimization under constraints */

    (void)_optimizeUnderConstraints(&crit);
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Initialize the system for Goulard algorithm
 **
 *****************************************************************************/
void ModelOptimSills::_initializeGoulard()
{
  MatrixSquareSymmetric aa(_ncova);
  VectorDouble bb(_ncova);
  MatrixRectangular Ae(_ncova, 1);
  VectorDouble be(1);
  MatrixRectangular Ai(_ncova, _ncova);
  VectorDouble bi(_ncova);
  VectorDouble res(_ncova);
  VectorDouble consSill = _constraints->getConstantSills();

  /* Initialize the constraints matrices */

  Ae.fill(1.);
  bi.fill(0.);
  for (int icov = 0; icov < _ncova; icov++)
    for (int jcov = 0; jcov < _ncova; jcov++)
      Ai.setValue(icov, jcov, (icov == jcov));

  /* Loop on the variables */

  int ijvar = 0;
  for (int ivar = 0; ivar < _nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
    {
      /* Reset the arrays */

      bb.fill(0.);
      aa.fill(0.);

      for (int ipadir = 0; ipadir < _npadir; ipadir++)
      {
        double wtloc = _WT(ijvar, ipadir);
        double ggloc = _GG(ijvar, ipadir);
        if (FFFF(wtloc) || FFFF(ggloc)) continue;
        for (int icov = 0; icov < _ncova; icov++)
        {
          double geloc1 = _ge[icov].getValue(ijvar, ipadir);
          if (FFFF(geloc1)) continue;
          bb[icov] += wtloc * ggloc * geloc1;
          for (int jcov = 0; jcov <= icov; jcov++)
          {
            double geloc2 = _ge[jcov].getValue(ijvar, ipadir);
            if (FFFF(geloc2)) continue;
            aa.updValue(icov, jcov, EOperator::ADD, wtloc * geloc1 * geloc2);
          }
        }
      }

      int retcode = 0;
      if (ivar == jvar && !consSill.empty() && !FFFF(consSill[ivar]))
      {
        be[0]   = consSill[ivar];
        retcode = aa.minimizeWithConstraintsInPlace(bb, Ae, be, Ai, bi, res);
      }
      else
      {
        retcode = aa.minimizeWithConstraintsInPlace(
          bb, MatrixRectangular(), VectorDouble(), MatrixRectangular(),
          VectorDouble(), res);
      }

      /* Update (taking into account possible constraints) */

      if (retcode)
      {
        for (int icov = 0; icov < _ncova; icov++) res[icov] = (ivar == jvar);
        if (ivar == jvar && !consSill.empty() && !FFFF(consSill[ivar]))
        {
          double temp = consSill[ivar] / _ncova;
          for (int icov = 0; icov < _ncova; icov++) res[icov] = temp;
        }
      }

      /* Store in the output matrix */

      for (int icov = 0; icov < _ncova; icov++)
        _sill[icov].setValue(ivar, jvar, res[icov]);
    }
}

/*****************************************************************************/
/*!
 ** Make sure the current matrix of sills is definite positive
 ** (diagonal unchanged)
 **
 ** \param[in]      icov0    Index of the target basic structure
 ** \param[in]      eps      Tolerance
 **
 *****************************************************************************/
int ModelOptimSills::_makeDefinitePositive(int icov0, double eps)
{
  VectorDouble muold(_nvar);
  VectorDouble norme1(_nvar);
  VectorDouble consSill = _constraints->getConstantSills();

  for (int ivar = 0; ivar < _nvar; ivar++)
    muold[ivar] = _sill[icov0].getValue(ivar, ivar);

  int flag_positive = _truncateNegativeEigen(icov0);

  if (flag_positive) return flag_positive;

  for (int ivar = 0; ivar < _nvar; ivar++)
  {
    if (FFFF(consSill[ivar]))
      norme1[ivar] = 1.;
    else
    {
      if (ABS(_sill[icov0].getValue(ivar, ivar)) > eps)
        norme1[ivar] = sqrt(muold[ivar] / _sill[icov0].getValue(ivar, ivar));
      else
        norme1[ivar] = (ABS(muold[ivar]) < eps) ? 1. : 0.;
    }
  }

  for (int ivar = 0; ivar < _nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++)
      _sill[icov0].setValue(ivar, jvar,
                            _sill[icov0].getValue(ivar, jvar) * norme1[ivar] *
                              norme1[jvar]);

  return flag_positive;
}

int ModelOptimSills::_optimizeUnderConstraints(double* score)
{
  double score_old, xrmax;

  /* Core allocation */

  VectorDouble xr(_nvar);
  std::vector<MatrixSquareSymmetric> alpha;
  VectorDouble consSill = _constraints->getConstantSills();
  alpha.reserve(_ncova);
  for (int icova = 0; icova < _ncova; icova++)
    alpha.push_back(MatrixSquareSymmetric(_nvar));
  int iter = 0;

  /* Calculate the initial score */

  double score_new = _score();

  for (int icov = 0; icov < _ncova; icov++) alpha[icov] = _sill[icov];

  /***********************/
  /* Iterative procedure */
  /***********************/

  /* Optional printout */

  for (int ivar = 0; ivar < _nvar; ivar++)
  {
    if (FFFF(consSill[ivar]))
      xr[ivar] = 1.;
    else
      xr[ivar] = sqrt(consSill[ivar] / _sumSills(ivar, alpha));
  }

  for (iter = 0; iter < _mauto->getMaxiter(); iter++)
  {
    for (int icov0 = 0; icov0 < _ncova; icov0++)
    {
      for (int ivar0 = 0; ivar0 < _nvar; ivar0++)
      {
        double denom = _sumSills(ivar0, alpha) - alpha[icov0].getValue(ivar0, ivar0);
        if (!FFFF(consSill[ivar0]) && denom < 1e-30) continue;
        if (!FFFF(consSill[ivar0]))
        {
          xrmax     = sqrt(consSill[ivar0] / denom);
          xr[ivar0] = _minimizeP4(icov0, ivar0, xrmax, xr, alpha);

          if (isZero(xr[ivar0]))
          {
            xr[ivar0] = 1.;
            for (int jcov = 0; jcov < _ncova; jcov++)
            {
              if (jcov == icov0) continue;
              for (int jvar = 0; jvar < _nvar; jvar++)
                alpha[jcov].setValue(ivar0, jvar, 0.);
            }
          }
        }

        /* Update 'alpha' (diagonal only) */

        if (!FFFF(consSill[ivar0]))
          _updateAlphaDiag(icov0, ivar0, xr, alpha);

        /* Update 'sills' for the structures other than the current one */

        _updateOtherSills(icov0, ivar0, alpha, xr);

        /* Update the sill matrix for the current structure */
        /* (except diagonal in the constrained case)        */

        _updateCurrentSillGoulard(icov0, ivar0);

        /* Update sill matrix for the current structure (for diagonal) */

        if (!FFFF(consSill[ivar0]))
          _updateCurrentSillDiag(icov0, ivar0, alpha, xr);

        /* Make sure the current matrix of sills if definite positive */
        /* (diagonal unchanged)                                       */

        (void)_makeDefinitePositive(icov0);

        /* Update 'alpha' for the current structure */

        for (int ivar = 0; ivar < _nvar; ivar++)
          _updateAlphaNoDiag(icov0, ivar, xr, alpha);
      }
    }

    /* Update the score */

    score_old = score_new;
    score_new = _score();
    if (ABS(score_new - score_old) / score_old < _mauto->getTolred()) break;
  }

  /* Optional printout */

  *score = score_new;
  return (iter);
}

int ModelOptimSills::_truncateNegativeEigen(int icov0)
{
  MatrixSquareSymmetric cc(_nvar);
  for (int ivar = 0; ivar < _nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++)
      cc.setValue(ivar, jvar, _sill[icov0].getValue(ivar, jvar));

  if (cc.computeEigen()) messageAbort("st_truncate_negative_eigen");

  VectorDouble valpro               = cc.getEigenValues();
  const MatrixSquareGeneral* vecpro = cc.getEigenVectors();

  /* Check positiveness */

  int flag_positive = 1;
  for (int ivar = 0; ivar < _nvar; ivar++)
    if (valpro[ivar] <= 0) flag_positive = 0;
  if (flag_positive) return flag_positive;

  /* Calculate the new coregionalization matrix */

  for (int ivar = 0; ivar < _nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++)
    {
      double sum = 0.;
      for (int kvar = 0; kvar < _nvar; kvar++)
        sum += MAX(valpro[kvar], 0.) * vecpro->getValue(ivar, kvar) *
               vecpro->getValue(jvar, kvar);
      _sill[icov0].setValue(ivar, jvar, sum);
    }

  return flag_positive;
}

double ModelOptimSills::_score()
{
  double score = 0.;
  int ijvar    = 0;
  for (int ivar = 0; ivar < _nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
    {
      double coeff = (ivar == jvar) ? 1 : 2;
      for (int ipadir = 0; ipadir < _npadir; ipadir++)
      {
        double dd = _GG(ijvar, ipadir);
        if (FFFF(dd)) continue;
        for (int icov = 0; icov < _ncova; icov++)
          dd -= _sill[icov].getValue(ivar, jvar) *
                _ge[icov].getValue(ijvar, ipadir);
        score += coeff * _WT(ijvar, ipadir) * dd * dd;
      }
    }
  return (score);
}

double ModelOptimSills::_sumSills(int ivar0,
                                  std::vector<MatrixSquareSymmetric>& alpha) const
{
  double Sr = 0;
  for (int icov = 0; icov < _ncova; icov++)
    Sr += alpha[icov].getValue(ivar0, ivar0);
  return Sr;
}

/****************************************************************************/
/*!
 ** Minimization of the order-4 polynomial
 **
 ** \return Value for which the fourth order polynomial is minimum.
 **
 ** \param[in] icov0    Index of the target basic structure
 ** \param[in] ivar0    Index of the target variable
 ** \param[in] xrmax    Maximum value for the solution
 ** \param[in] xr       Current vector of sqrt(constraint/(sum of the sills))
 ** \param[in] alpha    Current auxiliary matrices alpha
 **
 *****************************************************************************/
double ModelOptimSills::_minimizeP4(int icov0,
                                    int ivar0,
                                    double xrmax,
                                    VectorDouble& xr,
                                    std::vector<MatrixSquareSymmetric>& alpha)
{
  double retval, a, c, d, s;

  /* Core allocation */

  VectorDouble Nir_v(_nvar);
  VectorDouble Mrr_v(_npadir);
  VectorDouble Crr_v(_npadir);
  MatrixRectangular Airk_v(_npadir, _nvar);
  MatrixRectangular Birk_v(_npadir, _nvar);
  VectorDouble xx(2);
  VectorDouble xt(2);
  VectorDouble xest(2);
  VectorDouble x(3);
  VectorDouble consSill = _constraints->getConstantSills();

  int irr = _combineVariables(ivar0, ivar0);

  for (int ivar = 0; ivar < _nvar; ivar++)
  {
    int irl     = _combineVariables(ivar0, ivar);
    Nir_v[ivar] = 0.;
    for (int k = 0; k < _npadir; k++)
      Nir_v[ivar] +=
        _WT(irl, k) * _ge[icov0].getValue(0, k) * _ge[icov0].getValue(0, k);
  }

  for (int k = 0; k < _npadir; k++)
  {
    Mrr_v[k] = 0.;
    for (int jcov = 0; jcov < _ncova; jcov++)
    {
      if (jcov == icov0) continue;
      Mrr_v[k] += alpha[jcov].getValue(ivar0, ivar0) *
                  (_ge[jcov].getValue(0, k) - _ge[icov0].getValue(0, k));
    }
  }

  for (int k = 0; k < _npadir; k++)
    for (int ivar = 0; ivar < _nvar; ivar++)
    {
      int irl      = _combineVariables(ivar0, ivar);
      double value = 0.;
      for (int l = 0; l < _npadir; l++)
        value += _WT(irl, l) * _GG(irl, l) * _ge[icov0].getValue(0, l);
      value *= _ge[icov0].getValue(0, k) / Nir_v[ivar];
      Airk_v.setValue(k, ivar, _GG(irl, k) - value);
    }

  for (int k = 0; k < _npadir; k++)
    for (int ivar = 0; ivar < _nvar; ivar++)
    {
      int irl = _combineVariables(ivar0, ivar);
      Birk_v.setValue(k, ivar, 0.);
      for (int icov = 0; icov < _ncova; icov++)
      {
        if (icov == icov0) continue;
        double value = 0.;
        for (int l = 0; l < _npadir; l++)
          value += _WT(irl, l) * _ge[icov].getValue(0, l) * _ge[icov].getValue(0, l);
        Birk_v.setValue(
          k, ivar,
          Birk_v.getValue(k, ivar) +
            alpha[icov].getValue(ivar0, ivar) *
              (_ge[icov].getValue(0, k) - value * _ge[icov0].getValue(0, k)) /
              Nir_v[ivar]);
      }
    }

  for (int k = 0; k < _npadir; k++)
    Crr_v[k] = _GG(irr, k) - consSill[ivar0] * _ge[icov0].getValue(0, k);

  a = 0.;
  for (int k = 0; k < _npadir; k++) a += _WT(irr, k) * Mrr_v[k] * Mrr_v[k];

  c = 0.;
  for (int k = 0; k < _npadir; k++)
    for (int ivar = 0; ivar < _nvar; ivar++)
    {
      if (ivar != ivar0)
      {
        int irl = _combineVariables(ivar0, ivar);
        s       = xr[ivar] * Birk_v.getValue(k, ivar);
        c += _WT(irl, k) * s * s;
      }
      else
      {
        c -= _WT(irr, k) * Mrr_v[k] * Crr_v[k];
      }
    }
  c *= 2.;

  d = 0.;
  for (int k = 0; k < _npadir; k++)
    for (int ivar = 0; ivar < _nvar; ivar++)
    {
      if (ivar == ivar0) continue;
      int irl = _combineVariables(ivar0, ivar);
      d -= _WT(irl, k) * Airk_v.getValue(k, ivar) * xr[ivar] *
           Birk_v.getValue(k, ivar);
    }

  d *= 4.;

  int number = solve_P3(4. * a, 0., 2. * c, d, x);

  switch (number)
  {
    case 1: retval = MAX(0., MIN(xrmax, x[0])); break;

    case 3:
    {
      int nin = 0;
      xx[0]   = x[0];
      xx[1]   = x[2];
      for (int k = 0; k < 2; k++)
      {
        xt[k] = MAX(0., MIN(xrmax, xx[k]));
        xest[k] =
          (a * xt[k] * xt[k] * xt[k] * xt[k] + c * xt[k] * xt[k] + d * xt[k]) /
          2.;
        if (isEqual(xt[k], xx[k])) nin++;
      }
      if (nin == 1)
      {
        retval = (isEqual(xt[0], xx[0])) ? xx[0] : xx[1];
      }
      else
      {
        retval = (xest[0] < xest[1]) ? xt[0] : xt[1];
      }
    }
    break;

    default:
      retval = xr[ivar0];
      // mes_abort("st_minimize_P4: Number of extrema is impossible");
      break;
  }

  if (retval < 1.e-9) retval = 0.;

  return retval;
}

void ModelOptimSills::_updateAlphaDiag(int icov0,
                                       int ivar0,
                                       VectorDouble& xr,
                                       std::vector<MatrixSquareSymmetric>& alpha)
{
  VectorDouble consSill = _constraints->getConstantSills();
  double srm   = _sumSills(ivar0, alpha) - alpha[icov0].getValue(ivar0, ivar0);
  double value = consSill[ivar0] / (xr[ivar0] * xr[ivar0]) - srm;
  alpha[icov0].setValue(ivar0, ivar0, MAX(0., value));
}

void ModelOptimSills::_updateOtherSills(int icov0,
                                        int ivar0,
                                        std::vector<MatrixSquareSymmetric>& alpha,
                                        VectorDouble& xr)
{
  for (int jcov = 0; jcov < _ncova; jcov++)
  {
    if (jcov == icov0) continue;
    for (int ivar = 0; ivar < _nvar; ivar++)
    {
      double value = alpha[jcov].getValue(ivar0, ivar) * xr[ivar0] * xr[ivar];
      _sill[jcov].setValue(ivar0, ivar, value);
      _sill[jcov].setValue(ivar, ivar0, value);
    }
  }
}

void ModelOptimSills::_updateCurrentSillGoulard(int icov0, int ivar0)
{
  VectorDouble mv(_npadir);
  VectorDouble consSill = _constraints->getConstantSills();

  // Loop on the variables

  for (int ivar = 0; ivar < _nvar; ivar++)
  {
    if (ivar == ivar0 && !FFFF(consSill[ivar0])) continue;

    for (int ilagdir = 0; ilagdir < _npadir; ilagdir++)
    {
      mv[ilagdir] = 0.;
      for (int icov = 0; icov < _ncova; icov++)
      {
        if (icov == icov0) continue;
        mv[ilagdir] +=
          _sill[icov].getValue(ivar0, ivar) * _ge[icov].getValue(0, ilagdir);
      }
    }

    double tot1 = 0.;
    double tot2 = 0.;
    int ivs2    = _combineVariables(ivar0, ivar);

    for (int ilagdir = 0; ilagdir < _npadir; ilagdir++)
    {
      double wtloc = _WT(ivs2, ilagdir);
      double ggloc = _GG(ivs2, ilagdir);
      double geloc = _ge[icov0].getValue(0, ilagdir);

      if (!FFFF(ggloc))
      {
        tot1 += wtloc * geloc * (ggloc - mv[ilagdir]);
        tot2 += wtloc * geloc * geloc;
      }
    }
    double value = tot1 / tot2;
    _sill[icov0].setValue(ivar0, ivar, value);
    _sill[icov0].setValue(ivar, ivar0, value);
  }
}

void ModelOptimSills::_updateCurrentSillDiag(int icov0,
                                             int ivar0,
                                             std::vector<MatrixSquareSymmetric>& alpha,
                                             VectorDouble& xr)
{
  double value = xr[ivar0] * xr[ivar0] * alpha[icov0].getValue(ivar0, ivar0);
  if (value < 0.) value = 0.;
  _sill[icov0].setValue(ivar0, ivar0, value);
}

void ModelOptimSills::_updateAlphaNoDiag(int icov0,
                                         int ivar0,
                                         VectorDouble& xr,
                                         std::vector<MatrixSquareSymmetric>& alpha)
{
  VectorDouble consSill = _constraints->getConstantSills();
  for (int ivar = 0; ivar < _nvar; ivar++)
  {
    if (ivar == ivar0 && !FFFF(consSill[ivar0])) continue;
    double value = _sill[icov0].getValue(ivar0, ivar) / (xr[ivar0] * xr[ivar]);
    alpha[icov0].setValue(ivar0, ivar, value);
  }
}

int ModelOptimSills::_combineVariables(int ivar0, int jvar0)
{
  if (ivar0 > jvar0) return (ivar0 + jvar0 * (jvar0 + 1) / 2);
  return (jvar0 + ivar0 * (ivar0 + 1) / 2);
}

int ModelOptimSills::_sillFittingIntrinsic()
{
  int nvs2 = _nvar * (_nvar + 1) / 2;
  double crit = 0.;
  double crit_mem = 1.e30;
  for (int icov = 0; icov < _ncova; icov++) _alphau[icov] = 1. / (double)_ncova;

  /* Iterative procedure */

  Option_AutoFit mauto_new(*_mauto);
  mauto_new.setMaxiter(1);
  for (int iter = 0; iter < _mauto->getMaxiter(); iter++)
  {

    /* Initialize the arrays for the first pass */

    for (int ipadir = 0; ipadir < _npadir; ipadir++)
    {
      int ijvar = 0;
      for (int ivar = 0; ivar < _nvar; ivar++)
      {
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
        {
          double sum = 0.;
          for (int icov = 0; icov < _ncova; icov++)
            sum += _alphau[icov].getValue(0, 0) * _ge[icov].getValue(ijvar, ipadir);
          _ge1[0].setValue(ijvar, ipadir, sum);
        }
      }
    }

    /* Call Goulard with 1 structure (no constraint) */

    _resetSill(1, _sill1);
    if (_goulardWithoutConstraint(&mauto_new, _nvar, 1, _npadir, _wt, _gg, _ge1,
                                  _sill1, &crit)) return 1;

    /* Initialize the arrays for the second pass */

    int ijvar = 0;
    for (int ivar = 0; ivar < _nvar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
      {
        double pivot = _sill1[0].getValue(ivar, jvar);
        for (int ipadir = 0; ipadir < _npadir; ipadir++)
        {
          _GG2(ijvar, ipadir) = (isZero(pivot)) ? 0. : _GG(ijvar, ipadir) / pivot;
          _WT2(ijvar, ipadir) = _WT(ijvar, ipadir) * pivot * pivot;
          for (int icov = 0; icov < _ncova; icov++)
            _ge2[icov].setValue(ijvar, ipadir, _ge[icov].getValue(ijvar, ipadir));
        }
      }

    /* Call Goulard with 1 variable (no constraint) */

    if (_goulardWithoutConstraint(&mauto_new, 1, _ncova, _npadir * nvs2, _wt2,
                                  _gg2, _ge2, _alphau, &crit)) return 1;

    /* Stopping criterion */

    if (ABS(crit) < _mauto->getTolred() ||
        ABS(crit - crit_mem) / ABS(crit) < _mauto->getTolred())
      break;
    crit_mem = crit;
  }

  /* Patch the final model */

  for (int icov = 0; icov < _ncova; icov++)
  {
    int ijvar = 0;
    for (int ivar = 0; ivar < _nvar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
      {
        double newval = _alphau[icov].getValue(0, 0) * _sill1[0].getValue(ivar, jvar);
        _setSill(icov, ivar, jvar, newval);
        _setSill(icov, jvar, ivar, newval);
      }
  }
  return 0;
}

int ModelOptimSills::_goulardWithoutConstraint(
  const Option_AutoFit* mauto,
  int nvar,
  int ncova,
  int npadir,
  VectorDouble& wt,
  VectorDouble& gg,
  std::vector<MatrixRectangular>& ge,
  std::vector<MatrixSquareSymmetric>& sill,
  double* crit_arg) const
{
  int allpos;
  double temp, crit, crit_mem, value;
  VectorDouble valpro;
  const MatrixSquareGeneral* vecpro;

  /*******************/
  /* Initializations */
  /*******************/

  double sum  = 0.;
  double sum1 = 0.;
  double sum2 = 0.;
  int nvs2    = nvar * (nvar + 1) / 2;

  /* Core allocation */

  MatrixRectangular mp(nvs2, npadir);
  std::vector<MatrixRectangular> fk;
  fk.reserve(ncova);
  for (int icova = 0; icova < ncova; icova++)
    fk.push_back(MatrixRectangular(nvs2, npadir));
  MatrixSquareSymmetric cc(nvar);

  std::vector<MatrixSquareSymmetric> aic;
  aic.reserve(ncova);
  for (int icova = 0; icova < ncova; icova++)
    aic.push_back(MatrixSquareSymmetric(nvar));
  std::vector<MatrixSquareSymmetric> alphak;
  alphak.reserve(ncova);
  for (int icova = 0; icova < ncova; icova++)
    alphak.push_back(MatrixSquareSymmetric(nvar));

  /********************/
  /* Pre-calculations */
  /********************/

  int ijvar = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
      for (int ipadir = 0; ipadir < npadir; ipadir++)
      {
        mp.setValue(ijvar, ipadir, 0.);
        for (int icov = 0; icov < ncova; icov++)
          mp.setValue(ijvar, ipadir,
                      mp.getValue(ijvar, ipadir) +
                        sill[icov].getValue(ivar, jvar) *
                          ge[icov].getValue(ijvar, ipadir));
      }

  for (int icov = 0; icov < ncova; icov++)
  {
    ijvar = 0;
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
      {
        sum1 = sum2 = 0;
        aic[icov].setValue(ivar, jvar, 0.);
        for (int ipadir = 0; ipadir < npadir; ipadir++)
        {
          if (FFFF(WT(ijvar, ipadir))) continue;
          temp = WT(ijvar, ipadir) * ge[icov].getValue(ijvar, ipadir);
          fk[icov].setValue(ijvar, ipadir, temp);
          sum1 += temp * GG(ijvar, ipadir);
          sum2 += temp * ge[icov].getValue(ijvar, ipadir);
        }
        alphak[icov].setValue(ivar, jvar, 1. / sum2);
        aic[icov].setValue(ivar, jvar,
                           sum1 * alphak[icov].getValue(ivar, jvar));
      }
  }

  crit  = 0.;
  ijvar = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
      for (int ipadir = 0; ipadir < npadir; ipadir++)
      {
        if (FFFF(WT(ijvar, ipadir))) continue;
        temp  = GG(ijvar, ipadir) - mp.getValue(ijvar, ipadir);
        value = (ivar != jvar) ? 2. : 1.;
        crit += value * WT(ijvar, ipadir) * temp * temp;
      }

  /***********************/
  /* Iterative procedure */
  /***********************/

  int iter;
  for (iter = 0; iter < mauto->getMaxiter(); iter++)
  {

    /* Loop on the elementary structures */

    for (int icov = 0; icov < ncova; icov++)
    {

      /* Establish the coregionalization matrix */

      ijvar = 0;
      for (int ivar = 0; ivar < nvar; ivar++)
      {
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
        {
          sum = 0;
          for (int ipadir = 0; ipadir < npadir; ipadir++)
          {
            mp.setValue(ijvar, ipadir,
                        mp.getValue(ijvar, ipadir) -
                          sill[icov].getValue(ivar, jvar) *
                            ge[icov].getValue(ijvar, ipadir));
            sum +=
              fk[icov].getValue(ijvar, ipadir) * mp.getValue(ijvar, ipadir);
          }
          value = aic[icov].getValue(ivar, jvar) -
                  alphak[icov].getValue(ivar, jvar) * sum;
          cc.setValue(ivar, jvar, value);
          cc.setValue(jvar, ivar, value);
        }
      }
      /* Computing and sorting the eigen values and eigen vectors */

      if (cc.computeEigen()) return 1;
      valpro = cc.getEigenValues();
      vecpro = cc.getEigenVectors();

      int kvar = 0;
      allpos   = 1;
      while ((kvar < nvar) && allpos)
      {
        if (valpro[kvar++] < 0) allpos = 0;
      }

      /* Calculate the new coregionalization matrix */

      ijvar = 0;
      for (int ivar = 0; ivar < nvar; ivar++)
      {
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
        {
          if (allpos)
          {
            sill[icov].setValue(ivar, jvar, cc.getValue(ivar, jvar));
          }
          else
          {
            sum = 0.;
            for (kvar = 0; kvar < nvar; kvar++)
              sum += (MAX(valpro[kvar], 0.) * vecpro->getValue(ivar, kvar) *
                      vecpro->getValue(jvar, kvar));
            sill[icov].setValue(ivar, jvar, sum);
          }
          for (int ipadir = 0; ipadir < npadir; ipadir++)
            mp.setValue(ijvar, ipadir,
                        mp.getValue(ijvar, ipadir) +
                          sill[icov].getValue(ivar, jvar) *
                            ge[icov].getValue(ijvar, ipadir));
        }
      }
    }

    /* Update the global criterion */

    crit_mem = crit;
    crit     = 0.;
    for (int ipadir = 0; ipadir < npadir; ipadir++)
    {
      ijvar = 0;
      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
        {
          if (FFFF(WT(ijvar, ipadir))) continue;
          temp  = GG(ijvar, ipadir) - mp.getValue(ijvar, ipadir);
          value = (ivar != jvar) ? 2. : 1.;
          crit += value * WT(ijvar, ipadir) * temp * temp;
        }
    }

    /* Stopping criterion */

    if (ABS(crit) < mauto->getTolred() ||
        ABS(crit - crit_mem) / ABS(crit) < mauto->getTolred())
      break;
  }

  *crit_arg = crit;
  return (0);
}
