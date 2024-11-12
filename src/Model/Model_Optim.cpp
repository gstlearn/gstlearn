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
#include "Model/Model_Optim.hpp"

#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "geoslib_define.h"

#include "Model/Model.hpp"
#include "Variogram/Vario.hpp"

#include <nlopt.h>

#define IJDIR(ijvar, ipadir) ((ijvar)*npadir + (ipadir))
#define WT(ijvar, ipadir)    wt[IJDIR(ijvar, ipadir)]

Model_Optim::Model_Optim(const Vario* vario)
  : _vario(vario)
  , _wmode(1)
  , _tabval()
  , _tablow()
  , _tabupp()
  , _algorithm()
{
}

Model_Optim::Model_Optim(const Model_Optim& m)
  : _vario(m._vario)
  , _wmode(m._wmode)
  , _tabval(m._tabval)
  , _tablow(m._tablow)
  , _tabupp(m._tabupp)
  , _algorithm()
{
  _algorithm._model  = m._algorithm._model;
  _algorithm._lags   = m._algorithm._lags;
  _algorithm._params = m._algorithm._params;
}

Model_Optim& Model_Optim::operator=(const Model_Optim &m)
{
  if (this != &m)
  {
    _vario  = m._vario;
    _wmode  = m._wmode;
    _tabval = m._tabval;
    _tablow = m._tablow;
    _tabupp = m._tabupp;

    _algorithm._model  = m._algorithm._model;
    _algorithm._lags   = m._algorithm._lags;
    _algorithm._params = m._algorithm._params;
  }
  return (*this);
}

Model_Optim::~Model_Optim()
{
}

int Model_Optim::fit(Model* model, int wmode)
{
  _wmode = wmode;
  _algorithm._model = model;
  
  // Constitute the experimental material (using '_vario')
  if (_buildExperimental()) return 1;

  // Constitute the list of parameters
  if (_buildModelParamList()) return 1;

  // Check consistency
  if (_vario->getDimensionNumber() != _algorithm._model->getDimensionNumber())
  {
    messerr("'_vario'(%d) and '_model'(%d) should have same Space Dimension",
            _vario->getDimensionNumber(), _algorithm._model->getDimensionNumber());
    return 1;
  }
  if (_vario->getVariableNumber() != _algorithm._model->getVariableNumber())
  {
    messerr("'_vario'(%d) and '_model'(%d) should have same number of Variables",
            _vario->getVariableNumber(), _algorithm._model->getVariableNumber());
    return 1;
  }

  // Perform the optimization
  int npar = (int)_algorithm._params.size();
  nlopt_opt opt = nlopt_create(NLOPT_LN_NELDERMEAD, npar);

  // Define the bounds
  nlopt_set_lower_bounds(opt, _tablow.data());
  nlopt_set_upper_bounds(opt, _tabupp.data());

  // Define the cost function
  nlopt_set_min_objective(opt, evalCost, &_algorithm);

  // Perform the optimization (store the minimized value in 'minf')
  double minf;
  nlopt_optimize(opt, _tabval.data(), &minf);
  return 0;
}

int Model_Optim::_buildExperimental()
{
  if (_vario == nullptr)
  {
    messerr("Argument 'vario' must be defined beforehand");
    return 1;
  }

  // Clean previous contents
  _algorithm._lags.clear();

  int nvar = _vario->getVariableNumber();
  int ndim = _vario->getDimensionNumber();
  VectorDouble dd(ndim);

  for (int idir = 0, ndir = _vario->getDirectionNumber(); idir < ndir; idir++)
  {
    for (int ipas = 0, npas = _vario->getLagNumber(idir); ipas < npas; ipas++)
    {
      int ijvar = 0;
      for (int ivar = ijvar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
        {

          /* Calculate the variogram value */

          double dist = 0.;
          double gg   = TEST;
          if (_vario->getFlagAsym())
          {
            int iad = _vario->getDirAddress(idir, ivar, jvar, ipas, false, 1);
            int jad = _vario->getDirAddress(idir, ivar, jvar, ipas, false, -1);
            double c00 = _getC00(idir, ivar, jvar);
            double n1  = _vario->getSwByIndex(idir, iad);
            double n2  = _vario->getSwByIndex(idir, jad);
            if (n1 + n2 > 0)
            {
              double g1 = _vario->getGgByIndex(idir, iad);
              double g2 = _vario->getGgByIndex(idir, jad);
              if (_isLagCorrect(idir, iad) && _isLagCorrect(idir, jad))
              {
                gg   = c00 - (n1 * g1 + n2 * g2) / (n1 + n2);
                dist = (ABS(_vario->getHhByIndex(idir, iad)) +
                        ABS(_vario->getHhByIndex(idir, jad))) /
                       2.;
              }
            }
          }
          else
          {
            int iad = _vario->getDirAddress(idir, ivar, jvar, ipas, false, 1);
            if (_isLagCorrect(idir, iad))
            {
              gg   = _vario->getGgByIndex(idir, iad);
              dist = ABS(_vario->getHhByIndex(idir, iad));
            }
          }

          /* Define the item of the StrExp array (if defined) */

          if (FFFF(gg)) continue;
          OneLag onelag = _createOneLag(ndim, idir, ivar, jvar, gg, dist);
          _algorithm._lags.push_back(onelag);
        }
    }
  }

  // Update the weight
  VectorDouble wt = _computeWeight();
  int npadir      = _getTotalLagsPerDirection();
  int ecr         = 0;
  int ipadir      = 0;
  for (int idir = 0, ndir = _vario->getDirectionNumber(); idir < ndir; idir++)
  {
    for (int ipas = 0, npas = _vario->getLagNumber(idir); ipas < npas;
         ipas++, ipadir++)
    {
      int ijvar = 0;
      for (int ivar = ijvar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
          _algorithm._lags[ecr]._weight = WT(ijvar, ipadir);
    }
  }

  return 0;
}

OneLag Model_Optim::_createOneLag(
  int ndim, int idir, int ivar, int jvar, double gg, double dist)
{
  OneLag onelag;
  onelag._ivar   = ivar;
  onelag._jvar   = jvar;
  onelag._gg     = gg;
  onelag._weight = 1.;
  VectorDouble dd(ndim);
  for (int idim = 0; idim < ndim; idim++)
    dd[idim] = dist * _vario->getCodir(idir, idim);
  onelag._P.setCoords(dd);
  return onelag;
}

double Model_Optim::_getC00(int idir, int ivar, int jvar)
{
  int iad0   = _vario->getDirAddress(idir, ivar, jvar, 0, false, 0);
  int iad    = iad0;
  double c00 = _vario->getSwByIndex(idir, iad);
  if (!isZero(c00) || _vario->getSwByIndex(idir, iad) > 0) return c00;

  for (int ipas = 0, npas = _vario->getLagNumber(idir); ipas < npas; ipas++)
  {
    iad = _vario->getDirAddress(idir, ivar, jvar, ipas, false, 1);
    if (!isZero(_vario->getGgByIndex(idir, iad)))
      return _vario->getGgByIndex(idir, iad);
    iad = _vario->getDirAddress(idir, ivar, jvar, ipas, false, -1);
    if (!isZero(_vario->getGgByIndex(idir, iad)))
      return _vario->getGgByIndex(idir, iad);
  }
  iad = iad0;
  return (_vario->getGgByIndex(idir, iad));
}

int Model_Optim::_buildModelParamList()
{
  if (_algorithm._model == nullptr)
  {
    messerr("Argument '_model' must be defined beforehand");
    return 1;
  }

  // Clear previous contents
  _tabval.clear();
  _tablow.clear();
  _tabupp.clear();

  // Loop on the covariances
  int nvar = _algorithm._model->getVariableNumber();
  for (int icov = 0, ncov = _algorithm._model->getCovaNumber(); icov < ncov; icov++)
  {
    const CovAniso* cova = _algorithm._model->getCova(icov);

    if (cova->hasRange())
    {
      // Add the 'Range' (scalar) attribute
      _addOneModelParam(icov, 0, 0, EPSILON2, 1.e30);
    }
    if (cova->hasParam())
    {
      // Add the 'Param' (scalar) attribute
      _addOneModelParam(icov, 2, 0, 0., 2.);
    }

    // Add the 'Sill' (vectorial)attribute
    int ijvar = 0;
    for (int ivar = ijvar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
        _addOneModelParam(icov, 1, ijvar, -1.e30, 1.e30);
  }
  return 0;
}

void Model_Optim::_addOneModelParam(int icov,
                                    int type,
                                    int rank,
                                    double lbound,
                                    double ubound)
{
  OneParam param;
  param._icov = icov;
  param._type = type;
  param._rank = rank;
  _algorithm._params.push_back(param);
  _tabval.push_back(1.);
  _tablow.push_back(lbound);
  _tabupp.push_back(ubound);
}

bool Model_Optim::_isLagCorrect(int idir, int k)
{
  double hh = _vario->getHhByIndex(idir, k);
  if (isZero(hh) || FFFF(hh)) return false;
  double sw = _vario->getSwByIndex(idir, k);
  if (isZero(sw) || FFFF(sw)) return false;
  double gg = _vario->getGgByIndex(idir, k);
  return !FFFF(gg);
}

/*****************************************************************************/
/*!
 **  Calculates the weighting factors for each experimental variogram
 **  value of each directional variogram
 **
 *****************************************************************************/
VectorDouble Model_Optim::_computeWeight()
{
  int ipadir;

  /* Initializations */

  int ndir           = _vario->getDirectionNumber();
  int nvar           = _vario->getVariableNumber();
  int nvs2           = nvar * (nvar + 1) / 2;
  VectorDouble count = _computeWeightPerDirection();

  int npadir = _getTotalLagsPerDirection();
  VectorDouble wt(npadir * nvs2, 0.);

  /* Determine the count of significant directions */

  for (int idir = 0; idir < ndir; idir++)
  {
    count[idir] = 0.;
    int npas    = _vario->getLagNumber(idir);
    for (int ipas = 0; ipas < npas; ipas++)
      for (int ijvar = 0; ijvar < nvs2; ijvar++)
      {
        int shift = ijvar * _vario->getLagTotalNumber(idir);
        if (_vario->getFlagAsym())
        {
          int iad   = shift + _vario->getLagNumber(idir) + ipas + 1;
          int jad   = shift + _vario->getLagNumber(idir) - ipas - 1;
          double n1 = _vario->getSwByIndex(idir, iad);
          double n2 = _vario->getSwByIndex(idir, jad);
          if (_isLagCorrect(idir, iad)) count[idir] += n1;
          if (_isLagCorrect(idir, jad)) count[idir] += n2;
        }
        else
        {
          int iad   = shift + ipas;
          double nn = _vario->getSwByIndex(idir, iad);
          if (_isLagCorrect(idir, iad)) count[idir] += nn;
        }
      }
  }

  switch (_wmode)
  {
    case 1:
      ipadir = 0;
      for (int idir = 0; idir < ndir; idir++)
      {
        int npas = _vario->getLagNumber(idir);
        for (int ipas = 0; ipas < npas; ipas++, ipadir++)
        {
          if (isZero(count[idir])) continue;
          for (int ijvar = 0; ijvar < nvs2; ijvar++)
          {
            int shift = ijvar * _vario->getLagTotalNumber(idir);
            if (_vario->getFlagAsym())
            {
              int iad = shift + _vario->getLagNumber(idir) + ipas + 1;
              int jad = shift + _vario->getLagNumber(idir) - ipas - 1;
              if (_isLagCorrect(idir, iad) && _isLagCorrect(idir, jad))
                WT(ijvar, ipadir) = count[idir];
            }
            else
            {
              int iad = shift + ipas;
              if (_isLagCorrect(idir, iad)) WT(ijvar, ipadir) = count[idir];
            }
          }
        }
      }
      break;

    case 2:
      ipadir = 0;
      for (int idir = 0; idir < ndir; idir++)
      {
        int npas = _vario->getLagNumber(idir);
        for (int ipas = 0; ipas < npas; ipas++, ipadir++)
        {
          if (isZero(count[idir])) continue;
          for (int ijvar = 0; ijvar < nvs2; ijvar++)
          {
            int shift = ijvar * _vario->getLagTotalNumber(idir);
            if (_vario->getFlagAsym())
            {
              int iad = shift + _vario->getLagNumber(idir) + ipas + 1;
              int jad = shift + _vario->getLagNumber(idir) - ipas - 1;
              if (_isLagCorrect(idir, iad) || _isLagCorrect(idir, jad))
                continue;
              double n1 = _vario->getSwByIndex(idir, iad);
              double n2 = _vario->getSwByIndex(idir, jad);
              double d1 = ABS(_vario->getHhByIndex(idir, iad));
              double d2 = ABS(_vario->getHhByIndex(idir, jad));
              if (d1 > 0 && d2 > 0)
                WT(ijvar, ipadir) =
                  sqrt((n1 + n2) * (n1 + n2) / (n1 * d1 + n2 * d2) / 2.);
            }
            else
            {
              int iad = shift + ipas;
              if (_isLagCorrect(idir, iad)) continue;
              double nn = _vario->getSwByIndex(idir, iad);
              double dd = ABS(_vario->getHhByIndex(idir, iad));
              if (dd > 0) WT(ijvar, ipadir) = nn / dd;
            }
          }
        }
      }
      break;

    case 3:
      ipadir = 0;
      for (int idir = 0; idir < ndir; idir++)
      {
        int npas = _vario->getLagNumber(idir);
        for (int ipas = 0; ipas < npas; ipas++, ipadir++)
        {
          if (isZero(count[idir])) continue;
          for (int ijvar = 0; ijvar < nvs2; ijvar++)
          {
            int shift = ijvar * _vario->getLagTotalNumber(idir);
            if (_vario->getFlagAsym())
            {
              int iad = shift + _vario->getLagNumber(idir) + ipas + 1;
              int jad = shift + _vario->getLagNumber(idir) - ipas - 1;
              if (_isLagCorrect(idir, iad) && _isLagCorrect(idir, jad))
                WT(ijvar, ipadir) = 1. / _vario->getLagNumber(idir);
            }
            else
            {
              int iad = shift + ipas;
              if (_isLagCorrect(idir, iad))
                WT(ijvar, ipadir) = 1. / _vario->getLagNumber(idir);
            }
          }
        }
      }
      break;

    default:
      ipadir = 0;
      for (int idir = 0; idir < ndir; idir++)
      {
        int npas = _vario->getLagNumber(idir);
        for (int ipas = 0; ipas < npas; ipas++, ipadir++)
        {
          if (isZero(count[idir])) continue;
          for (int ijvar = 0; ijvar < nvs2; ijvar++)
          {
            int shift = ijvar * _vario->getLagTotalNumber(idir);
            if (_vario->getFlagAsym())
            {
              int iad = shift + _vario->getLagNumber(idir) + ipas + 1;
              int jad = shift + _vario->getLagNumber(idir) - ipas - 1;
              if (_isLagCorrect(idir, iad) && _isLagCorrect(idir, jad))
                WT(ijvar, ipadir) = 1.;
            }
            else
            {
              int iad = shift + ipas;
              if (_isLagCorrect(idir, iad)) WT(ijvar, ipadir) = 1.;
            }
          }
        }
      }
      break;
  }

  /* Scaling by direction and by variable */

  for (int ijvar = 0; ijvar < nvs2; ijvar++)
  {
    ipadir = 0;
    for (int idir = 0; idir < ndir; idir++)
    {
      double total = 0.;
      int npas     = _vario->getLagNumber(idir);
      for (int ipas = 0; ipas < npas; ipas++, ipadir++)
      {
        if (isZero(count[idir])) continue;
        if (WT(ijvar, ipadir) > 0 && !FFFF(WT(ijvar, ipadir)))
          total += WT(ijvar, ipadir);
      }
      if (isZero(total)) continue;
      ipadir -= _vario->getLagNumber(idir);
      for (int ipas = 0, npas = _vario->getLagNumber(idir); ipas < npas;
           ipas++, ipadir++)
      {
        if (isZero(count[idir])) continue;
        if (WT(ijvar, ipadir) > 0 && !FFFF(WT(ijvar, ipadir)))
          WT(ijvar, ipadir) /= total;
      }
    }
  }

  /* Scaling by variable variances */

  int ijvar0 = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++, ijvar0++)
    {
      double ratio =
        (_vario->getVar(ivar, jvar) > 0 && _vario->getVar(jvar, ivar) > 0)
          ? sqrt(_vario->getVar(ivar, jvar) * _vario->getVar(jvar, ivar))
          : 1.;
      ipadir = 0;
      for (int idir = 0; idir < ndir; idir++)
      {
        int npas = _vario->getLagNumber(idir);
        for (int ipas = 0; ipas < npas; ipas++, ipadir++)
          if (!FFFF(WT(ijvar0, ipadir))) WT(ijvar0, ipadir) /= ratio;
      }
    }

  return wt;
}

/*****************************************************************************/
/*!
 **  Calculates the sum of weighting factors per direction
 **
 *****************************************************************************/
VectorDouble Model_Optim::_computeWeightPerDirection()
{
  int ndir = _vario->getDirectionNumber();
  int nvar = _vario->getVariableNumber();
  int nvs2 = nvar * (nvar + 1) / 2;
  VectorDouble count(ndir);

  /* Determine the count of significant directions */

  for (int idir = 0; idir < ndir; idir++)
  {
    count[idir] = 0.;
    for (int ipas = 0, npas = _vario->getLagNumber(idir); ipas < npas; ipas++)
      for (int ijvar = 0; ijvar < nvs2; ijvar++)
      {
        int shift = ijvar * _vario->getLagTotalNumber(idir);
        if (_vario->getFlagAsym())
        {
          int iad   = shift + _vario->getLagNumber(idir) + ipas + 1;
          int jad   = shift + _vario->getLagNumber(idir) - ipas - 1;
          double n1 = _vario->getSwByIndex(idir, iad);
          double n2 = _vario->getSwByIndex(idir, jad);
          if (_isLagCorrect(idir, iad)) count[idir] += n1;
          if (_isLagCorrect(idir, jad)) count[idir] += n2;
        }
        else
        {
          int iad   = shift + ipas;
          double nn = _vario->getSwByIndex(idir, iad);
          if (_isLagCorrect(idir, iad)) count[idir] += nn;
        }
      }
  }

  return count;
}

int Model_Optim::_getTotalLagsPerDirection() const
{
  int npatot = 0;
  int ndir   = _vario->getDirectionNumber();
  for (int idir = 0; idir < ndir; idir++)
    npatot += _vario->getLagTotalNumber(idir);
  return npatot;
}

void Model_Optim::_patchModel(Algorithm* algo, const double* current)
{
  // Initializations
  int nvar = algo->_model->getVariableNumber();
  int nvs2 = nvar * (nvar + 1) / 2;
  int nparams = (int)algo->_params.size();

  // Loop on the parameters
  int iparam = 0;
  while (iparam < nparams)
  {
    const OneParam& param = algo->_params[iparam];
    int icov              = param._icov;
    CovAniso* cova        = algo->_model->getCova(icov);

    if (param._type == 0)
    {
      // Range
      cova->setRangeIsotropic(current[iparam++]);
    }

    if (param._type == 1)
    {
      // Sill (through the AIC matrix)
      MatrixSquareGeneral aic(nvar);
      int ijvar = 0;
      for (int ivar = ijvar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
          aic.setValue(ivar, jvar, current[iparam++]);
      if (ijvar == nvs2)
      {
        MatrixSquareSymmetric sills(nvar);
        sills.prodNormMatVecInPlace(aic, VectorDouble());
        cova->setSill(sills);
      }
    }

    if (param._type == 2)
    {
      // Set the 'param' attribute
      cova->setParam(current[iparam++]);
    }
  }
}

double Model_Optim::evalCost(unsigned int nparams,
                             const double* current,
                             double* grad,
                             void* my_func_data)
{
  DECLARE_UNUSED(nparams);
  DECLARE_UNUSED(grad);
  Algorithm* algorithm = (Algorithm*) my_func_data;
  
  // Update the Model
  _patchModel(algorithm, current);

  // Evaluate the Cost function
  int nlags = (int) algorithm->_lags.size();
  double total = 0.;
  SpacePoint origin;
  CovCalcMode calcmod;
  calcmod.setAsVario(true);
  for (int ipas = 0; ipas < nlags; ipas++)
  {
    const OneLag& lag = algorithm->_lags[ipas];
    double vexp        = lag._gg;
    double vtheo =
      algorithm->_model->eval(origin, lag._P, lag._ivar, lag._jvar, &calcmod);
    double delta = vexp - vtheo;
    total += lag._weight * delta * delta;
  }
  
  return total;
}
