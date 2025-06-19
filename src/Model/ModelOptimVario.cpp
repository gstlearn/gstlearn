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
#include "Model/ModelOptimVario.hpp"

#include "Model/AModelFitSills.hpp"
#include "Model/ModelFitSillsVario.hpp"
#include "geoslib_define.h"

#include "Model/Model.hpp"
#include "Variogram/Vario.hpp"
#include <cstddef>

#define IJDIR(ijvar, ipadir) ((ijvar) * npadir + (ipadir))
#define WT(ijvar, ipadir)    wt[IJDIR(ijvar, ipadir)]

ModelOptimVario::ModelOptimVario(ModelGeneric* model,
                                 const Constraints* constraints,
                                 const ModelOptimParam& mop)
  : AModelOptim(model)
  , _mop(mop)
  , _constraints(constraints)
  , _vario()
  , _lags()
{
  setAuthorizedAnalyticalGradients(false);
}

ModelOptimVario::ModelOptimVario(const ModelOptimVario& m)
  : AModelOptim(m)
  , _mop(m._mop)
  , _constraints(m._constraints)
  , _calcmode(m._calcmode)
  , _vario(m._vario)
  , _lags(m._lags)
{
  setAuthorizedAnalyticalGradients(m.getAuthorizedAnalyticalGradients());
}

ModelOptimVario& ModelOptimVario::operator=(const ModelOptimVario& m)
{
  if (this != &m)
  {
    AModelOptim::operator=(m);
    _mop         = m._mop;
    _constraints = m._constraints;
    _calcmode    = m._calcmode;
    _vario       = m._vario;
    _lags        = m._lags;

    setAuthorizedAnalyticalGradients(m.getAuthorizedAnalyticalGradients());
  }
  return (*this);
}

ModelOptimVario::~ModelOptimVario()
{
}

bool ModelOptimVario::_checkConsistency()
{
  if (_vario->getNDim() != (int)_model->getNDim())
  {
    messerr("'_vario'(%d) and '_model'(%d) should have same Space Dimension",
            _vario->getNDim(), _model->getNDim());
    return false;
  }
  if (_vario->getNVar() != _model->getNVar())
  {
    messerr("'_vario'(%d) and '_model'(%d) should have same number of Variables",
            _vario->getNVar(), _model->getNVar());
    return false;
  }
  return true;
}

int ModelOptimVario::_buildExperimental()
{
  if (_vario == nullptr)
  {
    messerr("Argument 'vario' must be defined beforehand");
    return 1;
  }

  // Clean previous contents
  _lags.clear();

  int nvar = _vario->getNVar();
  int ndim = _vario->getNDim();
  VectorDouble dd(ndim);

  for (int idir = 0, ndir = _vario->getNDir(); idir < ndir; idir++)
  {
    for (int ilag = 0, nlag = _vario->getNLag(idir); ilag < nlag; ilag++)
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
            int iad    = _vario->getDirAddress(idir, ivar, jvar, ilag, false, 1);
            int jad    = _vario->getDirAddress(idir, ivar, jvar, ilag, false, -1);
            double c00 = _vario->getC00(idir, ivar, jvar);
            double n1  = _vario->getSwByIndex(idir, iad);
            double n2  = _vario->getSwByIndex(idir, jad);
            if (n1 + n2 > 0)
            {
              double g1 = _vario->getGgByIndex(idir, iad);
              double g2 = _vario->getGgByIndex(idir, jad);
              if (_vario->isLagCorrect(idir, iad) && _vario->isLagCorrect(idir, jad))
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
            int iad = _vario->getDirAddress(idir, ivar, jvar, ilag, false, 1);
            if (_vario->isLagCorrect(idir, iad))
            {
              gg   = _vario->getGgByIndex(idir, iad);
              dist = ABS(_vario->getHhByIndex(idir, iad));
            }
          }

          /* Define the item of the StrExp array (if defined) */

          if (FFFF(gg)) continue;
          OneLag onelag = _createOneLag(ndim, idir, ivar, jvar, gg, dist);
          _lags.push_back(onelag);
        }
    }
  }

  // Update the weight
  VectorDouble wt = _vario->computeWeightsFromVario(_mop.getWmode());
  int npadir      = _vario->getTotalLagsPerDirection();
  int ecr         = 0;
  int ipadir      = 0;

  for (int idir = 0, ndir = _vario->getNDir(); idir < ndir; idir++)
    for (int ilag = 0, nlag = _vario->getNLag(idir); ilag < nlag; ilag++, ipadir++)
    {
      int ijvar = 0;
      for (int ivar = ijvar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
          _lags[ecr++]._weight = WT(ijvar, ipadir);
    }

  return 0;
}

ModelOptimVario::OneLag ModelOptimVario::_createOneLag(int ndim,
                                                       int idir,
                                                       int ivar,
                                                       int jvar,
                                                       double gg,
                                                       double dist) const
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

ModelOptimVario* ModelOptimVario::createForOptim(ModelGeneric* model,
                                                 const Vario* vario,
                                                 const Constraints* constraints,
                                                 const ModelOptimParam& mop)
{
  auto* optim = new ModelOptimVario(model, constraints, mop);

  optim->_vario = vario;

  // Constitute the experimental material (using '_vario')
  if (optim->_buildExperimental())
  {
    delete optim;
    return nullptr;
  }

  // Check consistency
  if (!optim->_checkConsistency())
  {
    delete optim;
    return nullptr;
  }

  // Instantiate Goulard algorithm (optional)
  if (mop.getFlagGoulard())
  {
    ModelCovList* mcv = dynamic_cast<ModelCovList*>(model);
    if (mcv != nullptr)
    {
      mcv->setFitSills(ModelFitSillsVario::createForOptim(vario, model, constraints, mop));
      if (mcv->getFitSills() == nullptr)
      {
        delete optim;
        return nullptr;
      }
    }
  }

  // Perform the Fitting in terms of variograms
  optim->_calcmode.setAsVario(true);

  return optim;
}

double ModelOptimVario::computeCost(bool verbose)
{
  DECLARE_UNUSED(verbose);

  // Evaluate the Cost function
  int nlags    = (int)_lags.size();
  double score = 0.;
  SpacePoint origin;
  _resid.resize(nlags);
  for (int ilag = 0; ilag < nlags; ilag++)
  {
    const OneLag& lag = _lags[ilag];
    double vtheo      = _model->evalCov(origin, lag._P, lag._ivar, lag._jvar, &_calcmode);
    
    double resid      = lag._gg - vtheo;
    score += lag._weight * resid * resid;
    _resid[ilag] = lag._weight * resid;
  }
  return score;
}

void ModelOptimVario::evalGrad(vect res)
{

  auto gradcov = _model->getGradients();
  int nlags    = (int)_lags.size();
  SpacePoint origin;

  for (size_t i = 0; i < gradcov.size(); i++)
  {
    res[i] = 0.;
    for (int ilag = 0; ilag < nlags; ilag++)
    {
      const OneLag& lag = _lags[ilag];
      double dvtheo = gradcov[i](origin, lag._P, lag._ivar, lag._jvar, &_calcmode);
      res[i] += -2. * _resid[ilag] * dvtheo;
    }
  }

}
