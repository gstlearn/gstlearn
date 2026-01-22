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
#include "Model/ModelFitSillsVMap.hpp"
#include "geoslib_define.h"
#include "geoslib_old_f.h"

#include "Db/DbGrid.hpp"
#include "Model/Model.hpp"
#include "Model/ModelOptimVMap.hpp"

#define IJDIR(ijvar, ipadir) ((ijvar) * npadir + (ipadir))
#define _WT(ijvar, ipadir)   _wt[IJDIR(ijvar, ipadir)]
#define _GG(ijvar, ipadir)   _gg[IJDIR(ijvar, ipadir)]

namespace gstlrn
{
ModelOptimVMap::ModelOptimVMap(ModelGeneric* model,
                               const Constraints* constraints,
                               const ModelOptimParam& mop)
  : AModelOptim(model)
  , _mop(mop)
  , _constraints(constraints)
  , _calcmode()
  , _dbmap(nullptr)
  , _indg1()
  , _indg2()
  , _ndim(0)
  , _nvar(0)
  , _nech(0)
  , _npadir(0)
{
  setAuthorizedAnalyticalGradients(true);
}

ModelOptimVMap::ModelOptimVMap(const ModelOptimVMap& m)
  : AModelOptim(m)
  , _mop(m._mop)
  , _constraints(m._constraints)
  , _calcmode(m._calcmode)
  , _dbmap(m._dbmap)
  , _indg1(m._indg1)
  , _indg2(m._indg2)
  , _ndim(m._ndim)
  , _nvar(m._nvar)
  , _nech(m._nech)
  , _npadir(m._npadir)
{
  setAuthorizedAnalyticalGradients(m.getAuthorizedAnalyticalGradients());
}

ModelOptimVMap& ModelOptimVMap::operator=(const ModelOptimVMap& m)
{
  if (this != &m)
  {
    AModelOptim::operator=(m);
    _mop         = m._mop;
    _constraints = m._constraints;
    _calcmode    = m._calcmode;
    _dbmap       = m._dbmap;
    _indg1       = m._indg1;
    _indg2       = m._indg2;
    _ndim        = m._ndim;
    _nvar        = m._nvar;
    _nech        = m._nech;
    _npadir      = m._npadir;

    setAuthorizedAnalyticalGradients(m.getAuthorizedAnalyticalGradients());
  }
  return (*this);
}

ModelOptimVMap::~ModelOptimVMap()
{
}

bool ModelOptimVMap::_checkConsistency()
{
  if (_dbmap == nullptr)
  {
    messerr("You must have defined 'dbmap' beforehand");
    return false;
  }
  Id nvar     = _dbmap->getNLoc(ELoc::Z);
  size_t ndim = _dbmap->getNLoc(ELoc::X);

  if (_model->getNVar() != nvar)
  {
    messerr("Number of variables in Dbmap (%d) must match the one in Model (%d)",
            nvar, _model->getNVar());
    return false;
  }
  if (_model->getNDim() != ndim)
  {
    messerr("'_dbmap'(%d) and '_model'(%d) should have same Space Dimensions",
            ndim, _model->getNDim());
    return false;
  }
  // if (_constraints->isConstraintSillDefined())
  // {
  //   if (!_optvar->getFlagGoulardUsed())
  //   {
  //     messerr("When Constraints on the sum of Sills are defined");
  //     messerr("The Goulard option must be switched ON");
  //     return false;
  //   }
  //   if (!FFFF(_constraints->getConstantSillValue()))
  //     _constraints->expandConstantSill(nvar);
  // }

  _indg1.fill(0., ndim);
  _indg2.fill(0., ndim);
  return true;
}

Id ModelOptimVMap::_getDimensions()
{
  if (_dbmap == nullptr)
  {
    messerr("You must have defined 'dbmap' beforehand");
    return 1;
  }
  Id nbexp  = 0;
  Id npadir = 0;
  Id nvar   = _dbmap->getNLoc(ELoc::Z);
  Id nech   = _dbmap->getNSample();
  Id nvs2   = nvar * (nvar + 1) / 2;

  /* Calculate the total number of lags */

  for (Id iech = 0; iech < nech; iech++)
  {
    Id ndef = 0;
    for (Id ijvar = 0; ijvar < nvs2; ijvar++)
      if (!FFFF(_dbmap->getZVariable(iech, ijvar))) ndef++;
    nbexp += ndef;
    if (ndef > 0) npadir++;
  }

  /* Setting the return arguments */

  _npadir = npadir;

  if (nbexp <= 0)
  {
    messerr("No active experimental variogram map samples");
    return 1;
  }
  return 0;
}

double ModelOptimVMap::computeCost(bool flagPrint, bool verbose)
{
  DECLARE_UNUSED(flagPrint);
  DECLARE_UNUSED(verbose);

  // Evaluate the Cost function
  VectorDouble d0(_ndim);
  _dbmap->rankToIndice(_nech / 2, _indg1);
  for (Id idim = 0; idim < _ndim; idim++)
    d0[idim] = _indg1[idim] * _dbmap->getDX(idim);
  SpacePoint origin(d0);
  SpacePoint P = origin;

  /* Loop on the experimental conditions */
  double total = 0.;
  for (Id iech = 0; iech < _nech; iech++)
  {
    _dbmap->rankToIndice(iech, _indg2);
    for (Id idim = 0; idim < _ndim; idim++)
      d0[idim] = _indg2[idim] * _dbmap->getDX(idim);
    P.setCoords(d0);

    double dist = distance_intra(_dbmap, _nech / 2, iech, NULL);
    double wgt  = (dist > 0) ? 1. / dist : 0.;

    Id ijvar = 0;
    for (Id ivar = 0; ivar < _nvar; ivar++)
      for (Id jvar = 0; jvar <= ivar; jvar++, ijvar++)
      {
        double vexp = _dbmap->getZVariable(iech, ijvar);
        if (FFFF(vexp)) continue;
        double vtheo = _model->evalCov(origin, P, ivar, jvar, &_calcmode);
        double delta = vexp - vtheo;
        total += wgt * delta * delta;
      }
  }
  return total;
}

void ModelOptimVMap::evalGrad(vect res)
{
  VectorDouble d0(_ndim);
  _dbmap->rankToIndice(_nech / 2, _indg1);
  for (Id idim = 0; idim < _ndim; idim++)
    d0[idim] = _indg1[idim] * _dbmap->getDX(idim);
  SpacePoint origin(d0);
  SpacePoint P = origin;

  /* Loop on the experimental conditions */
  auto gradcov = _model->getCovGradients();

  for (size_t i = 0, ngrad = gradcov.size(); i < ngrad; i++)
  {
    double total = 0.;
    for (Id iech = 0; iech < _nech; iech++)
    {
      _dbmap->rankToIndice(iech, _indg2);
      for (Id idim = 0; idim < _ndim; idim++)
        d0[idim] = _indg2[idim] * _dbmap->getDX(idim);
      P.setCoords(d0);

      double dist = distance_intra(_dbmap, _nech / 2, iech, NULL);
      double wgt  = (dist > 0) ? 1. / dist : 0.;

      Id ijvar = 0;
      for (Id ivar = 0; ivar < _nvar; ivar++)
        for (Id jvar = 0; jvar <= ivar; jvar++, ijvar++)
        {
          double vexp = _dbmap->getZVariable(iech, ijvar);
          if (FFFF(vexp)) continue;
          double vtheo  = _model->evalCov(origin, P, ivar, jvar, &_calcmode);
          const auto& func = gradcov[i];
          double dvtheo = func(origin, P, ivar, jvar, &_calcmode);
          double delta  = vexp - vtheo;
          total += -2. * wgt * delta * dvtheo;
        }
    }
    res[i] = total;
  }
}

ModelOptimVMap* ModelOptimVMap::createForOptim(ModelGeneric* model,
                                               const DbGrid* dbmap,
                                               const Constraints* constraints,
                                               const ModelOptimParam& mop)
{
  auto* optim = new ModelOptimVMap(model, constraints, mop);

  MatrixSymmetric vars(model->getNVar());
  double hmax = dbmap->getExtensionDiagonal();
  optim->setEnvironment(vars, hmax);
  optim->_dbmap = dbmap;

  // Get internal dimension
  if (optim->_getDimensions())
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
      mcv->setFitSills(ModelFitSillsVMap::createForOptim(dbmap, model, constraints, mop));
      if (mcv->getFitSills() == nullptr)
      {
        delete optim;
        return nullptr;
      }
    }
  }

  // Perform the Fitting in terms of variograms
  optim->_calcmode.setAsVario(true);

  // Local variables
  optim->_ndim = dbmap->getNLoc(ELoc::X);
  optim->_nvar = dbmap->getNLoc(ELoc::Z);
  optim->_nech = dbmap->getNSample();

  return optim;
}
} // namespace gstlrn