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
#include "Model/ModelOptimSillsVMap.hpp"
#include "geoslib_define.h"
#include "geoslib_old_f.h"

#include "Model/ModelOptimVMap.hpp"
#include "Model/Model.hpp"
#include "Db/DbGrid.hpp"

#define IJDIR(ijvar, ipadir) ((ijvar)*npadir + (ipadir))
#define _WT(ijvar, ipadir)   _wt[IJDIR(ijvar, ipadir)]
#define _GG(ijvar, ipadir)   _gg[IJDIR(ijvar, ipadir)]

ModelOptimVMap::ModelOptimVMap(Model* model,
                               Constraints* constraints,
                               const Option_AutoFit& mauto,
                               const Option_VarioFit& optvar)
  : AModelOptim(model, constraints, mauto, optvar)
  , _vmapPart()
  , _optGoulard(model)
{
}

ModelOptimVMap::ModelOptimVMap(const ModelOptimVMap& m)
  : AModelOptim(m)
  , _vmapPart()
  , _optGoulard(m._optGoulard)
{
  _copyVMapPart(m._vmapPart);
}

ModelOptimVMap& ModelOptimVMap::operator=(const ModelOptimVMap& m)
{
  if (this != &m)
  {
    AModelOptim::operator=(m);
    _optGoulard = m._optGoulard;
    _copyVMapPart(m._vmapPart);
  }
  return (*this);
}

ModelOptimVMap::~ModelOptimVMap()
{
}

void ModelOptimVMap::_copyVMapPart(const VMap_Part& vmapPart)
{
  _vmapPart._dbmap = vmapPart._dbmap;
}

bool ModelOptimVMap::_checkConsistency()
{
  const Model* model = _modelPart._model;

  if (_vmapPart._dbmap == nullptr)
  {
    messerr("You must have defined 'dbmap' beforehand");
    return false;
  }
  int nvar = _vmapPart._dbmap->getNLoc(ELoc::Z);
  int ndim = _vmapPart._dbmap->getNLoc(ELoc::X);

  if (model->getNVar() != nvar)
  {
    messerr("Number of variables in Dbmap (%d) must match the one in Model (%d)",
            nvar, model->getNVar());
    return false;
  }
  if (model->getDimensionNumber() != ndim)
  {
    messerr(
      "'_dbmap'(%d) and '_model'(%d) should have same Space Dimensions",
      ndim, model->getDimensionNumber());
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

  _vmapPart._indg1.fill(0., ndim);
  _vmapPart._indg2.fill(0., ndim);
  return true;
}

double ModelOptimVMap::evalCost(unsigned int nparams,
                                 const double* current,
                                 double* /*grad*/,
                                 void* my_func_data)
{
  DECLARE_UNUSED(nparams);
  AlgorithmVMap* algorithm = static_cast<AlgorithmVMap*>(my_func_data);
  if (algorithm == nullptr) return TEST;
  Model_Part& modelPart           = algorithm->_modelPart;
  VMap_Part& vmapPart             = algorithm->_vmapPart;
  ModelOptimSillsVMap& optGoulard = algorithm->_goulardPart;
  const DbGrid* dbmap   = vmapPart._dbmap;
  int ndim              = dbmap->getNLoc(ELoc::X);
  int nvar              = dbmap->getNLoc(ELoc::Z);
  int nech              = dbmap->getNSample();

  // Update the Model
  _patchModel(modelPart, current);

  // Perform sill fitting using Goulard (optional)
  if (modelPart._optvar.getFlagGoulardUsed())
  {
    optGoulard.updateFromModel();
    optGoulard.fitPerform();
  }

  // Evaluate the Cost function
  double total = 0.;
  VectorDouble d0(ndim);
  dbmap->rankToIndice(nech / 2, vmapPart._indg1);

  // mode.setOrderVario(strmod->norder);

  /* Loop on the experimental conditions */
  for (int iech = 0; iech < nech; iech++)
  {
    dbmap->rankToIndice(iech, vmapPart._indg2);
    for (int idim = 0; idim < ndim; idim++)
      d0[idim] = (vmapPart._indg2[idim] - vmapPart._indg1[idim]) * dbmap->getDX(idim);

    int ijvar = 0;
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
      {
        double vexp = dbmap->getZVariable(iech, ijvar);
        if (FFFF(vexp)) continue;
        double vtheo = modelPart._model->evalIvarIpas(1., d0, ivar, jvar, &modelPart._calcmode);
        double delta = vexp - vtheo;
        total += delta * delta;
      }
  }
  _printResult("Cost Function (VMap Fit)", modelPart, total);
  
  return total;
}

int ModelOptimVMap::loadEnvironment(const DbGrid* dbmap, bool flagGoulard, bool verbose)
{
  _modelPart._verbose = verbose;
  _modelPart._optvar.setFlagGoulardUsed(flagGoulard);
  _vmapPart._dbmap    = dbmap;

  // Get internal dimension
  if (_getDimensions()) return 1;

  // Constitute the list of parameters
  if (_buildModelParamList()) return 1;

  // Check consistency
  if (!_checkConsistency()) return 1;

  // Instantiate Goulard algorithm (optional)
  if (flagGoulard)
  {
    _optGoulard = ModelOptimSillsVMap(_modelPart._model, _constraints, _mauto,
                                      _modelPart._optvar);
    _optGoulard.loadEnvironment(dbmap, verbose);
  }

  return 0;
}

int ModelOptimVMap::fit(const DbGrid* dbmap, bool flagGoulard, bool verbose)
{
  // Load the environment
  if (loadEnvironment(dbmap, flagGoulard, verbose)) return 1;

  // Perform the optimization
  AlgorithmVMap algorithm {_modelPart, _vmapPart, _optGoulard};
  _performOptimization(evalCost, &algorithm);
    
  return 0;
}

int ModelOptimVMap::_getDimensions()
{
  if (_vmapPart._dbmap == nullptr)
  {
    messerr("You must have defined 'dbmap' beforehand");
    return 1;
  }
  int nbexp  = 0;
  int npadir = 0;
  int nvar   = _vmapPart._dbmap->getNLoc(ELoc::Z);
  int nech   = _vmapPart._dbmap->getNSample();
  int nvs2   = nvar * (nvar + 1) / 2;

  /* Calculate the total number of lags */

  for (int iech = 0; iech < nech; iech++)
  {
    int ndef = 0;
    for (int ijvar = 0; ijvar < nvs2; ijvar++)
      if (!FFFF(_vmapPart._dbmap->getZVariable(iech, ijvar))) ndef++;
    nbexp += ndef;
    if (ndef > 0) npadir++;
  }

  /* Setting the return arguments */

  _vmapPart._npadir = npadir;

  if (nbexp <= 0)
  {
    messerr("No active experimental variogram map samples");
    return 1;
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Fill the array of pointers on the experimental conditions
 **
 *****************************************************************************/
void ModelOptimVMap::_computeFromVMap()
{
  const DbGrid* dbmap = _vmapPart._dbmap;
  int nvar            = dbmap->getNLoc(ELoc::Z);
  int nech            = dbmap->getNSample();
  int nvs2            = nvar * (nvar + 1) / 2;
  int npadir          = _vmapPart._npadir;
  dbmap->rankToIndice(nech / 2, _vmapPart._indg1);
  VectorDouble _wt; // TODO: fictitious dimensions to let compiling correct
  VectorDouble _gg;

  /* Load the Experimental conditions structure */

  int ipadir = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    dbmap->rankToIndice(iech, _vmapPart._indg2);
    double dist = distance_intra(dbmap, nech / 2, iech, NULL);
    double wgt  = (dist > 0) ? 1. / dist : 0.;

    /* Check samples containing only undefined values */

    int ntest = 0;
    for (int ijvar = 0; ijvar < nvs2; ijvar++)
      if (!FFFF(dbmap->getZVariable(iech, ijvar))) ntest++;
    if (ntest <= 0) continue;

    for (int ijvar = 0; ijvar < nvs2; ijvar++)
    {
      _WT(ijvar, ipadir) = 0.;
      _GG(ijvar, ipadir) = 0.;

      double value = dbmap->getZVariable(iech, ijvar);
      if (FFFF(value)) continue;

      _WT(ijvar, ipadir) = wgt;
      _GG(ijvar, ipadir) = value;
    }
    ipadir++;
  }
}
