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
#include "Model/ModelOptimLikelihood.hpp"

#include "Enum/ELoc.hpp"
#include "geoslib_define.h"

#include "Model/Model.hpp"
#include "Db/Db.hpp"

#include <nlopt.h>

ModelOptimLikelihood::ModelOptimLikelihood()
  : ModelOptim()
  , _dbPart()
{
}

ModelOptimLikelihood::ModelOptimLikelihood(const ModelOptimLikelihood& m)
  : ModelOptim(m)
  , _dbPart()
{
   _copyDbPart(m._dbPart);
}

void ModelOptimLikelihood::_copyDbPart(const Db_Part& dbPart)
{
  _dbPart._db = dbPart._db;
}

ModelOptimLikelihood& ModelOptimLikelihood::operator=(const ModelOptimLikelihood& m)
{
  if (this != &m)
  {
    _copyDbPart(m._dbPart);
  }
  return (*this);
}

ModelOptimLikelihood::~ModelOptimLikelihood()
{
}

bool ModelOptimLikelihood::_checkConsistency()
{
  const Model* model = _modelPart._model;
  const Db* db       = _dbPart._db;

  if (model->getDimensionNumber() != db->getLocatorNumber(ELoc::X))
  {
    messerr("'_model'(%d) and 'db'(%d) should have same Space Dimension",
            model->getDimensionNumber(), db->getLocatorNumber(ELoc::X));
    return false;
  }
  if (model->getVariableNumber() != db->getLocatorNumber(ELoc::Z))
  {
    messerr("'_model'(%d) and '_db'(%d) should have same number of Variables",
            model->getVariableNumber(), db->getLocatorNumber(ELoc::Z));
    return false;
  }
  return true;
}

int ModelOptimLikelihood::fit(const Db* db, Model* model, bool verbose)
{
  _modelPart._model   = model;
  _modelPart._verbose = verbose;
  _dbPart._db         = db;

  // Constitute the list of parameters
  if (_buildModelParamList()) return 1;

  // Check consistency
  if (! _checkConsistency()) return 1;

  // Perform the optimization
  int npar = _getParamNumber();
  nlopt_opt opt = nlopt_create(NLOPT_LN_NELDERMEAD, npar);

  // Define the bounds
  nlopt_set_lower_bounds(opt, _modelPart._tablow.data());
  nlopt_set_upper_bounds(opt, _modelPart._tabupp.data());

  // Define the cost function
  AlgorithmLikelihood algorithm{_modelPart, _dbPart};
  nlopt_set_min_objective(opt, evalCost, &algorithm);

  // Perform the optimization (store the minimized value in 'minf')
  if (_modelPart._verbose) mestitle(1, "Model Fitting from Likelihood");
  double minf;
  nlopt_optimize(opt, _modelPart._tabval.data(), &minf);
  return 0;
}

double ModelOptimLikelihood::evalCost(unsigned int nparams,
                                      const double* current,
                                      double* grad,
                                      void* my_func_data)
{
  DECLARE_UNUSED(nparams);
  DECLARE_UNUSED(grad);
  AlgorithmLikelihood* algorithm = (AlgorithmLikelihood*)(my_func_data);
  if (algorithm == nullptr) return TEST;
  Model_Part modelPart = algorithm->_modelPart;
  Db_Part dbPart       = algorithm->_dbPart;

  // Update the Model
  _patchModel(modelPart, current);

  // Evaluate the Cost function
  double total = -modelPart._model->computeLogLikelihood(dbPart._db);

  if (modelPart._verbose) message("Cost Function = %lf\n", total);

  return total;
}
