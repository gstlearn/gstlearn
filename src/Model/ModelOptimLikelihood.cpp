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
#include "Stats/Classical.hpp"
#include "API/SPDE.hpp"

#include <nlopt.h>

typedef struct
{
  // Part of the structure dedicated to the Model
  ModelOptim::Model_Part& _modelPart;

  // Part relative to the Experimental variograms
  ModelOptimLikelihood::Db_Part& _dbPart;

} AlgorithmLikelihood;

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

int ModelOptimLikelihood::fit(Db* db, Model* model, bool flagSPDE, bool verbose)
{
  _modelPart._model   = model;
  _modelPart._verbose = verbose;
  _dbPart._db         = db;
  _dbPart._flagSPDE   = flagSPDE;

  // Constitute the list of parameters
  if (_buildModelParamList()) return 1;

  // Check consistency
  if (! _checkConsistency()) return 1;

  // Define the optimization criterion
  int npar = _getParamNumber();
  nlopt_opt opt = nlopt_create(NLOPT_LN_NELDERMEAD, npar);
  nlopt_set_lower_bounds(opt, _modelPart._tablow.data());
  nlopt_set_upper_bounds(opt, _modelPart._tabupp.data());
  nlopt_srand(12345);
  nlopt_set_ftol_rel(opt, 1e-6);

  // Update the initial optimization values (due to variogram)
  updateModelParamList(db->getExtensionDiagonal(), dbVarianceMatrix(db));

  // Define the cost function
  AlgorithmLikelihood algorithm {_modelPart, _dbPart};
  nlopt_set_min_objective(opt, evalCost, &algorithm);

  // Perform the optimization (store the minimized value in 'minf')
  if (_modelPart._verbose) mestitle(1, "Model Fitting using Likelihood");
  double minf;
  nlopt_optimize(opt, _modelPart._tabval.data(), &minf);
  nlopt_destroy(opt);
  return 0;
}

double ModelOptimLikelihood::evalCost(unsigned int nparams,
                                      const double* current,
                                      double* /*grad*/,
                                      void* my_func_data)
{
  DECLARE_UNUSED(nparams);
  AlgorithmLikelihood* algorithm = static_cast<AlgorithmLikelihood*>(my_func_data);
  if (algorithm == nullptr) return TEST;
  Model_Part& modelPart = algorithm->_modelPart;
  Db_Part& dbPart       = algorithm->_dbPart;

  // Update the Model
  _patchModel(modelPart, current);

  // Evaluate the Cost function
  double total;
  if (dbPart._flagSPDE)
    total = -logLikelihoodSPDE(dbPart._db, modelPart._model);
  else
    total = -modelPart._model->computeLogLikelihood(dbPart._db);

  _printResult("Cost Function (Likelihood)", modelPart, total);

  return total;
}
