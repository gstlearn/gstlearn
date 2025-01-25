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

ModelOptimLikelihood::ModelOptimLikelihood(Model* model)
  : AModelOptim(model)
  , _dbPart()
{
}

ModelOptimLikelihood::ModelOptimLikelihood(const ModelOptimLikelihood& m)
  : AModelOptim(m)
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

  if (model->getDimensionNumber() != db->getNLoc(ELoc::X))
  {
    messerr("'_model'(%d) and 'db'(%d) should have same Space Dimension",
            model->getDimensionNumber(), db->getNLoc(ELoc::X));
    return false;
  }
  if (model->getNVar() != db->getNLoc(ELoc::Z))
  {
    messerr("'_model'(%d) and '_db'(%d) should have same number of Variables",
            model->getNVar(), db->getNLoc(ELoc::Z));
    return false;
  }
  return true;
}

int ModelOptimLikelihood::loadEnvironment(Db* db, bool flagSPDE, bool verbose)
{
  _modelPart._verbose = verbose;
  _modelPart._optvar.setFlagGoulardUsed(false);
  _dbPart._db         = db;
  _dbPart._flagSPDE   = flagSPDE;

  // Constitute the list of parameters
  if (_buildModelParamList()) return 1;

  // Check consistency
  if (!_checkConsistency()) return 1;

  return 0;
}

int ModelOptimLikelihood::fit(Db* db, bool flagSPDE, bool verbose)
{
  // Load the environment
  if (loadEnvironment(db, flagSPDE, verbose)) return 1;

  // Perform the optimization
  AlgorithmLikelihood algorithm {_modelPart, _dbPart};
  _performOptimization(evalCost, &algorithm, db->getExtensionDiagonal(),
                       dbVarianceMatrix(db));

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
