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

#include "Model/Model.hpp"
#include "Db/Db.hpp"
#include "Stats/Classical.hpp"
#include "API/SPDE.hpp"

ModelOptimLikelihood::ModelOptimLikelihood(Model* model)
  : AModelOptim(model)
  , _flagSPDE(false)
  , _db(nullptr)
{
}

ModelOptimLikelihood::ModelOptimLikelihood(const ModelOptimLikelihood& m)
  : AModelOptim(m) 
  , _flagSPDE(m._flagSPDE)
  , _db(m._db)
{
}

ModelOptimLikelihood& ModelOptimLikelihood::operator=(const ModelOptimLikelihood& m)
{
  if (this != &m)
  {
    _flagSPDE = m._flagSPDE;
    _db = m._db;
  }
  return (*this);
}

ModelOptimLikelihood::~ModelOptimLikelihood()
{
}

bool ModelOptimLikelihood::_checkConsistency()
{
  const Model* model = _modelPart._model;

  if ((int)model->getNDim() != _db->getNLoc(ELoc::X))
  {
    messerr("'_model'(%d) and 'db'(%d) should have same Space Dimension",
            model->getNDim(), _db->getNLoc(ELoc::X));
    return false;
  }
  if (model->getNVar() != _db->getNLoc(ELoc::Z))
  {
    messerr("'_model'(%d) and '_db'(%d) should have same number of Variables",
            model->getNVar(), _db->getNLoc(ELoc::Z));
    return false;
  }
  return true;
}

int ModelOptimLikelihood::loadEnvironment(Db* db, bool flagSPDE, bool verbose)
{
  _modelPart._verbose = verbose;
  _modelPart._optvar.setFlagGoulardUsed(false);
  _db         = db;
  _flagSPDE   = flagSPDE;

  // Constitute the list of parameters
  if (_buildModelParamList()) return 1;

  // Check consistency
  if (!_checkConsistency()) return 1;

  return 0;
}

// double ModelOptimLikelihood::evalCost(unsigned int nparams,
//                                       const double* current,
//                                       double* /*grad*/,
//                                       void* my_func_data)
// {
// Sauvegarder pour garder la trace de loglikelihood par SPDE.
//   // Evaluate the Cost function
//   double total;
//   if (dbPart._flagSPDE)
//     total = -logLikelihoodSPDEOld(dbPart._db, modelPart._model);
//   else
//     total = -modelPart._model->computeLogLikelihood(dbPart._db);
