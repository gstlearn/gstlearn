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
#include "Model/ModelGeneric.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ListParams.hpp"
#include "Db/Db.hpp"
#include "Drifts/DriftFactory.hpp"
#include "Estimation/AModelOptim.hpp"
#include "Estimation/AModelOptimFactory.hpp"
#include "Estimation/Likelihood.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Model/AModelFitSills.hpp"
#include "Model/Model.hpp"
#include "Model/ModelOptimVMap.hpp"
#include "Model/ModelOptimVario.hpp"
#include <memory>

namespace gstlrn
{
ModelGeneric::ModelGeneric(const CovContext& ctxt)
  : _cova(nullptr)
  , _driftList(nullptr)
  , _ctxt(ctxt)
{
}

ModelGeneric::ModelGeneric(const ModelGeneric& r)
{
  _cova      = (r._cova != nullptr) ? std::dynamic_pointer_cast<ACov>(r._cova->cloneShared()) : nullptr;
  _driftList = (r._driftList != nullptr) ? r._driftList->clone() : nullptr;
  _ctxt      = r._ctxt;
}

ModelGeneric& ModelGeneric::operator=(const ModelGeneric& r)
{
  if (this != &r)
  {
    _cova      = (r._cova != nullptr) ? std::dynamic_pointer_cast<ACov>(r._cova->cloneShared()) : nullptr;
    _driftList = (r._driftList != nullptr) ? r._driftList->clone() : nullptr;
    _ctxt      = r._ctxt;
  }
  return *this;
}

ModelGeneric::~ModelGeneric()
{
  _cova = nullptr;
  delete _driftList;
  _driftList = nullptr;
}

void ModelGeneric::setField(double field)
{
  _ctxt.setField(field);
  setContext(_ctxt);
  copyCovContext(_ctxt);
}

void ModelGeneric::setContext(const CovContext& ctxt)
{
  _ctxt = ctxt;
  _cova->setContext(ctxt);
}

bool ModelGeneric::isValid() const
{
  return _isValid();
}

bool ModelGeneric::_isValid() const
{
  return true;
}

/**
 * Compute the log-likelihood (based on covariance)
 *
 * @param db  Db structure where variable are loaded from
 * @param verbose Verbose flag
 *
 * @remarks The calculation considers all the active samples.
 * @remarks It can work in multivariate case with or without drift conditions (linked or not)
 * @remarks The algorithm is stopped (with a message) in the heterotopic case
 * // TODO; improve for heterotopic case
 */
double ModelGeneric::computeLogLikelihood(const Db* db, bool verbose)
{
  auto* like = Likelihood::createForOptim(this, db);
  like->init(verbose);
  return like->computeLogLikelihood(verbose);
}

/**
 * Add a list of Drifts. This operation cleans any previously stored drift function
 * @param driftlist List of Drifts to be added
 *
 * @remark This method deletes any pre-existing drift functions
 */
void ModelGeneric::setDriftList(const DriftList* driftlist)
{
  if (driftlist == nullptr) return;
  delete _driftList;
  _driftList = driftlist->clone();

  // Check that the DriftList has the same type of CovContext as the Model
  _driftList->copyCovContext(_ctxt);
}

void ModelGeneric::setCov(const ACov* cova)
{
  if (cova == nullptr) return;

  // Set the Context of ModelGeneric (cross_check with DriftList)
  if (_driftList != nullptr)
  {
    if (!_driftList->getContext().isEqual(cova->getContext()))
    {
      messerr("Cova and DriftList do not share the same CovContext");
      messerr("Operation cancelled");
      return;
    }
  }
  _cova = (std::dynamic_pointer_cast<ACov>)(cova->cloneShared());

  _ctxt = cova->getContext();
}

/**
 * Define the list of drift functions for:
 * - a given degree of the IRF
 * - a given number of external drifts
 * @param order Order of the IRF
 * @param nfex  Number of External Drifts
 *
 * @remark This method deletes any pre-existing drift functions and replaces them by the new definition
 * @remark This replacement is performed accounting for information stored in 'model', such as:
 * - the space dimension
 * - the number of variables
 */
void ModelGeneric::setDriftIRF(Id order, Id nfex)
{
  delete _driftList;
  _driftList = DriftFactory::createDriftListFromIRF(order, nfex, _ctxt);
}

void ModelGeneric::addDrift(const ADrift* drift)
{
  if (drift == nullptr) return;
  if (_driftList == nullptr) _driftList = new DriftList(_ctxt);
  ADrift* drift_loc = dynamic_cast<ADrift*>(drift->clone());
  _driftList->addDrift(drift_loc);

  // Check that the DriftList has the same type of CovContext as the Model
  _driftList->copyCovContext(_ctxt);
}

void ModelGeneric::setDrifts(const VectorString& driftSymbols)
{
  if (_driftList == nullptr)
    _driftList = new DriftList();
  else
    delAllDrifts();

  for (Id i = 0; i < static_cast<Id>(driftSymbols.size()); i++)
  {
    ADrift* drift = DriftFactory::createDriftBySymbol(driftSymbols[i]);
    addDrift(drift);
  }
}

static MatrixDense _transformF(const MatrixDense& F1, Id type, Id idx)
{
  MatrixDense F1loc;
  switch (type)
  {
    case 1:
    case 2:
      F1loc = F1;
      break;
    case 3:
      F1loc = F1;
      F1loc.fill(0.);
      break;
    case 4:
      F1loc = F1;
      F1loc.fill(0.);
      for (Id i = 0; i < F1.getNRows(); i++)
        F1loc.setValue(i, idx, 1.);
      break;
  }
  return (F1loc);
}

Id computeCovMatSVCLHSInPlace(MatrixSymmetric& cov,
                              const MatrixSymmetric& Sigma,
                              const MatrixDense& F1,
                              Id type,
                              Id idx)
{
  MatrixDense F1loc = _transformF(F1, type, idx);
  auto nech         = F1.getNRows();
  auto nvar         = Sigma.getNRows() / nech;
  cov.resize(nech, nech);

  for (Id iech = 0; iech < nech; iech++)
  {
    for (Id jech = 0; jech < nech; jech++)
    {
      if (iech > jech) continue;
      double value = 0.;
      for (Id lvar = 0; lvar < nvar; lvar++)
      {
        for (Id pvar = 0; pvar < nvar; pvar++)
        {
          Id shifti = lvar * nech;
          Id shiftj = pvar * nech;
          value += Sigma.getValue(shifti + iech, shiftj + jech) *
                   F1loc.getValue(iech, lvar) *
                   F1loc.getValue(jech, pvar);
        }
      }
      cov.setValue(iech, jech, value);
    }
  }
  return 0;
}

Id computeCovMatSVCRHSInPlace(MatrixDense& cov,
                              const MatrixSymmetric& Sigma,
                              const MatrixDense& F1,
                              const MatrixDense& F2,
                              Id type1,
                              Id idx1,
                              Id type2,
                              Id idx2)
{
  MatrixDense F1loc = _transformF(F1, type1, idx1);
  MatrixDense F2loc = _transformF(F2, type2, idx2);
  auto nech1        = F1.getNRows();
  auto nech2        = F2.getNRows();
  auto nvar         = Sigma.getNCols();
  cov.resize(nech1, nech2);

  for (Id iech = 0; iech < nech1; iech++)
  {
    for (Id jech = 0; jech < nech2; jech++)
    {
      double value = 0.;
      for (Id lvar = 0; lvar < nvar; lvar++)
      {
        for (Id pvar = 0; pvar < nvar; pvar++)
        {
          Id shifti = lvar * nech1;
          Id shiftj = pvar * nech2;
          value += Sigma.getValue(shifti + iech, shiftj + jech) *
                   F1loc.getValue(iech, lvar) *
                   F2loc.getValue(jech, pvar);
        }
      }
      cov.setValue(iech, jech, value);
    }
  }
  return 0;
}

Id computeDriftMatSVCRHSInPlace(MatrixDense& mat,
                                const MatrixDense& F,
                                Id type,
                                Id idx,
                                bool flagCenteredFactors)
{
  if (flagCenteredFactors)
  {
    mat.resize(1, 1);
    switch (type)
    {
      case 1:
      case 3:
        mat.setRowToConstant(0, 1.);
        break;
      case 2:
        mat.setRowToConstant(0, 0.);
        break;
      case 4:
        mat.setValue(0, 0, (idx == 0) ? 1. : 0.);
        break;
    }
  }
  else
  {
    mat.resize(1, F.getNCols());
    switch (type)
    {
      case 1:
      case 3:
        mat.setRow(0, F.getRow(0));
        break;

      case 2:
        mat.setRowToConstant(0, 0.);
        break;

      case 4:
        mat.setRowToConstant(0, 0.);
        mat.setValue(0, idx, 1.);
        break;
    }
  }
  return 0;
}

std::shared_ptr<ListParams> ModelGeneric::generateListParams() const
{
  auto listParams = std::make_shared<ListParams>();
  _gradFuncs.clear();

  // Add Covariance parameters
  if (_cova != nullptr)
  {
    _cova->appendParams(*listParams, &_gradFuncs);
  }

  // Add Drift parameters
  if (_driftList != nullptr)
  {
    _driftList->appendParams(*listParams);
  }
  listParams->updateDispatch();

  return listParams;
}

ListParams* ModelGeneric::createListParams(std::shared_ptr<ListParams>& lp)
{
  return lp.get();
}

void ModelGeneric::updateModel()
{
  // Update the Covariance
  if (_cova != nullptr)
  {
    _cova->updateCov();
  }

  // Update the DriftList
  if (_driftList != nullptr)
  {
    _driftList->updateDriftList();
  }
}

void ModelGeneric::initParams(const MatrixSymmetric& vars, double href)
{
  // Initialize the parameters in the Covariance
  if (_cova != nullptr)
  {
    _cova->initParams(vars, href);
  }

  // Initialize the parameters in the DriftList
  if (_driftList != nullptr)
  {
    gstlrn::DriftList::initParams(vars, href);
  }
}

void ModelGeneric::fitNew(const Db* db,
                          Vario* vario,
                          const DbGrid* dbmap,
                          Constraints* constraints,
                          const ModelOptimParam& mop,
                          Id nb_neighVecchia,
                          bool verbose,
                          bool trace,
                          bool reml)
{
  AModelOptim* amopt = AModelOptimFactory::create(this, db, vario, dbmap,
                                                  constraints, mop,
                                                  nb_neighVecchia,
                                                  reml);
  amopt->setVerbose(verbose, trace);
  amopt->resetIter();
  amopt->run();
  delete amopt;

  // Cancel the structure possibly used for Goulard (to be improved)
  ModelCovList* mcv = dynamic_cast<ModelCovList*>(this);
  if (mcv != nullptr)
    mcv->deleteFitSills();
}
} // namespace gstlrn