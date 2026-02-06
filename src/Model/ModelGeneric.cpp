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
#include "Basic/ListParams.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/CovAniso.hpp"
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
#include "Transform/ATransform.hpp"
#include <memory>

namespace gstlrn
{
ModelGeneric::ModelGeneric(const CovContext& ctxt)
  : ASerializable()
  , _cova(nullptr)
  , _driftList(nullptr)
  , _ctxt(ctxt)
  , _transform(nullptr)
{
}

ModelGeneric::ModelGeneric(const ModelGeneric& r)
  : ASerializable(r)
{
  _cova      = (r._cova != nullptr) ? std::dynamic_pointer_cast<ACov>(r._cova->cloneShared()) : nullptr;
  _driftList = (r._driftList != nullptr) ? r._driftList->clone() : nullptr;
  _ctxt      = r._ctxt;
  _transform = (r._transform != nullptr) ? std::dynamic_pointer_cast<ATransform>(r._transform->cloneShared()) : nullptr;
}

ModelGeneric& ModelGeneric::operator=(const ModelGeneric& r)
{
  if (this != &r)
  {
    ASerializable::operator=(r);
    _cova      = (r._cova != nullptr) ? std::dynamic_pointer_cast<ACov>(r._cova->cloneShared()) : nullptr;
    _driftList = (r._driftList != nullptr) ? r._driftList->clone() : nullptr;
    _ctxt      = r._ctxt;
    _transform = (r._transform != nullptr) ? std::dynamic_pointer_cast<ATransform>(r._transform->cloneShared()) : nullptr;
  }
  return *this;
}

ModelGeneric::~ModelGeneric()
{
  _clear();
  _transform = nullptr;
}

void ModelGeneric::_clear()
{
  _cova = nullptr;
  delete _driftList;
  _driftList = nullptr;
}

void ModelGeneric::_create()
{
  // Manage the DriftList part
  delete _driftList;
  _driftList = new DriftList(_ctxt);

  // Manage the ATransform part
  _transform = nullptr;
}

ModelGeneric* ModelGeneric::create(const CovContext& ctxt)
{
  return new ModelGeneric(ctxt);
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
  auto like = Likelihood(this, db, false);
  like.initLikelihood(verbose);
  return like.computeLogLikelihood(verbose);
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
  auto* drift_loc = dynamic_cast<ADrift*>(drift->clone());
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
  _gradCovFuncs.clear();

  // Add Covariance parameters
  if (_cova != nullptr)
  {
    _cova->appendParams(*listParams, &_gradCovFuncs);
  }

  // Add Drift parameters
  if (_driftList != nullptr)
  {
    _driftList->appendParams(*listParams);
  }

  if (_transform != nullptr)
  {
    _transform->appendParams(*listParams);
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

  if (_transform != nullptr)
  {
    _transform->updateTransform();
  }
}

void ModelGeneric::initParams(const MatrixSymmetric& vars, double href, double min, double max)
{
  // Initialize the parameters in the Covariance
  if (_cova != nullptr)
  {
    _cova->initParams(vars, href);
  }

  // Initialize the parameters in the DriftList
  if (_driftList != nullptr)
  {
    DriftList::initParams(vars, href);
  }
  if (_transform != nullptr)
  {
    _transform->initParams(min, max);
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

  if (amopt == nullptr)
  {
    messerr("No Optimizer could be created");
    return;
  }

  amopt->setVerbose(verbose, trace);
  amopt->resetIter();
  amopt->run();
  delete amopt;

  // Cancel the structure possibly used for Goulard (to be improved)
  auto* mcv = dynamic_cast<ModelCovList*>(this);
  if (mcv != nullptr)
    mcv->deleteFitSills();
}

void ModelGeneric::setTransform(const ATransform* transform)
{
  if (transform == nullptr)
  {
    _transform = nullptr;
    return;
  }
  _transform = std::dynamic_pointer_cast<ATransform>(transform->cloneShared());
}

#ifdef HDF5
bool ModelGeneric::deserializeH5(H5::Group& grp)
{
  auto modelG = SerializeHDF5::getGroup(grp, "Model");
  if (!modelG) return false;

  bool ret       = true;
  Id varioNumber = 0;
  ret            = ret && SerializeHDF5::readValue(*modelG, "Version Number", varioNumber);
  if (varioNumber != 2)
  {
    messerr("ModelGeneric::deserializeH5: Unsupported version %d", varioNumber);
    return false;
  }

  // Deserialize the CovContext characteristics
  Id ndim      = 0;
  Id nvar      = 0;
  double field = 0.;

  ret = ret && SerializeHDF5::readValue(*modelG, "NDim", ndim);
  ret = ret && SerializeHDF5::readValue(*modelG, "NVar", nvar);
  ret = ret && SerializeHDF5::readValue(*modelG, "Field", field);
  if (!ret) return ret;

  _ctxt = CovContext(nvar, ndim);
  _ctxt.setField(field);

  _clear();
  _create(); // Requires the context to exist

  // Process the Covariances
  ret = ret && _cova->deserializeH5(*modelG);

  // Process the drift part
  ret = ret && _driftList->deserializeH5(*modelG);

  // Process the covariance matrix
  VectorDouble covar0s;
  ret = ret && SerializeHDF5::readVec(*modelG, "VarCov0", covar0s);
  setCovar0s(covar0s);

  return ret;
}

bool ModelGeneric::serializeH5(H5::Group& grp) const
{
  bool ret = true;

  auto modelG    = grp.createGroup("Model");
  Id varioNumber = 2;
  ret            = ret && SerializeHDF5::writeValue(modelG, "Version Number", static_cast<Id>(varioNumber));

  // Serialize the CovContext characteristics
  ret = ret && SerializeHDF5::writeValue(modelG, "NDim", static_cast<Id>(getNDim()));
  ret = ret && SerializeHDF5::writeValue(modelG, "NVar", getNVar());
  ret = ret && SerializeHDF5::writeValue(modelG, "Field", getField());

  // Serialize the covariance part
  ret = ret && _cova->serializeH5(modelG);

  // Serialize the drift part
  ret = ret && _driftList->serializeH5(modelG);

  /// Serialize the variance-covariance at the origin
  ret = ret && SerializeHDF5::writeVec(modelG, "VarCov0", getCovar0());

  return ret;
}
#endif

} // namespace gstlrn