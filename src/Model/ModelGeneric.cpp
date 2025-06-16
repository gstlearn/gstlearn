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
#include "Estimation/AModelOptim.hpp"
#include "Model/AModelFitSills.hpp"
#include "Model/Model.hpp"
#include "Model/ModelOptimVario.hpp"
#include "Model/ModelOptimVMap.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Db/Db.hpp"
#include "Estimation/AModelOptimFactory.hpp"
#include "Drifts/DriftFactory.hpp"
#include "Estimation/Likelihood.hpp"
#include "geoslib_define.h"
#include <memory>

ModelGeneric::ModelGeneric(const CovContext& ctxt)
  : _cova(nullptr)
  , _driftList(nullptr)
  , _ctxt(ctxt)
{
}

ModelGeneric::ModelGeneric(const ModelGeneric& r)
{
  _cova      = (r._cova != nullptr) ? (ACov*)r._cova->clone() : nullptr;
  _driftList = (r._driftList != nullptr) ? r._driftList->clone() : nullptr;
  _ctxt      = r._ctxt;
}

ModelGeneric& ModelGeneric::operator=(const ModelGeneric& r)
{
  if (this != &r)
  {
    _cova      = (r._cova != nullptr) ? (ACov*)r._cova->clone() : nullptr;
    _driftList = (r._driftList != nullptr) ? r._driftList->clone() : nullptr;
    _ctxt      = r._ctxt;
  }
  return *this;
}

ModelGeneric::~ModelGeneric()
{
  delete _cova;
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
  return like->computeCost();
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

void ModelGeneric::setCov(ACov* cova)
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
  delete _cova;
  _cova = (ACov*)cova->clone();

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
void ModelGeneric::setDriftIRF(int order, int nfex)
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

  for (int i = 0; i < (int)driftSymbols.size(); i++)
  {
    ADrift* drift = DriftFactory::createDriftBySymbol(driftSymbols[i]);
    addDrift(drift);
  }
}

static MatrixDense _transformF(const MatrixDense& F1, int type, int idx)
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
      for (int i = 0; i < F1.getNRows(); i++)
        F1loc.setValue(i, idx, 1.);
      break;
  }
  return (F1loc);
}

int computeCovMatSVCLHSInPlace(MatrixSymmetric& cov,
                               const MatrixSymmetric& Sigma,
                               const MatrixDense& F1,
                               int type,
                               int idx)
{
  MatrixDense F1loc = _transformF(F1, type, idx);
  int nech          = F1.getNRows();
  int nvar          = Sigma.getNRows() / nech;
  cov.resize(nech, nech);

  for (int iech = 0; iech < nech; iech++)
  {
    for (int jech = 0; jech < nech; jech++)
    {
      if (iech > jech) continue;
      double value = 0.;
      for (int lvar = 0; lvar < nvar; lvar++)
      {
        for (int pvar = 0; pvar < nvar; pvar++)
        {
          int shifti = lvar * nech;
          int shiftj = pvar * nech;
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

int computeCovMatSVCRHSInPlace(MatrixDense& cov,
                               const MatrixSymmetric& Sigma,
                               const MatrixDense& F1,
                               const MatrixDense& F2,
                               int type1,
                               int idx1,
                               int type2,
                               int idx2)
{
  MatrixDense F1loc = _transformF(F1, type1, idx1);
  MatrixDense F2loc = _transformF(F2, type2, idx2);
  int nech1         = F1.getNRows();
  int nech2         = F2.getNRows();
  int nvar          = Sigma.getNCols();
  cov.resize(nech1, nech2);

  for (int iech = 0; iech < nech1; iech++)
  {
    for (int jech = 0; jech < nech2; jech++)
    {
      double value = 0.;
      for (int lvar = 0; lvar < nvar; lvar++)
      {
        for (int pvar = 0; pvar < nvar; pvar++)
        {
          int shifti = lvar * nech1;
          int shiftj = pvar * nech2;
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

int computeDriftMatSVCRHSInPlace(MatrixDense& mat,
                                 const MatrixDense& F,
                                 int type,
                                 int idx,
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

  // Add Covariance parameters
  if (_cova != nullptr)
  {
    _cova->appendParams(*listParams);
  }

  // Add Drift parameters
  if (_driftList != nullptr)
  {
    _driftList->appendParams(*listParams);
  }

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

void ModelGeneric::initParams()
{
  // Initialize the parameters in the Covariance
  if (_cova != nullptr)
  {
    _cova->initParams();
  }

  // Initialize the parameters in the DriftList
  if (_driftList != nullptr)
  {
    _driftList->initParams();
  }
}

void ModelGeneric::fitNew(const Db* db,
                          Vario* vario,
                          const DbGrid* dbmap,
                          Constraints* constraints,
                          const ModelOptimParam& mop,
                          int nb_neighVecchia,
                          bool verbose,
                          bool trace)
{
  AModelOptim* amopt = AModelOptimFactory::create(this, db, vario, dbmap,
                                                  constraints, mop,
                                                  nb_neighVecchia);
  amopt->setVerbose(verbose, trace);
  amopt->resetIter();
  amopt->run();
  delete amopt;

  // Cancel the structure possibly used for Goulard (to be improved)
  ModelCovList* mcv = dynamic_cast<ModelCovList*>(this);
  if (mcv != nullptr)
  {
    mcv->deleteFitSills();
  }
}
