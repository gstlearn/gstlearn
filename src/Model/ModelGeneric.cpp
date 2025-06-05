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
#include "Basic/Optim.hpp"
#include "Estimation/Vecchia.hpp"
#include "Model/Model.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Db/Db.hpp"
#include "Basic/VectorHelper.hpp"
#include "Drifts/DriftFactory.hpp"
#include <memory>
#include <nlopt.h>

ModelGeneric::ModelGeneric(const CovContext& ctxt)
  : _cova(nullptr)
  , _driftList(nullptr)
  , _ctxt(ctxt)
{
  _driftList = new DriftList(_ctxt);
}

ModelGeneric::ModelGeneric(const ModelGeneric& r)
{
  if (r._cova != nullptr) _cova = (ACov*)r._cova->clone();
  if (r._driftList != nullptr) _driftList = r._driftList->clone();

  _ctxt = r._ctxt;
}

ModelGeneric& ModelGeneric::operator=(const ModelGeneric& r)
{
  if (this != &r)
  {
    _cova      = (ACov*)r._cova->clone();
    _driftList = r._driftList->clone();
    _ctxt      = r._ctxt;
  }
  return *this;
}

ModelGeneric::~ModelGeneric()
{
  delete _cova;
  delete _driftList;
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
  int nvar = db->getNLoc(ELoc::Z);
  if (nvar < 1)
  {
    messerr("The 'db' should have at least one variable defined");
    return TEST;
  }
  int nDrift = getNDriftEquation();

  // Calculate the covariance matrix C and perform its Cholesky decomposition
  MatrixSymmetric cov = evalCovMatSym(db);
  CholeskyDense covChol(&cov);
  if (!covChol.isReady())
  {
    messerr("Cholesky decomposition of Covariance matrix failed");
    return TEST;
  }

  // Establish the vector of multivariate data
  VectorDouble Z;
  if (nDrift > 0)
    Z = db->getColumnsByLocator(ELoc::Z, true, true);
  else
    Z = db->getColumnsByLocator(ELoc::Z, true, true, getMeans());

  int size = (int)Z.size();
  if (verbose)
  {
    message("Likelihood calculation:\n");
    message("- Number of active samples     = %d\n", db->getNSample(true));
    message("- Number of variables          = %d\n", nvar);
    message("- Length of Information Vector = %d\n", size);
    if (nDrift > 0)
      message("- Number of drift conditions = %d\n", getNDriftEquation());
    else
      VH::dump("Constant Mean(s)", getMeans());
  }

  // If Drift functions are present, evaluate the optimal Drift coefficients
  if (nDrift > 0)
  {
    // Extract the matrix of drifts at samples X
    MatrixDense X = evalDriftMat(db);

    // Calculate Cm1X = Cm1 * X
    MatrixDense Cm1X;
    if (covChol.solveMatrix(X, Cm1X))
    {
      messerr("Problem when solving a Linear System after Cholesky decomposition");
      return TEST;
    }

    // Calculate XtCm1X = Xt * Cm1 * X
    MatrixSymmetric* XtCm1X =
      MatrixFactory::prodMatMat<MatrixSymmetric>(&X, &Cm1X, true, false);

    // Construct ZtCm1X = Zt * Cm1 * X and perform its Cholesky decomposition
    VectorDouble ZtCm1X = Cm1X.prodVecMat(Z);
    CholeskyDense XtCm1XChol(XtCm1X);
    if (!XtCm1XChol.isReady())
    {
      messerr("Cholesky decomposition of XtCm1X matrix failed");
      delete XtCm1X;
      return TEST;
    }

    // Calculate beta = (XtCm1X)-1 * ZtCm1X
    VectorDouble beta(nDrift);
    if (XtCm1XChol.solve(ZtCm1X, beta))
    {
      messerr("Error when calculating Likelihood");
      delete XtCm1X;
      return TEST;
    }
    setBetaHat(beta);
    delete XtCm1X;

    if (verbose)
    {
      VH::dump("Optimal Drift coefficients = ", beta);
    }

    // Center the data by the optimal drift: Z = Z - beta * X
    VH::subtractInPlace(Z, X.prodMatVec(beta));
  }

  // Calculate Cm1Z = Cm1 * Z
  VectorDouble Cm1Z(Z.size());
  if (covChol.solve(Z, Cm1Z))
  {
    messerr("Error when calculating Cm1Z");
    return TEST;
  }

  // Calculate the log-determinant
  double logdet = covChol.computeLogDeterminant();

  // Calculate quad = Zt * Cm1Z
  double quad = VH::innerProduct(Z, Cm1Z);

  // Derive the log-likelihood
  double loglike = -0.5 * (logdet + quad + size * log(2. * GV_PI));

  // Optional printout
  if (verbose)
  {
    message("Log-Determinant = %lf\n", logdet);
    message("Quadratic term  = %lf\n", quad);
    message("Log-likelihood  = %lf\n", loglike);
  }
  return loglike;
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
void ModelGeneric::fitLikelihood(const Db* db, bool useVecchia, bool verbose)
{
  auto params = generateListParams();
  initParams();
  std::vector<double> x    = params->getValues();
  std::vector<double> xmin = params->getMinValues();
  std::vector<double> xmax = params->getMaxValues();
  updateModel();
  if (verbose)
  {
    params->display();
    _cova->display();
  }

  Optim opt(NLOPT_LN_NELDERMEAD, x.size());

  AModelOptimNew* amopt = nullptr;

  if (useVecchia)
  {
    int nbneigh = std::min(30, db->getNSample(true));
    amopt = Vecchia::createForOptim(this, db, nbneigh);
  }
  auto func = [amopt,db, params, useVecchia, verbose, this](const std::vector<double>& x) -> double
  {
    static int iter = 1;
    params->setValues(x);
    this->updateModel();

    double result;
    if (!useVecchia)
      result = computeLogLikelihood(db);
    else
    {
      result      = amopt->computeCost(false); 
    }
    
    if (verbose)
      message("Iteration %3d - Cost Function (Likelihood) = %lf\n", iter++, result);
    return -result;
  };

  opt.setObjective(func);
  opt.setLowerBounds(xmin);
  opt.setUpperBounds(xmax);
  opt.setXtolRel(EPSILON6);

  opt.optimize(x);

  delete amopt;
}
