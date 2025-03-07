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
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Db/Db.hpp"
#include "Basic/VectorHelper.hpp"
#include "Drifts/DriftFactory.hpp"

ModelGeneric::ModelGeneric(const CovContext &ctxt)
    : _cova(nullptr),
      _driftList(nullptr),
      _ctxt(ctxt)
{
}

ModelGeneric::ModelGeneric(const ModelGeneric& r)
{
  if (r._cova != nullptr) _cova = (ACov*)r._cova->clone();
  if (r._driftList != nullptr) _driftList = (DriftList*)r._driftList->clone();
  _ctxt = r._ctxt;
}

ModelGeneric& ModelGeneric::operator=(const ModelGeneric& r)
{
  if (this != &r)
  {
    _cova      = (ACov*)r._cova->clone();
    _driftList = (DriftList*)r._driftList->clone();
    _ctxt      = r._ctxt;
  }
  return *this;
}

ModelGeneric::~ModelGeneric()
{
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
  MatrixSquareSymmetric cov = evalCovMatSym(db);
  CholeskyDense covChol(&cov);
  if (! covChol.isReady())
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
    MatrixRectangular X = evalDriftMat(db);

    // Calculate Cm1X = Cm1 * X
    MatrixRectangular Cm1X;
    if (covChol.solveMatrix(X, Cm1X))
    {
      messerr("Problem when solving a Linear System after Cholesky decomposition");
      return TEST;
    }

    // Calculate XtCm1X = Xt * Cm1 * X
    MatrixSquareSymmetric* XtCm1X =
      MatrixFactory::prodMatMat<MatrixSquareSymmetric>(&X, &Cm1X, true, false);
   
    // Construct ZtCm1X = Zt * Cm1 * X and perform its Cholesky decomposition
    VectorDouble ZtCm1X = Cm1X.prodVecMat(Z);
    CholeskyDense XtCm1XChol(XtCm1X);
    if (! XtCm1XChol.isReady())
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
    message("Quadratic term = %lf\n", quad);
    message("Log-likelihood = %lf\n", loglike);
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
  delete _cova;
  _cova = cova;

  // Set the Context of ModelGeneric (cross_check with DriftList)
  if (_driftList != nullptr)
  {
    if (! _driftList->getContext().isEqual(cova->getContext()))
    {
      messerr("Cova and DriftList do not share the same CovContext");
      messerr("Operation cancelled");
      return;
    }
  }

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


void ModelGeneric::addDrift(const ADrift *drift)
{
  if (drift == nullptr) return;
  if (_driftList == nullptr) _driftList = new DriftList(_ctxt);
  ADrift* drift_loc = dynamic_cast<ADrift*>(drift->clone());
  _driftList->addDrift(drift_loc);

  // Check that the DriftList has the same type of CovContext as the Model
  _driftList->copyCovContext(_ctxt);
}

void ModelGeneric::setDrifts(const VectorString &driftSymbols)
{
  if (_driftList == nullptr)
    _driftList = new DriftList();
  else
    delAllDrifts();

  for (int i = 0; i < (int) driftSymbols.size(); i++)
  {
    ADrift *drift = DriftFactory::createDriftBySymbol(driftSymbols[i]);
    addDrift(drift);
  }
}

