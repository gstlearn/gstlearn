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

#include "Estimation/ALikelihood.hpp"
#include "Db/Db.hpp"
#include "Model/ModelGeneric.hpp"
#include "Basic/VectorHelper.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "LinearOp/CholeskyDense.hpp"
ALikelihood::ALikelihood(const ModelGeneric* model,
                         const Db* db)
  : AModelOptimNew(model)
  , _db(db)
{
}

ALikelihood::ALikelihood(const ALikelihood& r)
  : AModelOptimNew(r)
  , _db(r._db)
  , _Y(r._Y)
  , _X(r._X)
  , _beta(r._beta)
  , _Cm1X(r._Cm1X)
  , _Cm1Y(r._Cm1Y) {};

ALikelihood& ALikelihood::operator=(const ALikelihood& r)
{
  if (this != &r)
  {
    AModelOptimNew::operator=(r);
    _db   = r._db;
    _Y    = r._Y;
    _X    = r._X;
    _beta = r._beta;
    _Cm1X = r._Cm1X;
    _Cm1Y = r._Cm1Y;
  }
  return *this;
}

ALikelihood::~ALikelihood()
{
}

void ALikelihood::init(bool verbose)
{
  int nvar = _db->getNLoc(ELoc::Z);
  if (nvar < 1)
  {
    messerr("The 'db' should have at least one variable defined");
  }

  // Establish the vector of multivariate data
  int nDrift = _model->getNDriftEquation();
  if (nDrift > 0)
    _Y = _db->getColumnsByLocator(ELoc::Z, true, true);
  else
    _Y = _db->getColumnsByLocator(ELoc::Z, true, true, _model->getMeans());

  int size = (int)_Y.size();
  if (verbose)
  {
    message("Likelihood calculation:\n");
    message("- Number of active samples     = %d\n", _db->getNSample(true));
    message("- Number of variables          = %d\n", nvar);
    message("- Length of Information Vector = %d\n", size);
    if (nDrift > 0)
      message("- Number of drift conditions = %d\n", _model->getNDriftEquation());
    else
      VH::dump("Constant Mean(s)", _model->getMeans());
  }

  // If Drift function is present, evaluate the optimal Drift coefficients
  if (nDrift > 0)
  {
    // Extract the matrix of drifts at samples X
    _X = _model->evalDriftMat(_db);

    _beta.resize(nDrift);
  }
  _init();
}

double ALikelihood::computeCost(bool verbose)
{
  _updateModel(verbose);
  int nDrift = _model->getNDriftEquation();

  if (nDrift > 0)
  {
    // Calculate t(L-1) %*% D-1 %*% L-1 applied to X (L and D from Vecchia)
    _computeCm1X();

    // Calculate XtCm1X = Xt * Cm1 * X
    MatrixSymmetric* XtCm1X =
      MatrixFactory::prodMatMat<MatrixSymmetric>(&_X, &_Cm1X, true, false);

    // Construct ZtCm1X = Zt * Cm1 * X and perform its Cholesky decomposition
    VectorDouble ZtCm1X = _Cm1X.prodVecMat(_Y);
    CholeskyDense XtCm1XChol(XtCm1X);
    if (!XtCm1XChol.isReady())
    {
      messerr("Cholesky decomposition of XtCm1X matrix failed");
      delete XtCm1X;
      return TEST;
    }

    // Calculate beta = (XtCm1X)-1 * ZtCm1X
    if (XtCm1XChol.solve(ZtCm1X, _beta))
    {
      messerr("Error when calculating Likelihood");
      delete XtCm1X;
      return TEST;
    }
    // model->setBetaHat(beta);
    delete XtCm1X;

    if (verbose)
    {
      VH::dump("Optimal Drift coefficients = ", _beta);
    }

    // Center the data by the optimal drift: Y = Y - beta * X
    VH::subtractInPlace(_Y, _X.prodMatVec(_beta));
  }

  // Calculate t(L-1) %*% D-1 %*% L-1 applied to Y (L and D from Vecchia)

  _computeCm1Z();

  // Calculate the log-determinant

  double logdet = _computeLogDet();
  // Calculate quad = Zt * Cm1Z
  double quad = VH::innerProduct(_Y, _Cm1Y);

  // Derive the log-likelihood
  int size       = (int)_Y.size();
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