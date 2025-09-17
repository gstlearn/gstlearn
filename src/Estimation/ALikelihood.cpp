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
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Model/ModelGeneric.hpp"

namespace gstlrn
{
ALikelihood::ALikelihood(ModelGeneric* model,
                         const Db* db,
                         bool reml)
  : AModelOptim(model)
  , _db(db)
  , _reml(reml)
{
}

ALikelihood::ALikelihood(const ALikelihood& r)
  : AModelOptim(r)
  , _db(r._db)
  , _Y(r._Y)
  , _X(r._X)
  , _beta(r._beta)
  , _Cm1X(r._Cm1X)
  , _Cm1Y(r._Cm1Y)
  , _XtCm1X(r._XtCm1X)
  , _reml(r._reml) {};

ALikelihood& ALikelihood::operator=(const ALikelihood& r)
{
  if (this != &r)
  {
    AModelOptim::operator=(r);
    _db     = r._db;
    _Y      = r._Y;
    _X      = r._X;
    _beta   = r._beta;
    _Cm1X   = r._Cm1X;
    _Cm1Y   = r._Cm1Y;
    _XtCm1X = r._XtCm1X;
    _reml   = r._reml;
  }
  return *this;
}

ALikelihood::~ALikelihood()
{
}

void ALikelihood::init(bool verbose)
{
  Id nvar = _db->getNLoc(ELoc::Z);
  if (nvar < 1)
  {
    messerr("The 'db' should have at least one variable defined");
  }

  // Establish the vector of multivariate data
  Id nDrift = _model->getNDriftEquation();
  if (nDrift > 0)
    _Y = _db->getColumnsByLocator(ELoc::Z, true, true);
  else
    _Y = _db->getColumnsByLocator(ELoc::Z, true, true, _model->getMeans());

  Id size = static_cast<Id>(_Y.size());
  if (verbose)
  {
    message("Likelihood calculation:\n");
    message("- Number of active samples     = %d\n", _db->getNSample(true));
    message("- Number of variables          = %d\n", nvar);
    message("- Length of Information Vector = %d\n", size);
    if (nDrift > 0)
      message("- Number of drift conditions = %d\n", nDrift);
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

double ALikelihood::computeLogLikelihood(bool verbose)
{
  _updateModel(verbose);

  if (_model->getNDriftEquation() > 0)
  {
    // Calculate t(L-1) %*% D-1 %*% L-1 applied to X (L and D from Vecchia)
    _computeCm1X();

    // Calculate XtCm1X = Xt * Cm1 * X
    _XtCm1X.resize(_X.getNCols(), _X.getNCols());
    _XtCm1X.prodMatMatInPlace(&_X, &_Cm1X, true, false);

    // Construct ZtCm1X = Zt * Cm1 * X and perform its Cholesky decomposition
    // workaround to create a shared_ptr which is not deleted at the end of the scope
    VectorDouble ZtCm1X = _Cm1X.prodVecMat(_Y);
    CholeskyDense XtCm1XChol(_XtCm1X);
    if (!XtCm1XChol.isReady())
    {
      messerr("Cholesky decomposition of XtCm1X matrix failed");
      return TEST;
    }

    // Calculate beta = (XtCm1X)-1 * ZtCm1X
    if (XtCm1XChol.solve(ZtCm1X, _beta))
    {
      messerr("Error when calculating Likelihood");
      return TEST;
    }
    // model->setBetaHat(beta);

    if (verbose)
      VH::dump("Optimal Drift coefficients = ", _beta);

    // Center the data by the optimal drift: Y = Y - beta * X
    VH::subtractInPlace(_Y, _X.prodMatVec(_beta));
  }

  // Calculate t(L-1) %*% D-1 %*% L-1 applied to Y (L and D from Vecchia)

  _computeCm1Y();

  // Calculate the log-determinant

  double logdet = _computeLogDet();
  // Calculate quad = Zt * Cm1Z
  double quad = VH::innerProduct(_Y, _Cm1Y);

  // Derive the log-likelihood
  Id size        = static_cast<Id>(_Y.size());
  double loglike = -0.5 * (logdet + quad + size * log(2. * GV_PI));
  if (_reml && _model->getNDriftEquation() > 0)
  {
    CholeskyDense XtCm1XChol(_XtCm1X);
    loglike -= 0.5 * XtCm1XChol.computeLogDeterminant();
  }

  // Optional printout
  if (verbose)
  {
    message("Log-Determinant = %lf\n", logdet);
    message("Quadratic term  = %lf\n", quad);
    message("Log-likelihood  = %lf\n", loglike);
  }
  return loglike;
}

double ALikelihood::computeCost(bool verbose)
{
  return -computeLogLikelihood(verbose);
}
} // namespace gstlrn