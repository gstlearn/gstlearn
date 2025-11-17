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
#include "Stats/Classical.hpp"

namespace gstlrn
{
ALikelihood::ALikelihood(ModelGeneric* model,
                         const Db* db,
                         bool reml)
  : AModelOptim(model)
  , _db(db)
  , _reml(reml)
{
  _nDrift = _model->getNDriftEquation();
}

ALikelihood::ALikelihood(const ALikelihood& r)
  : AModelOptim(r)
  , _db(r._db)
  , _Z(r._Z)
  , _Y(r._Y)
  , _Yc(r._Yc)
  , _X(r._X)
  , _beta(r._beta)
  , _Cm1X(r._Cm1X)
  , _Cm1Yc(r._Cm1Yc)
  , _XtCm1X(r._XtCm1X)
  , _reml(r._reml)
  , _nDrift(r._nDrift) {};

ALikelihood& ALikelihood::operator=(const ALikelihood& r)
{
  if (this != &r)
  {
    AModelOptim::operator=(r);
    _db     = r._db;
    _Z      = r._Z;
    _Y      = r._Y;
    _Yc     = r._Yc;
    _X      = r._X;
    _beta   = r._beta;
    _Cm1X   = r._Cm1X;
    _Cm1Yc  = r._Cm1Yc;
    _XtCm1X = r._XtCm1X;
    _reml   = r._reml;
    _nDrift = r._nDrift;
  }
  return *this;
}

ALikelihood::~ALikelihood()
{
}

void ALikelihood::_initLikelihood(bool verbose)
{
  MatrixSymmetric vars = dbVarianceMatrix(_db);
  double hmax          = _db->getExtensionDiagonal();
  double vmax = VH::maximum(_db->getColumnByLocator(ELoc::Z, 0));
  setEnvironment(vars, hmax, EPSILON6, 0., vmax);
  Id nvar = _db->getNLoc(ELoc::Z);
  if (nvar < 1)
  {
    messerr("The 'db' should have at least one variable defined");
  }

  // Establish the vector of multivariate data
  if (_model->getTransform() == nullptr)
  {
    if (_nDrift > 0)
    {
      _Y = _db->getColumnsActiveAndDefined(ELoc::Z);
      _Yc.resize(_Y.size());
    }
    else
      _Yc = _db->getColumnsActiveAndDefined(ELoc::Z, _model->getMeans());
  }
  else
  {
    _Z = _db->getColumnsActiveAndDefined(ELoc::Z);
    _Y.resize(_Z.size());
    _Yc.resize(_Z.size());
  }

  Id size = static_cast<Id>(_Yc.size());
  if (verbose)
  {
    message("Likelihood calculation:\n");
    message("- Number of active samples     = %d\n", _db->getNSample(true));
    message("- Number of variables          = %d\n", nvar);
    message("- Length of Information Vector = %d\n", size);
    if (_nDrift > 0)
      message("- Number of drift conditions   = %d\n", _nDrift);
    else
      VH::dump("Constant Mean(s)", _model->getMeans());
  }

  // If Drift function is present, evaluate the optimal Drift coefficients
  if (_nDrift > 0)
  {
    // Extract the matrix of drifts at samples X
    _X = _model->evalDriftMat(_db);

    _beta.resize(_nDrift);
  }
}

double ALikelihood::computeLogLikelihood(bool flagPrint, bool verbose)
{
  _updateModel(verbose);

  if (_model->getTransform() != nullptr)//TODO do it only in init if no parameters in transform (e.g logNormal)
  {
    // Apply the transformation to data
    _model->getTransform()->inverseTransformVec(_Z, _Y);
    if (_nDrift == 0)
    {
      // Center the data by the means TODO do it!
      //_Yc = _Y - _model->getMeans();
    }
  }
  if (_nDrift > 0)
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
    _model->setBetaHat(_beta);

    if (verbose)
      VH::dump("Optimal Drift coefficients = ", _beta);

    // Center the data by the optimal drift: Yc = Y - beta * X
    VH::subtractInPlace(_X.prodMatVec(_beta), _Y, _Yc);
  }

  // Calculate t(L-1) %*% D-1 %*% L-1 applied to Y (L and D from Vecchia)

  _computeCm1Yc();

  // Calculate the log-determinant

  double logdet = _computeLogDet();
  // Calculate quad = Zt * Cm1Z
  double quad = _Yc.innerProduct(_Cm1Yc);

  // Derive the log-likelihood
  Id size        = static_cast<Id>(_Yc.size());
  double loglike = -0.5 * (logdet + quad + size * log(2. * GV_PI));
  if (_reml && _nDrift > 0)
  {
    CholeskyDense XtCm1XChol(_XtCm1X);
    loglike -= 0.5 * XtCm1XChol.computeLogDeterminant();
  }
  if (_model->getTransform() != nullptr)
  {
    // Add the Jacobian term
    double logjac = _model->getTransform()->evalLogJacobianVec(_Y);
    loglike -= logjac;
  }
  // Optional printout
  if (flagPrint)
  {
    message("Log-Determinant = %lf\n", logdet);
    message("Quadratic term  = %lf\n", quad);
    message("Log-likelihood  = %lf\n", loglike);
  }
  return loglike;
}

double ALikelihood::computeCost(bool flagPrint, bool verbose)
{
  return -computeLogLikelihood(flagPrint, verbose);
}
} // namespace gstlrn