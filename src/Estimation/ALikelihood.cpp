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
#include "LinearOp/ASimulableMatrix.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Model/ModelGeneric.hpp"
#include "Stats/Classical.hpp"

namespace gstlrn
{
  ALikelihood::ALikelihood(ModelGeneric* model, const Db* db, bool reml)
    : AModelOptim(model)
    , ASimulableMatrix()
    , _db(db)
    , _reml(reml)
  {
    _nDrift = _model->getNDriftEquation();
  }

  double ALikelihood::computeLogDet(Id nMC) const
  {
    DECLARE_UNUSED(nMC);
    return -_computeLogDet();
  }

  Id ALikelihood::getSize() const
  {
    return _Y.size();
  }

  Id
    ALikelihood::_addSimulateToDest(const constvect whitenoise, vect outv) const
  {
    DECLARE_UNUSED(whitenoise);
    DECLARE_UNUSED(outv);
    messerr(
      "To perform non conditional simulation with this covariance, use simtub "
      "instead");
    return 0;
  }

  Id ALikelihood::_addToDest(constvect inv, vect outv) const
  {
    _solveQ(inv, outv);
    return 0;
  }

  void ALikelihood::initLikelihood(bool verbose)
  {
    MatrixSymmetric vars = dbVarianceMatrix(_db);
    double hmax = _db->getExtensionDiagonal();
    double vmax = _db->getColumnByLocator(ELoc::Z, 0).maximum();
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
        printVector(_model->getMeans(), "Constant Mean(s)", true, true);
    }

    // If Drift function is present, evaluate the optimal Drift coefficients
    if (_nDrift > 0)
    {
      // Extract the matrix of drifts at samples X
      _X = _model->evalDriftMat(_db);

      _beta.resize(_nDrift);
    }
  }

  void ALikelihood::_initLikelihoodForOptim(bool verbose)
  {
    initLikelihood(verbose);
    MatrixSymmetric vars = dbVarianceMatrix(_db);
    double hmax = _db->getExtensionDiagonal();
    double vmax = _db->getColumnByLocator(ELoc::Z, 0).maximum();
    setEnvironment(vars, hmax, EPSILON6, 0., vmax);
  }

  void ALikelihood::updateModel(bool verbose)
  {
    _updateModel(verbose);
  }

  bool ALikelihood::_calculateBeta(bool verbose)
  {
    // Calculate t(L-1) %*% D-1 %*% L-1 applied to X (L and D from Vecchia)
    _computeCm1X();

    // Calculate XtCm1X = Xt * Cm1 * X
    _XtCm1X.resize(_X.getNCols(), _X.getNCols());
    AMatrix::prodMatMatInPlace(_XtCm1X, _X, _Cm1X, true, false);

    // Construct ZtCm1X = Zt * Cm1 * X and perform its Cholesky decomposition
    // workaround to create a shared_ptr which is not deleted at the end of the scope
    VectorDouble ZtCm1X = AMatrix::product(_Y, _Cm1X);
    CholeskyDense XtCm1XChol(_XtCm1X);
    if (!XtCm1XChol.isReady())
    {
      messerr("Cholesky decomposition of XtCm1X matrix failed");
      return false;
    }

    // Calculate beta = (XtCm1X)-1 * ZtCm1X
    if (XtCm1XChol.solve(ZtCm1X, _beta))
    {
      messerr("Error when calculating Likelihood");
      return false;
    }
    _model->setBetaHat(_beta);

    if (verbose)
      printVector(_beta, "Optimal Drift coefficients = ", true, true);

    return true;
  }

  double ALikelihood::computeLogLikelihood(bool flagPrint, bool verbose)
  {
    updateModel(verbose);

    if (
      _model->getTransform()
      != nullptr) // TODO do it only in init if no parameters in transform (e.g logNormal)
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
      // Calculate the set of Drift Coefficients
      if (!_calculateBeta(verbose)) return TEST;

      // Center the data by the optimal drift: Yc = Y - beta * X
      VH::subtractInPlace(AMatrix::product(_X, _beta), _Y, _Yc);
    }

    // Calculate t(L-1) %*% D-1 %*% L-1 applied to Y (L and D from Vecchia)

    _computeCm1Yc();

    // Calculate the log-determinant

    double logdet = _computeLogDet();
    // Calculate quad = Zt * Cm1Z
    double quad = _Yc.innerProduct(_Cm1Yc);

    // Derive the log-likelihood
    Id size = static_cast<Id>(_Yc.size());
    double loglike = -0.5 * (logdet + quad + size * log(2. * GV_PI));
    if (_reml && _nDrift > 0)
    {
      CholeskyDense XtCm1XChol(_XtCm1X);
      loglike -= 0.5 * XtCm1XChol.computeLogDeterminant();
    }

    double logjac = _model->getTransform() == nullptr
                    ? 0.
                    : _model->getTransform()->evalLogJacobianVec(_Y);

    loglike -= logjac;
    // Optional printout
    if (flagPrint)
    {
      message("Log-Determinant = %lf\n", logdet);
      message("Quadratic term  = %lf\n", quad);
      if (_model->getTransform() != nullptr)
        message("Jacobian term  = %lf\n", logjac);
      message("Log-likelihood  = %lf\n", loglike);
    }
    return loglike;
  }

  double ALikelihood::computeCost(bool flagPrint, bool verbose)
  {
    return -computeLogLikelihood(flagPrint, verbose);
  }
} // namespace gstlrn
