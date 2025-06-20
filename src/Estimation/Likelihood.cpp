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
#include "Estimation/Likelihood.hpp"
#include "Estimation/ALikelihood.hpp"
#include "Tree/Ball.hpp"
#include "Db/Db.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Model/ModelGeneric.hpp"
#include "Stats/Classical.hpp"
#include "geoslib_define.h"

Likelihood::Likelihood(ModelGeneric* model,
                       const Db* db)
  : ALikelihood(model, db)
{
}

Likelihood::Likelihood(const Likelihood& r)
  : ALikelihood(r)
{
}

Likelihood& Likelihood::operator=(const Likelihood& r)
{
  if (this != &r)
  {
    ALikelihood::operator=(r);
  }
  return *this;
}

Likelihood::~Likelihood()
{
}

double logLikelihood(const Db* db,
                     ModelGeneric* model,
                     bool verbose)
{
  Likelihood* vec = Likelihood::createForOptim(model, db);
  double result   = vec->computeLogLikelihood(verbose);
  delete vec;
  return result;
}

Likelihood* Likelihood::createForOptim(ModelGeneric* model,
                                       const Db* db)
{
  auto* vec = new Likelihood(model, db);
  MatrixSymmetric vars = dbVarianceMatrix(db);
  double hmax          = db->getExtensionDiagonal();
  vec->setEnvironment(vars, hmax);
  vec->init();
  return vec;
}

void Likelihood::_computeCm1X()
{
  if (_covChol.solveMatrix(_X, _Cm1X))
  {
    messerr("Problem when solving a Linear System after Cholesky decomposition");
  }
}

void Likelihood::_computeCm1Z()
{
  _Cm1Y.resize(_Y.size());
  if (_covChol.solve(_Y, _Cm1Y))
  {
    messerr("Error when calculating Cm1Z");
  }
}

double Likelihood::_computeLogDet() const
{
  return _covChol.computeLogDeterminant();
}

void Likelihood::_updateModel(bool verbose)
{
  DECLARE_UNUSED(verbose);
  MatrixSymmetric cov = _model->evalCovMatSymInPlace(_cov, _db);
  _covChol.setMatrix(&_cov);
  if (!_covChol.isReady())
  {
    messerr("Cholesky decomposition of Covariance matrix failed");
  }
}