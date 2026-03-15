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
#include "Db/Db.hpp"
#include "Db/RankHandler.hpp"
#include "Estimation/ALikelihood.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Model/ModelGeneric.hpp"
#include "Space/SpacePoint.hpp"
#include "Stats/Classical.hpp"
#include "Tree/Ball.hpp"
#include "geoslib_define.h"

namespace gstlrn
{
Likelihood::Likelihood(ModelGeneric* model,
                       const Db* db,
                       bool reml)
  : ALikelihood(model, db, reml)
  , _cov(std::make_shared<MatrixSymmetric>(0))
{
  setAuthorizedAnalyticalGradients(_model->getTransform() == nullptr);
  setAuthorizedAnalyticalGradients(_model->getTransform() == nullptr);
}

Likelihood::Likelihood(const Likelihood& r)
  : ALikelihood(r)
  , _cov(r._cov)
{
}

Likelihood& Likelihood::operator=(const Likelihood& r)
{
  if (this != &r)
  {
    ALikelihood::operator=(r);
    _cov = r._cov;
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
  Likelihood vec(model, db);
  double result = vec.computeLogLikelihood(verbose);
  return result;
}

Likelihood* Likelihood::createForOptim(ModelGeneric* model,
                                       const Db* db,
                                       bool reml,
                                       bool verbose)
{
  auto* vec = new Likelihood(model, db, reml);
  vec->_initLikelihoodForOptim(verbose);
  return vec;
}

void Likelihood::_computeCm1X()
{
  if (_covChol.solveMatrix(_X, _Cm1X))
  {
    messerr("Problem when solving a Linear System after Cholesky decomposition");
  }
}

void Likelihood::_computeCm1Yc()
{
  _Cm1Yc.resize(_Yc.size());
  if (_covChol.solve(_Yc, _Cm1Yc))
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
  _model->manage(_db, nullptr);
  _model->evalCovMatSymInPlace(*_cov, _db);
  _covChol.setMatrix(*_cov);
}

void Likelihood::evalGrad(vect res)
{
  _temp.resize(_Yc.size());
  _gradCovMatTimesInvCov.resize(static_cast<Id>(_Yc.size()), static_cast<Id>(_Yc.size()));
  const auto invcov = _covChol.inverse();
  RankHandler rkh(_db);
  rkh.defineSampleRanks();
  const auto& gradcov = _model->getCovGradients();
  _gradCovMat.resize(static_cast<Id>(_Yc.size()), static_cast<Id>(_Yc.size()));
  CholeskyDense XtCm1XChol;
  MatrixSymmetric invXtCm1X;
  if (_reml && _model->getNDriftEquation() > 0)
  {
    XtCm1XChol.setMatrix(_XtCm1X);
    invXtCm1X = XtCm1XChol.inverse();
  }
  for (size_t iparam = 0; iparam < gradcov.size(); iparam++)
  {
    const auto& func = gradcov[iparam];
    _fillGradCovMat(rkh, func);
    AMatrix::prodVec(_temp, _gradCovMat, _Cm1Yc);
    // _gradCovMat.prodMatVecInPlace(_Cm1Yc, _temp);
    double dquad = -_Cm1Yc.innerProduct(_temp);
    res[iparam]  = 0.0;
    if (_reml && _model->getNDriftEquation() > 0)
    {
      MatrixSymmetric temp(_X.getNCols());
      temp.prodNormMatMatInPlace(&_Cm1X, &_gradCovMat, true);
      double dlogdetreml = MatrixDense::traceProd(invXtCm1X, temp);
      res[iparam] += 0.5 * dlogdetreml;
    }
    double dlogdet = MatrixDense::traceProd(invcov, _gradCovMat); // Warning: _gradCovMat is modified so the line
    // has to be after _gradCovMat.prodMatVecInPlace(_Cm1Y, _temp);
    res[iparam] += 0.5 * (dlogdet + dquad);
  }
}

void Likelihood::_fillGradCovMat(RankHandler& rkh, const covmaptype& gradcov)
{
  Id icur, jcur = 0;

  SpacePoint p1(_model->getSpace());
  SpacePoint p2(_model->getSpace());
  rkh.defineSampleRanks();

  for (Id jvar = 0; jvar < _model->getNVar(); jvar++)
  {
    const auto& indsj = rkh.getSampleRanksByVariable(jvar);

    for (const auto j: indsj)
    {
      icur = 0;
      _db->getSampleAsSPInPlace(p1, j);

      for (Id ivar = 0; ivar < _model->getNVar(); ivar++)
      {
        const auto& indsi = rkh.getSampleRanksByVariable(ivar);
        for (const auto i: indsi)
        {
          _db->getSampleAsSPInPlace(p2, i);

          _gradCovMat.setValue(icur, jcur, gradcov(p1, p2, ivar, jvar, nullptr));
          icur++;
        }
      }
      jcur++;
    }
  }
}
} // namespace gstlrn
