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
#include "Covariances/CorMatern.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/MathFunc.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/CorAniso.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Covariances/CovContext.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpaceComposite.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"

#include <cmath>
#include <memory>
#include <vector>

namespace gstlrn
{
CorMatern::CorMatern(const VectorDouble& ranges,
                     const VectorDouble& angles,
                     const VectorDouble& coeffScales,
                     const VectorDouble& params,
                     const MatrixSymmetric& sigma,
                     bool flagRange)
  : ACov()
  , _nVar(static_cast<Id>(params.size()))
  , _corRef(std::shared_ptr<const CorAniso>(CorAniso::createAnisotropic(CovContext(1, static_cast<Id>(ranges.size())), ECov::MATERN, ranges, params[0], angles, flagRange)))
  , _corMatern(*_corRef)
  , _coeffScales(coeffScales)
  , _params(params)
  , _C0(sigma)
  , _Rr(_nVar)
  , _Nu(_nVar)
  , _Tau(_nVar)
{
  if (static_cast<Id>(_coeffScales.size()) != _nVar)
  {
    messerr("CorMatern: inconsistent size between coeffScales and params");
    messerr("CorMatern: coeffScales size = %d, params size = %d", _coeffScales.size(), _nVar);
    _nVar        = 0;
    _coeffScales = VectorDouble();
    _params      = VectorDouble();
    _C0          = MatrixSymmetric(0);
    _Rr          = MatrixSymmetric(0);
    _Nu          = MatrixSymmetric(0);
    _Tau         = MatrixSymmetric(0);
    return;
  }
  //_coeffScales.push_front(1.);
  setContext(_corMatern.getContext());
  _ctxt.setNVar(_nVar);
  for (Id ivar = 0; ivar < _nVar; ivar++)
  {
    _Tau.setValue(ivar, ivar, 1.);
    for (Id jvar = ivar + 1; jvar < _nVar; jvar++)
    {
      double scalei  = _coeffScales[ivar];
      double scalej  = _coeffScales[jvar];
      double scaleij = computeScale(ivar, jvar);
      double nui     = _params[ivar];
      double nuj     = _params[jvar];
      double nuij    = computeParam(ivar, jvar);
      double gni     = exp(loggamma(nui));
      double gnj     = exp(loggamma(nuj));
      double gnij    = exp(loggamma(nuij));
      double ratioi  = gni / pow(scalei, 2. * nui);
      double ratioj  = gnj / pow(scalej, 2. * nuj);
      double ratioij = gnij / pow(scaleij, 2. * nuij);
      double val     = ratioij / sqrt(ratioi * ratioj);
      _Rr.setValue(ivar, jvar, scaleij);
      _Nu.setValue(ivar, jvar, nuij);
      _Tau.setValue(ivar, jvar, val);
    }
  }
  // checking the consistancy of the model
  MatrixSymmetric qcMat(_nVar);

  for (Id ivar = 0; ivar < _nVar; ivar++)
  {
    for (Id jvar = ivar; jvar < _nVar; jvar++)
    {
      double c0  = _C0.getValue(ivar, jvar);
      double tau = _Tau.getValue(ivar, jvar);
      qcMat.setValue(ivar, jvar, c0 / tau);
    }
  }
  if (!qcMat.isDefinitePositive())
  {
    messerr("CorMatern: inconsistent covariance matrix");
    _nVar        = 0;
    _coeffScales = VectorDouble();
    _params      = VectorDouble();
    _C0          = MatrixSymmetric(0);
    _Rr          = MatrixSymmetric(0);
    _Nu          = MatrixSymmetric(0);
    _Tau         = MatrixSymmetric(0);
    return;
  }
}

CorMatern::CorMatern(const CorMatern& r)
  : ACov(r)
  , _corMatern(r._corMatern)
{
  _nVar        = r._nVar;
  _corRef      = r._corRef;
  _coeffScales = r._coeffScales;
  _params      = r._params;
  _C0          = r._C0;
  _Rr          = r._Rr;
  _Nu          = r._Nu;
  _Tau         = r._Tau;
}

CorMatern& CorMatern::operator=(const CorMatern& r)
{
  if (this != &r)
  {
    ACov::operator=(r);
    _nVar        = r._nVar;
    _corMatern   = r._corMatern;
    _coeffScales = r._coeffScales;
    _params      = r._params;
    _C0          = r._C0;
    _Rr          = r._Rr;
    _Nu          = r._Nu;
    _Tau         = r._Tau;
  }
  return *this;
}

double CorMatern::computeScale(Id ivar, Id jvar) const
{
  if (ivar == jvar)
  {
    return _coeffScales[ivar];
  }
  double ci2 = _coeffScales[ivar] * _coeffScales[ivar];
  double cj2 = _coeffScales[jvar] * _coeffScales[jvar];

  return sqrt(0.5 * (ci2 + cj2));
}

double CorMatern::computeParam(Id ivar, Id jvar) const
{
  if (ivar == jvar)
  {
    return _params[ivar];
  }

  return 0.5 * (_params[ivar] + _params[jvar]);
}
double CorMatern::evalSpectrum(const VectorDouble& freq, Id ivar, Id jvar) const
{
  _corMatern.setParam(computeParam(ivar, jvar));

  VectorDouble scales = _corRef->getScales();
  VectorDouble angles = _corRef->getAnisoAngles();
  double correcScale  = computeScale(ivar, jvar);
  for (size_t idim = 0; idim < getSpace()->getNDim(); idim++)
  {
    scales[idim] /= correcScale;
  }
  _corMatern.setRotationAnglesAndRadius(angles, VectorDouble(), scales);
  return _C0.getValue(ivar, jvar) * _corMatern.evalSpectrum(freq, ivar, jvar);
}

double CorMatern::evalSpectrumRatio(const VectorDouble& freq, Id ivar, Id jvar, const ACov* cov0) const
{
  DECLARE_UNUSED(cov0)
  return evalSpectrum(freq, ivar, jvar) / _corRef->evalSpectrum(freq);
}

MatrixDense CorMatern::simulateSpectralOmega(Id nb) const
{
  return _corRef->simulateSpectralOmega(nb);
}

CorMatern::~CorMatern()
{
}

void CorMatern::_optimizationSetTarget(SpacePoint& pt) const
{
  DECLARE_UNUSED(pt)
}

void CorMatern::_optimizationPreProcess(Id mode, const std::vector<SpacePoint>& ps) const
{
  DECLARE_UNUSED(mode)
  DECLARE_UNUSED(ps)
  // _covS->_optimizationPreProcess(p);
  // _covTemp->_optimizationPreProcess(p);
}

void CorMatern::_optimizationPostProcess() const
{
  //_covS->optimizationPostProcess();
  //_covTemp->optimizationPostProcess();
}

double CorMatern::_eval(const SpacePoint& p1,
                        const SpacePoint& p2,
                        Id ivar,
                        Id jvar,
                        const CovCalcMode* mode) const
{
  _corMatern.setParam(computeParam(ivar, jvar));

  VectorDouble scales = _corRef->getScales();
  VectorDouble angles = _corRef->getAnisoAngles();
  double correcScale  = computeScale(ivar, jvar);
  for (size_t idim = 0; idim < getSpace()->getNDim(); idim++)
  {
    scales[idim] /= correcScale;
  }
  _corMatern.setRotationAnglesAndRadius(angles, VectorDouble(), scales);
  return _C0.getValue(ivar, jvar) * _corMatern.evalCov(p1, p2, 0, 0, mode);
}

} // namespace gstlrn