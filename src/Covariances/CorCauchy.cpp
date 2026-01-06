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
#include "Covariances/CorCauchy.hpp"
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
CorCauchy::CorCauchy(const VectorDouble& ranges,
                     const VectorDouble& angles,
                     const VectorDouble& coeffScales,
                     const VectorDouble& params,
                     bool flagRange)
  : ACov()
  , _nVar(static_cast<Id>(params.size()))
  , _corRef(std::shared_ptr<const CorAniso>(CorAniso::createAnisotropic(CovContext(1, static_cast<Id>(ranges.size())), ECov::CAUCHY, ranges, params[0], angles, flagRange)))
  , _corCauchy(*_corRef)
  , _coeffScales(coeffScales)
  , _params(params)
  , _C0(_nVar)
  , _Rr(_nVar)
  , _Nu(_nVar)
{
  if (static_cast<Id>(_coeffScales.size()) != _nVar)
  {
    messerr("CorCauchy: inconsistent size between coeffScales and params");
    messerr("CorCauchy: coeffScales size = %d, params size = %d", _coeffScales.size(), _nVar);
    _nVar        = 0;
    _coeffScales = VectorDouble();
    _params      = VectorDouble();
    _C0          = MatrixSymmetric(0);
    _Rr          = MatrixSymmetric(0);
    _Nu          = MatrixSymmetric(0);
    return;
  }
  //_coeffScales.push_front(1.);
  setContext(_corCauchy.getContext());
  _ctxt.setNVar(_nVar);
  for (Id ivar = 0; ivar < _nVar; ivar++)
  {
    for (Id jvar = ivar; jvar < _nVar; jvar++)
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
      _C0.setValue(ivar, jvar, val);
    }
  }
}

CorCauchy::CorCauchy(const CorCauchy& r)
  : ACov(r)
  , _corCauchy(r._corCauchy)
{
  _nVar        = r._nVar;
  _corRef      = r._corRef;
  _coeffScales = r._coeffScales;
  _params      = r._params;
  _C0          = r._C0;
  _Rr          = r._Rr;
  _Nu          = r._Nu;
}

CorCauchy& CorCauchy::operator=(const CorCauchy& r)
{
  if (this != &r)
  {
    ACov::operator=(r);
    _nVar        = r._nVar;
    _corCauchy   = r._corCauchy;
    _coeffScales = r._coeffScales;
    _params      = r._params;
    _C0          = r._C0;
    _Rr          = r._Rr;
    _Nu          = r._Nu;
  }
  return *this;
}

double CorCauchy::computeScale(Id ivar, Id jvar) const
{
  if (ivar == jvar)
  {
    return _coeffScales[ivar];
  }
  double ci2 = _coeffScales[ivar] * _coeffScales[ivar];
  double cj2 = _coeffScales[jvar] * _coeffScales[jvar];

  return sqrt(0.5 * (ci2 + cj2));
}

double CorCauchy::computeParam(Id ivar, Id jvar) const
{
  if (ivar == jvar)
  {
    return _params[ivar];
  }

  return 0.5 * (_params[ivar] + _params[jvar]);
}
double CorCauchy::evalSpectrum(const VectorDouble& freq, Id ivar, Id jvar) const
{
  _corCauchy.setParam(computeParam(ivar, jvar));

  VectorDouble scales = _corRef->getScales();
  VectorDouble angles = _corRef->getAnisoAngles();
  double correcScale  = computeScale(ivar, jvar);
  for (size_t idim = 0; idim < getSpace()->getNDim(); idim++)
  {
    scales[idim] /= correcScale;
  }
  _corCauchy.setRotationAnglesAndRadius(angles, VectorDouble(), scales);
  return _C0.getValue(ivar, jvar) * _corCauchy.evalSpectrum(freq, ivar, jvar);
}

double CorCauchy::evalSpectrumRatio(const VectorDouble& freq, Id ivar, Id jvar, const ACov* cov0) const
{
  DECLARE_UNUSED(cov0)
  return evalSpectrum(freq, ivar, jvar) / _corRef->evalSpectrum(freq);
}

MatrixDense CorCauchy::simulateSpectralOmega(Id nb) const
{
  return _corRef->simulateSpectralOmega(nb);
}

CorCauchy::~CorCauchy()
{
}

CorCauchy* CorCauchy::create(const ECov& type,
                                 const CovContext& ctxt,
                                 const VectorDouble& params, // = (nu_i)
                                 const VectorDouble& coeffScales, // = (r_i)
                                 const VectorDouble& ranges,
                                 const VectorDouble& angles,
                                 bool flagRange)
{
  if (type != ECov::CAUCHY)
  {
    messerr("This function implements the Cauchy model");
    return nullptr;
  }
  Id ndim = ctxt.getNDim(); // dimension of the space
  Id nvar = ctxt.getNVar(); // number of variables
  if (params.length() != nvar)
  {
    messerr("Cauchy model requires %d parameters (%d provided)", nvar, params.length());
    return nullptr;
  }
  if (coeffScales.length() != nvar)
  {
    messerr("Cauchy model requires %d scaling factors (%d provided)", nvar, coeffScales.length());
    return nullptr;
  }
  if (ranges.length() != ndim)
  {
    messerr("Inconsistent number of ranges (%d vs. ndim = %d)", ranges.length(), ndim);
    return nullptr;
  }
  Id nang = angles.length();
  if (nang > 0)
  {
    if (nang != ndim)
    {
      messerr("Inconsistent number of angles (%d)", nang);
      return nullptr;
    }
  }
  return new CorCauchy(ranges, angles, coeffScales, params, flagRange);
}

void CorCauchy::_optimizationSetTarget(SpacePoint& pt) const
{
  DECLARE_UNUSED(pt)
}

void CorCauchy::_optimizationPreProcess(Id mode, const std::vector<SpacePoint>& ps) const
{
  DECLARE_UNUSED(mode)
  DECLARE_UNUSED(ps)
  // _covS->_optimizationPreProcess(p);
  // _covTemp->_optimizationPreProcess(p);
}

void CorCauchy::_optimizationPostProcess() const
{
  //_covS->optimizationPostProcess();
  //_covTemp->optimizationPostProcess();
}

double CorCauchy::_eval(const SpacePoint& p1,
                        const SpacePoint& p2,
                        Id ivar,
                        Id jvar,
                        const CovCalcMode* mode) const
{
  _corCauchy.setParam(computeParam(ivar, jvar));

  VectorDouble scales = _corRef->getScales();
  VectorDouble angles = _corRef->getAnisoAngles();
  double correcScale  = computeScale(ivar, jvar);
  for (size_t idim = 0; idim < getSpace()->getNDim(); idim++)
  {
    scales[idim] /= correcScale;
  }
  _corCauchy.setRotationAnglesAndRadius(angles, VectorDouble(), scales);
  return _C0.getValue(ivar, jvar) * _corCauchy.evalCov(p1, p2, 0, 0, mode);
}

} // namespace gstlrn