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
#include "Covariances/KernelMatern.hpp"
#include "Basic/Law.hpp"
#include "Basic/MathFunc.hpp"
#include "Basic/Utilities.hpp"
#include "Covariances/CovContext.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Simulation/TurningBandOperate.hpp"
#include "geoslib_define.h"

#include <cmath>

#define MAXTAB 100
namespace gstlrn
{
static bool bessel_Old_Style = false;

KernelMatern::KernelMatern(const CovContext& ctxt)
  : AKernel(ECov::MATERN, ctxt)
  , _correc(1.)
  , _markovCoeffs(VectorDouble())
{
  setParam(1);
  computeMarkovCoeffs(2);
  // TODO compute blin (rapatrier de PrecisionOp.cpp
}

KernelMatern::KernelMatern(const KernelMatern& r)
  : AKernel(r)
  , _correc(r._correc)
  , _markovCoeffs(r._markovCoeffs)
{
}

KernelMatern& KernelMatern::operator=(const KernelMatern& r)
{
  if (this != &r)
  {
    AKernel::operator=(r);
  }
  return *this;
}

KernelMatern::~KernelMatern()
{
}

double KernelMatern::getScadef() const
{
  return sqrt(12. * getParam());
}

double KernelMatern::_evaluateCov(double h) const
{
  if (bessel_Old_Style)
  {
    return _oldMatern(h);
  }
  return _newMatern(h);
}

double KernelMatern::_newMatern(double h) const
{
  if (h == 0) return 1;
  double third = getParam();
  return 2. * pow(h / 2., third) * _besselK(getParam(), h) / exp(loggamma(third));
}

double KernelMatern::_evaluateCovFirstDerivative(double h) const
{
  if (h == 0) return 1;
  double nu    = getParam();
  double ratio = pow(2., 1. - nu) / exp(loggamma(nu));
  double term1 = nu * pow(h, nu - 1) * _besselK(nu, h);
  double term2 = 0.5 * pow(h, nu) * (_besselK(nu - 1, h) + _besselK(nu + 1, h));
  return ratio * (term1 - term2);
}

double KernelMatern::_besselK(double nu, double h)
{
#if !defined(__cpp_lib_math_special_functions)
  double TAB[MAXTAB];
  Id nb = static_cast<Id>(floor(nu));
  if (nu <= 0 || nb >= MAXTAB) return (0.);
  double alpha = nu - nb;
  if (besselk(h, alpha, nb + 1, TAB) < nb + 1) return 0.;
  return TAB[nb];
#else
  return std::cyl_bessel_k(nu, h);
#endif
}

double KernelMatern::_oldMatern(double h) const
{
  double TAB[MAXTAB];
  double cov   = 0.;
  double third = getParam();
  Id nb        = static_cast<Id>(floor(third));
  double alpha = third - nb;
  if (third <= 0 || nb >= MAXTAB) return (0.);
  double coeff = (h > 0) ? pow(h / 2., third) : 1.;
  cov          = 1.;
  if (h > 0)
  {
    if (besselk(h, alpha, nb + 1, TAB) < nb + 1) return 0.;
    cov = 2. * coeff * TAB[nb] / exp(loggamma(third));
  }
  return (cov);
}

String KernelMatern::getFormula() const
{
  return "C(h)=\\frac{2^{1-\\nu}}{\\Gamma(\\nu)} h^\\nu K_{\\nu}( h )";
}

double KernelMatern::evaluateSpectrum(double freq) const
{

 

  size_t ndim = getContext().getNDim();
  double param = getParam();
  size_t ndims2  = 0.5 * ndim;
  double alpha = param + ndims2;
  double val = pow(2, ndim) / getCorrec() / pow (1 + (freq * freq), alpha);
  return val;

  /*
  Id ndim     = getContext().getNDim();
  double alpha = (double)ndim / 2. + getParam();
  return 1. / pow(1. + freq, alpha);
  */
}

void KernelMatern::computeMarkovCoeffs(Id ndim)
{
  double param  = getParam();
  double ndims2 = (static_cast<double>(ndim)) / 2.;
  double alpha  = param + ndims2;
  auto p        = getClosestInteger(alpha);
  Id ndimp      = p + 1;
  _markovCoeffs.resize(ndimp);
  for (Id i = 0; i < ndimp; i++)
  {
    _markovCoeffs[i] = ut_cnp(p, i);
  }
  computeCorrec(ndim);
}

void KernelMatern::computeCorrec(Id ndim)
{
  double g0, ndims2, gammap, gammaa;
  ndims2  = (static_cast<double>(ndim)) / 2.;
  gammap  = exp(loggamma(getParam()));
  gammaa  = exp(loggamma(getParam() + ndims2));
  g0      = pow(4. * GV_PI, ndims2);
  _correc = gammap / (g0 * gammaa);
}

VectorDouble KernelMatern::getMarkovCoeffs() const
{
  return _markovCoeffs;
}

double KernelMatern::simulateTurningBand(double t0, TurningBandOperate& operTB) const
{
  if (getParam() > 0.5)
    return operTB.cosineOne(t0);
  return operTB.spectralOne(t0);
}

MatrixDense KernelMatern::simulateSpectralOmega(Id nb) const
{
  auto ndim    = static_cast<Id>(getContext().getNDim());
  double param = getParam();
  MatrixDense mat(nb, ndim);

  for (Id irow = 0; irow < nb; irow++)
  {
    double scale = sqrt(1 / law_gamma(param) / 2);
    for (Id icol = 0; icol < ndim; icol++)
      mat.setValue(irow, icol, scale * law_gaussian());
  }
  return mat;
}

VectorDouble KernelMatern::_evaluateSpectrumOnSphere(Id n, double scale) const
{
  double scale2 = scale * scale;
  double mu     = getParam();
  double alpha  = mu + 1.;

  VectorDouble sp(1 + n, 0.);

  for (Id k = 0; k <= n; k++)
    sp[k] = (2. * k + 1.) / (4. * GV_PI) / pow(1. + (scale2 * k * (k + 1.)), alpha);

  sp.normalizeInPlace(1);
  return sp;
}

void bessel_set_old_style(bool style)
{
  bessel_Old_Style = style;
}
} // namespace gstlrn
