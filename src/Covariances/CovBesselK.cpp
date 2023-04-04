/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "geoslib_old_f.h"

#include "Covariances/CovBesselK.hpp"
#include "Covariances/CovContext.hpp"
#include "Basic/MathFunc.hpp"

#include "math.h"

#define MAXTAB 100

CovBesselK::CovBesselK(const CovContext& ctxt)
: ACovFunc(ECov::BESSEL_K, ctxt)
  ,_correc(1.)
  ,_markovCoeffs(VectorDouble())
{
  setParam(1);
  computeMarkovCoeffs(2);
  //TODO compute blin (rapatrier de PrecisionOp.cpp
}

CovBesselK::CovBesselK(const CovBesselK &r)
: ACovFunc(r)
  ,_correc(r._correc)
  ,_markovCoeffs(r._markovCoeffs)
{
}

CovBesselK& CovBesselK::operator=(const CovBesselK &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovBesselK::~CovBesselK()
{
}

double CovBesselK::getScadef() const
{
  return sqrt(12. * getParam());
}

double CovBesselK::_evaluateCov(double h) const
{
  static double TAB[MAXTAB];

  double cov = 0.;
  double third = getParam();
  int nb = (int) floor(third);
  double alpha = third - nb;
  if (third <= 0 || nb >= MAXTAB) return (0.);
  double coeff = (h > 0) ? pow(h / 2., third) : 1.;
  cov = 1.;
  if (h > 0)
  {
    if (bessel_k(h, alpha, nb + 1, TAB) < nb + 1) return (cov);
    cov = 2. * coeff * TAB[nb] / exp(loggamma(third));
  }
  return (cov);
}

String CovBesselK::getFormula() const
{
  return "C(h)=\\frac{ \\left( \\frac{h}{a_t} \\right)^\\alpha}{2^{\\alpha-1}\\Gamma(\\alpha)}K_{-\\alpha} \\left( \\frac{h}{a_t} \\right)";
}

double CovBesselK::_evaluateCovOnSphere(double scale, int degree) const
{
  double kappa2 = 1. / ( scale * scale );
  double cons = 1. / (4 * GV_PI);
  return  cons * (2. * degree + 1.) / pow(kappa2 + degree * (degree + 1), 1. + getParam());
}

double CovBesselK::evaluateSpectrum(double freq, int ndim) const
{
  double alpha = (double) ndim / 2. + getParam();
  return 1. /  pow(1. + freq, alpha);
}

void CovBesselK::computeMarkovCoeffs(int ndim)
{
  double param = getParam();
  double ndims2 = ((double) ndim) / 2.;
  double alpha = param + ndims2;
  int p = getClosestInteger(alpha);
  int ndimp = p + 1;
  _markovCoeffs.resize(ndimp);
  for (int i = 0; i < ndimp; i++)
  {
    _markovCoeffs[i] = (double)ut_cnp(p, i);
  }
  computeCorrec(ndim);

}

void CovBesselK::computeCorrec(int ndim)
{
  double g0, ndims2, gammap, gammaa;
  ndims2 = ((double) ndim) / 2.;
  gammap = exp(loggamma(getParam()));
  gammaa = exp(loggamma(getParam() + ndims2));
  g0 = pow(4. * GV_PI, ndims2);
  _correc = gammap / (g0 * gammaa);
}

VectorDouble CovBesselK::getMarkovCoeffs()const
{
  return _markovCoeffs;
}
