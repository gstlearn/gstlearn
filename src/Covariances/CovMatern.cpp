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
#include "Covariances/CovMatern.hpp"
#include "Covariances/CovContext.hpp"
#include "Simulation/TurningBandOperate.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Basic/Law.hpp"
#include "Basic/MathFunc.hpp"
#include "Basic/Utilities.hpp"

#include "math.h"

#define MAXTAB 100

CovMatern::CovMatern(const CovContext& ctxt)
: ACovFunc(ECov::MATERN, ctxt)
  ,_correc(1.)
  ,_markovCoeffs(VectorDouble())
{
  setParam(1);
  computeMarkovCoeffs(2);
  //TODO compute blin (rapatrier de PrecisionOp.cpp
}

CovMatern::CovMatern(const CovMatern &r)
: ACovFunc(r)
  ,_correc(r._correc)
  ,_markovCoeffs(r._markovCoeffs)
{
}

CovMatern& CovMatern::operator=(const CovMatern &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovMatern::~CovMatern()
{
}

double CovMatern::getScadef() const
{
  return sqrt(12. * getParam());
}

double CovMatern::_evaluateCov(double h) const
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
    if (matern(h, alpha, nb + 1, TAB) < nb + 1) return (cov);
    cov = 2. * coeff * TAB[nb] / exp(loggamma(third));
  }
  return (cov);
}

String CovMatern::getFormula() const
{
  return "C(h)=\\frac{ \\left( \\frac{h}{a_t} \\right)^\\alpha}{2^{\\alpha-1}\\Gamma(\\alpha)}K_{-\\alpha} \\left( \\frac{h}{a_t} \\right)";
}

double CovMatern::evaluateSpectrum(double freq) const
{
  int ndim = getContext().getNDim();
  double alpha = (double) ndim / 2. + getParam();
  return 1. /  pow(1. + freq, alpha);
}

void CovMatern::computeMarkovCoeffs(int ndim)
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

void CovMatern::computeCorrec(int ndim)
{
  double g0, ndims2, gammap, gammaa;
  ndims2 = ((double) ndim) / 2.;
  gammap = exp(loggamma(getParam()));
  gammaa = exp(loggamma(getParam() + ndims2));
  g0 = pow(4. * GV_PI, ndims2);
  _correc = gammap / (g0 * gammaa);
}

VectorDouble CovMatern::getMarkovCoeffs()const
{
  return _markovCoeffs;
}

double CovMatern::simulateTurningBand(double t0, TurningBandOperate &operTB) const
{
  if (getParam() > 0.5)
    return operTB.cosineOne(t0);
  return operTB.spectralOne(t0);
}

MatrixRectangular CovMatern::simulateSpectralOmega(int nb) const
{
  int ndim = getContext().getNDim();
  double param = getParam();
  MatrixRectangular mat(nb, ndim);

  for (int irow = 0; irow < nb; irow++)
  {
    double scale = sqrt(param / law_gamma(param));
    for (int icol = 0; icol < ndim; icol++)
      mat.setValue(irow, icol, scale * law_gaussian());
  }
  return mat;
}

VectorDouble CovMatern::_evaluateSpectrumOnSphere(int n, double scale) const
{
  double scale2 = scale * scale;
  double mu = getParam();
  double alpha = mu + 1.;

  VectorDouble sp(1+n, 0.);

  for (int k = 0; k <= n; k++)
    sp[k] = (2. * k + 1.) / (4. * GV_PI) / pow(1. + scale2 * k * (k + 1.), alpha);

  VH::normalize(sp,1);
  return sp;
}