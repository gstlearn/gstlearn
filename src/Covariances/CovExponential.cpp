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
#include "Covariances/CovExponential.hpp"

#include "Covariances/CovContext.hpp"
#include "Simulation/TurningBandOperate.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Basic/Law.hpp"
#include "Basic/VectorHelper.hpp"

#include "math.h"

CovExponential::CovExponential(const CovContext& ctxt)
: ACovFunc(ECov::EXPONENTIAL, ctxt)
{
}

CovExponential::CovExponential(const CovExponential &r)
: ACovFunc(r)
{
}

CovExponential& CovExponential::operator=(const CovExponential &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovExponential::~CovExponential()
{
}

double CovExponential::getScadef() const
{
  return (2.995732);
}

double CovExponential::_evaluateCov(double h) const
{
  if (h > 100) return (0.);
  double cov = exp(-h);
  return (cov);
}

String CovExponential::getFormula() const
{
  return "C(h)=exp \\left( -\\frac{h}{a_t} \\right)";
}

double CovExponential::simulateTurningBand(double t0, TurningBandOperate &operTB) const
{
  return operTB.spectralOne(t0);
}

MatrixRectangular CovExponential::simulateSpectralOmega(int nb) const
{
  int ndim = getContext().getNDim();
  double param = 0.5;
  MatrixRectangular mat(nb, ndim);

  for (int irow = 0; irow < nb; irow++)
  {
    double scale = sqrt(param / law_gamma(param));
    for (int icol = 0; icol < ndim; icol++)
      mat.setValue(irow, icol, scale * law_gaussian());
  }
  return mat;
}

double CovExponential::_evaluateCovOnSphere(double alpha,
                                            double scale,
                                            int degree) const
{
  double nu = scale * getScadef();
  return exp(-nu * alpha);
}

VectorDouble CovExponential::_evaluateSpectrumOnSphere(int n, double scale) const
{
  double nu = scale * getScadef();
  double nu2 = nu * nu;
  double expnu = exp(-nu * GV_PI);

  VectorDouble sp(1+n, 0.);
  int k;

  k = 0;
  sp[k] = 1. / 2. * (1. + expnu) / (1. + nu2);
  k = 1;
  sp[k] = 3. / 2. * (1. - expnu) / (4. + nu2);

  while(1)
  {
    k++;
    if (k >= n + 1) break;
    sp[k] = (2. * k + 1.) / (2. * k - 3.) * (nu2 + (k - 2.) * (k - 2.))
        / (nu2 + (k + 1.) * (k + 1.)) * sp[k - 2];
  }

  VH::normalize(sp, 1);

  return sp;
}
