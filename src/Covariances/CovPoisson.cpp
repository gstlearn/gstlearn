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
#include "Covariances/CovPoisson.hpp"

#include "math.h"
#include "Basic/VectorHelper.hpp"
#include "Basic/MathFunc.hpp"
#include "Basic/Law.hpp"
#include "Covariances/CovContext.hpp"

CovPoisson::CovPoisson(const CovContext& ctxt)
: ACovFunc(ECov::POISSON, ctxt)
{
  setParam(1);
}

CovPoisson::CovPoisson(const CovPoisson &r)
: ACovFunc(r)
{
}

CovPoisson& CovPoisson::operator=(const CovPoisson &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovPoisson::~CovPoisson()
{
}

double CovPoisson::_evaluateCovOnSphere(double alpha, double scale) const
{
  double lambda = getParam();
  double valbes = bessel_j(lambda * sin(alpha), 0);
  return exp(lambda * (cos(alpha) - 1.)) * valbes;
}

VectorDouble CovPoisson::_evaluateSpectrumOnSphere(double scale) const
{
  double lambda = getParam();
  int degree = getDegree();
  VectorInt x = VH::sequence(degree+1);
  VectorDouble sp = law_df_poisson_vec(x, lambda);
  if (isFlagNormalizeSpectrum()) VH::normalize(sp, 1);

  return sp;
}
