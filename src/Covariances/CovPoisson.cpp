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
#include "Basic/Law.hpp"
#include "Basic/MathFunc.hpp"
#include "Basic/VectorHelper.hpp"
#include "Covariances/CovContext.hpp"

#include <cmath>

namespace gstlrn
{
CovPoisson::CovPoisson(const CovContext& ctxt)
  : ACovFunc(ECov::POISSON, ctxt)
{
  setParam(1);
}

CovPoisson::CovPoisson(const CovPoisson& r)
  : ACovFunc(r)
{
}

CovPoisson& CovPoisson::operator=(const CovPoisson& r)
{
  if (this != &r)
  {
    ACovFunc::operator=(r);
  }
  return *this;
}

CovPoisson::~CovPoisson()
{
}

double CovPoisson::_evaluateCovOnSphere(double alpha,
                                        double scale,
                                        Id degree) const
{
  DECLARE_UNUSED(scale);
  DECLARE_UNUSED(degree);
  double lambda = getParam();
  double valbes = besselj(lambda * sin(alpha), 0);
  return exp(lambda * (cos(alpha) - 1.)) * valbes;
}

VectorDouble CovPoisson::_evaluateSpectrumOnSphere(Id n, double scale) const
{
  DECLARE_UNUSED(scale);
  double lambda   = getParam();
  VectorInt x     = VH::sequence(n + 1);
  VectorDouble sp = law_df_poisson_vec(x, lambda);
  VH::normalize(sp, 1);

  return sp;
}
}
