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
#include "Covariances/KernelPoisson.hpp"
#include "Basic/Law.hpp"
#include "Basic/MathFunc.hpp"
#include "Basic/VectorHelper.hpp"
#include "Covariances/CovContext.hpp"

#include <cmath>

namespace gstlrn
{
KernelPoisson::KernelPoisson(const CovContext& ctxt)
  : AKernel(ECov::POISSON, ctxt)
{
  setParam(1);
}

KernelPoisson::KernelPoisson(const KernelPoisson& r)
  : AKernel(r)
{
}

KernelPoisson& KernelPoisson::operator=(const KernelPoisson& r)
{
  if (this != &r)
  {
    AKernel::operator=(r);
  }
  return *this;
}

KernelPoisson::~KernelPoisson()
{
}

double KernelPoisson::_evaluateCovOnSphere(double alpha,
                                        double scale,
                                        Id degree) const
{
  DECLARE_UNUSED(scale);
  DECLARE_UNUSED(degree);
  double lambda = getParam();
  double valbes = besselj(lambda * sin(alpha), 0);
  return exp(lambda * (cos(alpha) - 1.)) * valbes;
}

VectorDouble KernelPoisson::_evaluateSpectrumOnSphere(Id n, double scale) const
{
  DECLARE_UNUSED(scale);
  double lambda   = getParam();
  VectorInt x     = VH::sequence(n + 1);
  VectorDouble sp = law_df_poisson_vec(x, lambda);
  sp.normalizeInPlace(1);

  return sp;
}
}
