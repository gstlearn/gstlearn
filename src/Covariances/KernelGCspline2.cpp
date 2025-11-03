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
#include "Covariances/KernelGCspline2.hpp"
#include "Covariances/CovContext.hpp"

#include <cmath>

namespace gstlrn
{
KernelGCspline2::KernelGCspline2(const CovContext& ctxt)
  : AKernel(ECov::SPLINE2_GC, ctxt)
{
}

KernelGCspline2::KernelGCspline2(const KernelGCspline2& r)
  : AKernel(r)
{
}

KernelGCspline2& KernelGCspline2::operator=(const KernelGCspline2& r)
{
  if (this != &r)
  {
    AKernel::operator=(r);
  }
  return *this;
}

KernelGCspline2::~KernelGCspline2()
{
}

double KernelGCspline2::_evaluateCov(double h) const
{
  double B = 1.;
  double A = (7. - 10. * B) / 12.;
  double C = (-7. - 2. * B) / 12.;

  double h2     = h * h;
  double logval = (h < 10.e-5) ? 0. : log(h);
  double cov    = A + h2 * (B + h2 * (C + logval));

  return (cov);
}

double KernelGCspline2::_evaluateCovFirstDerivativeOverH(double h) const
{
  double B      = 1.;
  double C      = (-7. - 2. * B) / 12.;
  double h2     = h * h;
  double logval = (h < EPSILON5) ? 0. : log(h);

  double cov = 2. * B + h2 * (1. + 4. * (C + logval));
  return cov;
}

double KernelGCspline2::_evaluateCovDerivative(Id degree, double h) const
{
  double B      = 1.;
  double C      = (-7. - 2. * B) / 12.;
  double h2     = h * h;
  double logval = (h < EPSILON5) ? 0. : log(h);

  double cov = 0.;
  switch (degree)
  {
    case 1:
      cov = 2. * B + h2 * (1. + 4. * (C + logval));
      break;

    case 2:
      cov = 2. * B + h2 * (7. + 12. * (C + logval));
      break;
  }

  return (cov);
}
} // namespace gstlrn