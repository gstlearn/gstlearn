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
#include "Covariances/KernelGCspline.hpp"
#include "Covariances/CovContext.hpp"
#include "Simulation/TurningBandOperate.hpp"

#include <cmath>

namespace gstlrn
{
  KernelGCspline::KernelGCspline(const CovContext& ctxt)
    : AKernel(ECov::SPLINE_GC, ctxt)
  {
  }

  KernelGCspline::KernelGCspline(const KernelGCspline& r)
    : AKernel(r)
  {
  }

  KernelGCspline& KernelGCspline::operator=(const KernelGCspline& r)
  {
    if (this != &r)
    {
      AKernel::operator=(r);
    }
    return *this;
  }

  KernelGCspline::~KernelGCspline() {}

  double KernelGCspline::_evaluateCov(double h) const
  {
    auto ndim = getContext().getNDim();
    double r = getContext().getField();
    double r2 = r * r;
    double h2 = h * h;
    double logval = (r < 10.e-5 || h < r * 10.e-5) ? 0. : log(h / r);

    // Code for the first 3 Space dimensions
    double cov = 0.;
    if (ndim == 1)
      cov = 0.5 * r2 - h2 * (1.5 - log(2) - logval);
    else if (ndim == 2)
      cov = r2 - h2 * (1 - logval);
    else
      cov = 1.5 * r2 - h2 * (11. / 6. - log(2) - logval);

    return (cov);
  }

  double KernelGCspline::simulateTurningBand(double t0,
                                             TurningBandOperate& operTB) const
  {
    return operTB.cosineOne(t0);
  }
} // namespace gstlrn
