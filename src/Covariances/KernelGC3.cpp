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
#include "Covariances/KernelGC3.hpp"
#include "Covariances/CovContext.hpp"
#include "Simulation/TurningBandOperate.hpp"

namespace gstlrn
{
  KernelGC3::KernelGC3(const CovContext& ctxt)
    : AKernel(ECov::ORDER3_GC, ctxt)
  {
  }

  KernelGC3::KernelGC3(const KernelGC3& r)
    : AKernel(r)
  {
  }

  KernelGC3& KernelGC3::operator=(const KernelGC3& r)
  {
    if (this != &r)
    {
      AKernel::operator=(r);
    }
    return *this;
  }

  KernelGC3::~KernelGC3() {}

  double KernelGC3::_evaluateCov(double h) const
  {
    double cov;
    double r = getContext().getField();
    auto ndim = getContext().getNDim();

    double h2 = h * h;
    double r3 = r * r * r;

    if (ndim == 1)
      cov = h2 * (h - 3. * r) + 2. * r3;
    else if (ndim == 2)
      cov = h2 * (h - 9. * GV_PI * r / 8. * r) + 3. * GV_PI * r3 / 2;
    else
      cov = h2 * (h - 4. * r) + 8. * r3;

    return (cov);
  }

  double
    KernelGC3::simulateTurningBand(double t0, TurningBandOperate& operTB) const
  {
    return operTB.IRFProcessOne(t0);
  }
} // namespace gstlrn
