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
#include "Covariances/KernelLinear.hpp"
#include "Covariances/CovContext.hpp"
#include "Simulation/TurningBandOperate.hpp"

namespace gstlrn
{
  KernelLinear::KernelLinear(const CovContext& ctxt)
    : AKernel(ECov::LINEAR, ctxt)
  {
  }

  KernelLinear::KernelLinear(const KernelLinear& r)
    : AKernel(r)
  {
  }

  KernelLinear& KernelLinear::operator=(const KernelLinear& r)
  {
    if (this != &r)
    {
      AKernel::operator=(r);
    }
    return *this;
  }

  KernelLinear::~KernelLinear() {}

  double KernelLinear::_evaluateCov(double h) const
  {
    double cov, r;

    r = getContext().getField();

    auto ndim = getContext().getNDim();
    if (ndim == 1)
      cov = r - h;
    else if (ndim == 2)
      cov = r * GV_PI / 2 - h;
    else
      cov = r * 2 - h;

    return (cov);
  }

  double KernelLinear::simulateTurningBand(double t0,
                                           TurningBandOperate& operTB) const
  {
    return operTB.IRFProcessOne(t0);
  }
} // namespace gstlrn
