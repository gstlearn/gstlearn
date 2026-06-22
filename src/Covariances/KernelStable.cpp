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
#include "Covariances/KernelStable.hpp"
#include "Covariances/CovContext.hpp"
#include "Simulation/TurningBandOperate.hpp"

#include <cmath>

namespace gstlrn
{
  KernelStable::KernelStable(const CovContext& ctxt)
    : AKernel(ECov::STABLE, ctxt)
  {
    setParam(1);
  }

  KernelStable::KernelStable(const KernelStable& r)
    : AKernel(r)
  {
  }

  KernelStable& KernelStable::operator=(const KernelStable& r)
  {
    if (this != &r)
    {
      AKernel::operator=(r);
    }
    return *this;
  }

  KernelStable::~KernelStable() {}

  double KernelStable::getScadef() const
  {
    return pow(3., 1. / getParam());
  }

  double KernelStable::_evaluateCov(double h) const
  {
    double cov = 1.;
    if (h > 0) cov = exp(-pow(h, getParam()));

    return (cov);
  }

  double KernelStable::simulateTurningBand(
    double t0,
    TurningBandOperate& operTB) const
  {
    if (getParam() > 1) return operTB.cosineOne(t0);
    return operTB.spectralOne(t0);
  }
} // namespace gstlrn
