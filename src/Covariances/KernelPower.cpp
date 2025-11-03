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
#include "Covariances/KernelPower.hpp"
#include "Basic/MathFunc.hpp"
#include "Covariances/CovContext.hpp"
#include "Simulation/TurningBandOperate.hpp"

#include <cmath>

namespace gstlrn
{
KernelPower::KernelPower(const CovContext& ctxt)
  : AKernel(ECov::POWER, ctxt)
{
  setParam(1);
}

KernelPower::KernelPower(const KernelPower& r)
  : AKernel(r)
{
}

KernelPower& KernelPower::operator=(const KernelPower& r)
{
  if (this != &r)
  {
    AKernel::operator=(r);
  }
  return *this;
}

KernelPower::~KernelPower()
{
}

double KernelPower::_evaluateCov(double h) const
{
  double a;

  auto ndim    = getContext().getNDim();
  double alpha = getParam();
  double r     = getContext().getField();
  double ra    = pow(r, alpha);
  double g1    = loggamma((alpha + 1.) / 2.) + loggamma(1. - alpha / 2.);

  if (ndim == 1)
    a = exp(g1) * ra / sqrt(GV_PI);
  else if (ndim == 2)
    a = (exp(loggamma(alpha + 1.5) - loggamma(alpha + 1.) + g1) * 2. * ra / GV_PI);
  else
    a = (alpha + 1.) * exp(g1) * ra / sqrt(GV_PI);

  double cov = a;
  if (h > 0) cov -= pow(h, getParam());
  return (cov);
}

double KernelPower::simulateTurningBand(double t0, TurningBandOperate& operTB) const
{
  return operTB.cosineOne(t0);
}
}