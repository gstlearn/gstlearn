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
#include "Covariances/KernelGC1.hpp"
#include "Covariances/CovContext.hpp"
#include "Simulation/TurningBandOperate.hpp"

namespace gstlrn
{
KernelGC1::KernelGC1(const CovContext& ctxt)
  : AKernel(ECov::ORDER1_GC, ctxt)
{
}

KernelGC1::KernelGC1(const KernelGC1& r)
  : AKernel(r)
{
}

KernelGC1& KernelGC1::operator=(const KernelGC1& r)
{
  if (this != &r)
  {
    AKernel::operator=(r);
  }
  return *this;
}

KernelGC1::~KernelGC1()
{
}

double KernelGC1::_evaluateCov(double h) const
{
  double cov;
  double r = getContext().getField();
  auto ndim = getContext().getNDim();

  if (ndim == 1)
    cov = r - h;
  else if (ndim == 2)
    cov = r * GV_PI / 2 - h;
  else
    cov = r * 2 - h;

  return (cov);
}

double KernelGC1::simulateTurningBand(double t0, TurningBandOperate& operTB) const
{
  return operTB.IRFProcessOne(t0);
}
}