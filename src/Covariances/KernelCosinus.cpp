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
#include "Covariances/KernelCosinus.hpp"
#include "Covariances/CovContext.hpp"

#include <cmath>

namespace gstlrn
{
KernelCosinus::KernelCosinus(const CovContext& ctxt)
  : AKernel(ECov::COSINUS, ctxt)
{
}

KernelCosinus::KernelCosinus(const KernelCosinus& r)
  : AKernel(r)
{
}

KernelCosinus& KernelCosinus::operator=(const KernelCosinus& r)
{
  if (this != &r)
  {
    AKernel::operator=(r);
  }
  return *this;
}

KernelCosinus::~KernelCosinus()
{
}

double KernelCosinus::_evaluateCov(double h) const
{
  double cov = cos(2. * GV_PI * h);
  return (cov);
}
} // namespace gstlrn