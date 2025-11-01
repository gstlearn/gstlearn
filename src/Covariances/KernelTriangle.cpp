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
#include "Covariances/KernelTriangle.hpp"
#include "Covariances/CovContext.hpp"

namespace gstlrn
{
KernelTriangle::KernelTriangle(const CovContext& ctxt)
  : AKernel(ECov::TRIANGLE, ctxt)
{
}

KernelTriangle::KernelTriangle(const KernelTriangle& r)
  : AKernel(r)
{
}

KernelTriangle& KernelTriangle::operator=(const KernelTriangle& r)
{
  if (this != &r)
  {
    AKernel::operator=(r);
  }
  return *this;
}

KernelTriangle::~KernelTriangle()
{
}

double KernelTriangle::_evaluateCov(double h) const
{
  double cov = MAX(0, 1. - h);
  return (cov);
}

}