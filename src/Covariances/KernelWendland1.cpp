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
#include "Covariances/KernelWendland1.hpp"
#include "Covariances/CovContext.hpp"

namespace gstlrn
{
KernelWendland1::KernelWendland1(const CovContext& ctxt)
  : AKernel(ECov::WENDLAND1, ctxt)
{
}

KernelWendland1::KernelWendland1(const KernelWendland1& r)
  : AKernel(r)
{
}

KernelWendland1& KernelWendland1::operator=(const KernelWendland1& r)
{
  if (this != &r)
  {
    AKernel::operator=(r);
  }
  return *this;
}

KernelWendland1::~KernelWendland1()
{
}

double KernelWendland1::_evaluateCov(double h) const
{
  // From "Computed Supported Correlation Functions" by T. Gneiting with n=3
  double cov = 0.;
  if (h < 1)
    cov = (1. + 4. * h) * pow(1. - h, 4.);
  return (cov);
}

}