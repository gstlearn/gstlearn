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
#include "Covariances/KernelPenta.hpp"
#include "Covariances/CovContext.hpp"

namespace gstlrn
{
  KernelPenta::KernelPenta(const CovContext& ctxt)
    : AKernel(ECov::PENTA, ctxt)
  {
  }

  KernelPenta::KernelPenta(const KernelPenta& r)
    : AKernel(r)
  {
  }

  KernelPenta& KernelPenta::operator=(const KernelPenta& r)
  {
    if (this != &r)
    {
      AKernel::operator=(r);
    }
    return *this;
  }

  KernelPenta::~KernelPenta() {}

  double KernelPenta::_evaluateCov(double h) const
  {
    double cov = 0.;
    if (h < 1.)
    {
      cov = 1. - 3. * h * (1. - h / 2. * (1. + h / 6.));
    }
    else if (h < 2.)
    {
      cov = -2. + 3. * h * (1. - h / 2. * (1. - h / 6.));
    }
    return (cov);
  }

} // namespace gstlrn
