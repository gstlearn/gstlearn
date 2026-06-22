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
#include "Covariances/KernelWendland0.hpp"
#include "Covariances/CovContext.hpp"

namespace gstlrn
{
  KernelWendland0::KernelWendland0(const CovContext& ctxt)
    : AKernel(ECov::WENDLAND0, ctxt)
  {
  }

  KernelWendland0::KernelWendland0(const KernelWendland0& r)
    : AKernel(r)
  {
  }

  KernelWendland0& KernelWendland0::operator=(const KernelWendland0& r)
  {
    if (this != &r)
    {
      AKernel::operator=(r);
    }
    return *this;
  }

  KernelWendland0::~KernelWendland0() {}

  double KernelWendland0::_evaluateCov(double h) const
  {
    double cov = 0.;
    double h2 = h * h;
    if (h < 1) cov = 1. - 2 * h + h2;
    return (cov);
  }

} // namespace gstlrn
