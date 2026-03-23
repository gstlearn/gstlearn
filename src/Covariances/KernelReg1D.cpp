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
#include "Covariances/KernelReg1D.hpp"
#include "Covariances/CovContext.hpp"

namespace gstlrn
{
  KernelReg1D::KernelReg1D(const CovContext& ctxt)
    : AKernel(ECov::REG1D, ctxt)
  {
  }

  KernelReg1D::KernelReg1D(const KernelReg1D& r)
    : AKernel(r)
  {
  }

  KernelReg1D& KernelReg1D::operator=(const KernelReg1D& r)
  {
    if (this != &r)
    {
      AKernel::operator=(r);
    }
    return *this;
  }

  KernelReg1D::~KernelReg1D() {}

  double KernelReg1D::getScadef() const
  {
    return (2.);
  }

  double KernelReg1D::_evaluateCov(double h) const
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
