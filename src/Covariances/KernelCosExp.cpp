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
#include "Covariances/KernelCosExp.hpp"
#include "Covariances/CovContext.hpp"

#include <cmath>

namespace gstlrn
{
  KernelCosExp::KernelCosExp(const CovContext& ctxt)
    : AKernel(ECov::COSEXP, ctxt)
  {
    setParam(1);
  }

  KernelCosExp::KernelCosExp(const KernelCosExp& r)
    : AKernel(r)
  {
  }

  KernelCosExp& KernelCosExp::operator=(const KernelCosExp& r)
  {
    if (this != &r)
    {
      AKernel::operator=(r);
    }
    return *this;
  }

  KernelCosExp::~KernelCosExp() {}

  double KernelCosExp::getScadef() const
  {
    return (2.995732);
  }

  double KernelCosExp::_evaluateCov(double h) const
  {
    double cov = 1.;
    if (h > 100) return (0.);
    cov = exp(-h);
    double h2 = h / getParam();
    cov *= cos(2. * GV_PI * h2);
    return (cov);
  }
} // namespace gstlrn
