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
#include "Covariances/KernelStorkey.hpp"
#include "Covariances/CovContext.hpp"

#include <cmath>

namespace gstlrn
{
  KernelStorkey::KernelStorkey(const CovContext& ctxt)
    : AKernel(ECov::STORKEY, ctxt)
  {
  }

  KernelStorkey::KernelStorkey(const KernelStorkey& r)
    : AKernel(r)
  {
  }

  KernelStorkey& KernelStorkey::operator=(const KernelStorkey& r)
  {
    if (this != &r)
    {
      AKernel::operator=(r);
    }
    return *this;
  }

  KernelStorkey::~KernelStorkey() {}

  double KernelStorkey::_evaluateCov(double h) const
  {
    double cov = 0.;
    double pi2 = 2. * GV_PI;
    if (h < 1)
      cov = (2. * (1. - h) * (1. + cos(pi2 * h) / 2.) + 3 / pi2 * sin(pi2 * h))
          / 3.;
    return (cov);
  }

} // namespace gstlrn
