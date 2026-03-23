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
#include "Covariances/KernelGeometric.hpp"
#include "Covariances/CovContext.hpp"

#include <cmath>

namespace gstlrn
{
  KernelGeometric::KernelGeometric(const CovContext& ctxt)
    : AKernel(ECov::GEOMETRIC, ctxt)
  {
    setParam(1);
  }

  KernelGeometric::KernelGeometric(const KernelGeometric& r)
    : AKernel(r)
  {
  }

  KernelGeometric& KernelGeometric::operator=(const KernelGeometric& r)
  {
    if (this != &r)
    {
      AKernel::operator=(r);
    }
    return *this;
  }

  KernelGeometric::~KernelGeometric() {}

  double KernelGeometric::_evaluateCovOnSphere(double alpha,
                                               double scale,
                                               Id degree) const
  {
    DECLARE_UNUSED(degree);
    double rho = scale;
    return ((1. - rho) / sqrt(1. - 2. * rho * cos(alpha) + rho * rho));
  }

  VectorDouble KernelGeometric::_evaluateSpectrumOnSphere(Id n,
                                                          double scale,
                                                          bool flagScale) const
  {
    double rho = scale;
    VectorDouble sp(1 + n, 0.);

    double rhoprod = 1.;
    for (Id k = 0; k <= n; k++)
    {
      sp[k] = rhoprod;
      rhoprod *= rho;
    }

    if (flagScale) sp.normalizeInPlace(1);

    return sp;
  }
} // namespace gstlrn
