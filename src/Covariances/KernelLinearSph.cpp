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
#include "Covariances/KernelLinearSph.hpp"
#include "Covariances/CovContext.hpp"

#include <cmath>

namespace gstlrn
{
  KernelLinearSph::KernelLinearSph(const CovContext& ctxt)
    : AKernel(ECov::LINEARSPH, ctxt)
  {
    setParam(1);
  }

  KernelLinearSph::KernelLinearSph(const KernelLinearSph& r)
    : AKernel(r)
  {
  }

  KernelLinearSph& KernelLinearSph::operator=(const KernelLinearSph& r)
  {
    if (this != &r)
    {
      AKernel::operator=(r);
    }
    return *this;
  }

  KernelLinearSph::~KernelLinearSph() {}

  double KernelLinearSph::_evaluateCovOnSphere(double alpha,
                                               double scale,
                                               Id degree) const
  {
    DECLARE_UNUSED(scale);
    DECLARE_UNUSED(degree);
    return 1. - 2. * alpha / GV_PI;
  }

  VectorDouble KernelLinearSph::_evaluateSpectrumOnSphere(Id n,
                                                          double scale,
                                                          bool flagScale) const
  {
    DECLARE_UNUSED(scale);
    VectorDouble sp(n + 1, 0.);

    Id k = 1;
    sp[k] = 3. / 4.;
    while (1)
    {
      k += 2;
      if (k >= n + 1) break;
      double v = (k - 2.) / (k + 1.);
      sp[k] = (2. * k + 1.) / (2. * k - 3.) * v * v * sp[k - 2];
    }

    if (flagScale) sp.normalizeInPlace(1);

    return sp;
  }
} // namespace gstlrn
