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
#pragma once

#include "Covariances/AKernel.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
/* Be careful ! This is not a real covariance */

class CovContext;

class GSTLEARN_EXPORT KernelGCspline2: public AKernel
{
public:
  KernelGCspline2(const CovContext& ctx);
  KernelGCspline2(const KernelGCspline2& r);
  KernelGCspline2& operator=(const KernelGCspline2& r);
  virtual ~KernelGCspline2();

  Id getMinOrder() const override { return 1; }
  size_t getMaxNDim() const override { return 3; }
  String getCovName() const override { return "Spline-2 G.C."; }
  bool getCompatibleSpaceR() const override { return true; }
  bool hasCovDerivative() const override { return true; }

protected:
  double _evaluateCov(double h) const override;
  double _evaluateCovDerivative(Id degree, double h) const override;
  double _evaluateCovFirstDerivativeOverH(double h) const override;
};

} // namespace gstlrn