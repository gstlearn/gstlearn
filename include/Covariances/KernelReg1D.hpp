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
class CovContext;

class GSTLEARN_EXPORT KernelReg1D: public AKernel
{
public:
  KernelReg1D(const CovContext& ctx);
  KernelReg1D(const KernelReg1D& r);
  KernelReg1D& operator=(const KernelReg1D& r);
  virtual ~KernelReg1D();

  double getScadef() const override;
  String getCovName() const override { return "1-D Regularized"; }
  Id getMinOrder() const override { return -1; }
  bool getCompatibleSpaceR() const override { return true; }

  size_t getMaxNDim() const override { return 1; }

protected:
  double _evaluateCov(double h) const override;
};
} // namespace gstlrn
