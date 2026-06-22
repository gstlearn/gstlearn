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

  class GSTLEARN_EXPORT KernelCosExp: public AKernel
  {
  public:
    KernelCosExp(const CovContext& ctx);
    KernelCosExp(const KernelCosExp& r);
    KernelCosExp& operator=(const KernelCosExp& r);
    virtual ~KernelCosExp();

    double getParMax() const override { return TEST; }

    bool hasParam() const override { return true; }

    double getScadef() const override;

    String getCovName() const override { return "Cosexp"; }

    Id getMinOrder() const override { return -1; }

    bool getCompatibleSpaceR() const override { return true; }

  protected:
    double _evaluateCov(double h) const override;
  };
} // namespace gstlrn
