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

  class GSTLEARN_EXPORT KernelTriangle: public AKernel
  {
  public:
    KernelTriangle(const CovContext& ctx);
    KernelTriangle(const KernelTriangle& r);
    KernelTriangle& operator=(const KernelTriangle& r);
    virtual ~KernelTriangle();

    size_t getMaxNDim() const override { return 1; }

    String getCovName() const override { return "Triangle"; }

    Id getMinOrder() const override { return -1; }

    bool getCompatibleSpaceR() const override { return true; }

  protected:
    double _evaluateCov(double h) const override;
  };
} // namespace gstlrn
