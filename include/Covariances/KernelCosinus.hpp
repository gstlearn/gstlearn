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

#include "gstlearn_export.hpp"
#include "Covariances/AKernel.hpp"

namespace gstlrn
{
  // Forward declaration
class CovContext;

class GSTLEARN_EXPORT KernelCosinus : public AKernel
{
public:
  KernelCosinus(const CovContext& ctx);
  KernelCosinus(const KernelCosinus &r);
  KernelCosinus& operator= (const KernelCosinus &r);
  virtual ~KernelCosinus();

  size_t getMaxNDim() const override { return 1; }

  String         getCovName() const override { return "Cosinus"; }
  Id            getMinOrder() const override { return -1; }
  bool           getCompatibleSpaceR() const override { return true; }

protected:
  double _evaluateCov(double h) const override;
};

}