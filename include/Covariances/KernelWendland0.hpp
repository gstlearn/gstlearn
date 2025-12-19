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

// In Piecewise polynomial, positive definite and compactly supported
// radial functions of minimal degree, by H. Wendland
// Advances in Computational Mathematics, Vol. 4 (389-396), 1995
// It corresponds to Wendland \psi_{2,0}

namespace gstlrn
{
class CovContext;

class GSTLEARN_EXPORT KernelWendland0: public AKernel
{
public:
  KernelWendland0(const CovContext& ctx);
  KernelWendland0(const KernelWendland0& r);
  KernelWendland0& operator=(const KernelWendland0& r);
  virtual ~KernelWendland0();

  size_t getMaxNDim() const override { return 3; }

  String getCovName() const override { return "Wendland-2,0"; }
  Id getMinOrder() const override { return -1; }
  bool getCompatibleSpaceR() const override { return true; }

protected:
  double _evaluateCov(double h) const override;
};

} // namespace gstlrn