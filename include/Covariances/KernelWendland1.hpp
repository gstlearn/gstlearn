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

// In Piecewise polynomial, positive definite and compactly supported
// radial functions of minimal degree, by H. Wendland
// Advances in Computational Mathematics, Vol. 4 (389-396), 1995
// It corresponds to Wendland \psi_{3,1}

namespace gstlrn
{
class CovContext;

class GSTLEARN_EXPORT KernelWendland1 : public AKernel
{
public:
  KernelWendland1(const CovContext& ctx);
  KernelWendland1(const KernelWendland1 &r);
  KernelWendland1& operator= (const KernelWendland1 &r);
  virtual ~KernelWendland1();

  size_t getMaxNDim() const override { return 3; }

  String         getCovName() const override { return "Wendland-3,1"; }
  Id            getMinOrder() const override { return -1; }
  bool           getCompatibleSpaceR() const override { return true; }

protected:
  double _evaluateCov(double h) const override;
};

}