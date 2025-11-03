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

/* Be careful ! This is not a real covariance */

namespace gstlrn
{
class CovContext;
class TurningBandOperate;

class GSTLEARN_EXPORT KernelLinear : public AKernel
{
public:
  KernelLinear(const CovContext& ctx);
  KernelLinear(const KernelLinear &r);
  KernelLinear& operator= (const KernelLinear &r);
  virtual ~KernelLinear();

  Id    hasRange() const override { return -1; }
  Id    getMinOrder()  const override { return 0; }
  String getCovName() const override { return "Linear"; }
  bool   getCompatibleSpaceR() const override { return true; }

  bool isValidForTurningBand() const override { return true; }
  double simulateTurningBand(double t0, TurningBandOperate &operTB) const override;

protected:
  double _evaluateCov(double h) const override;
};

}