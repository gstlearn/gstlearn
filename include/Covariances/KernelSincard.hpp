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
class CovContext;
class TurningBandOperate;

class GSTLEARN_EXPORT KernelSincard : public AKernel
{
public:
  KernelSincard(const CovContext& ctx);
  KernelSincard(const KernelSincard &r);
  KernelSincard& operator= (const KernelSincard &r);
  virtual ~KernelSincard();

  double         getScadef() const override;
  String getFormula() const override;
  String         getCovName() const override { return "Cardinal Sine"; }
  Id            getMinOrder() const override { return -1; }
  bool           getCompatibleSpaceR() const override { return true; }

  bool isValidForTurningBand() const override { return true; }
  double simulateTurningBand(double t0, TurningBandOperate &operTB) const override;

protected:
  double _evaluateCov(double h) const override;
};

}