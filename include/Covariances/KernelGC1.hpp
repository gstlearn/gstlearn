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
#include "Covariances/CovContext.hpp"
#include "gstlearn_export.hpp"

/* Be careful ! This is not a real covariance */

namespace gstlrn
{
class CovContext;
class TurningBandOperate;

class GSTLEARN_EXPORT KernelGC1: public AKernel
{
public:
  KernelGC1(const CovContext& ctx);
  KernelGC1(const KernelGC1& r);
  KernelGC1& operator=(const KernelGC1& r);
  virtual ~KernelGC1();

  Id hasRange() const override { return -1; }
  Id getMinOrder() const override { return 0; }
  String getCovName() const override { return "Order-1 G.C."; }
  bool getCompatibleSpaceR() const override { return true; }

  bool isValidForSimulation(const ESimuType& simuType) const override
  {
    return (simuType == ESimuType::TB);
  }

  double simulateTurningBand(double t0, TurningBandOperate& operTB) const override;

protected:
  double _evaluateCov(double h) const override;
};

} // namespace gstlrn