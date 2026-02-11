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
class TurningBandOperate;

class GSTLEARN_EXPORT KernelGCspline: public AKernel
{
public:
  KernelGCspline(const CovContext& ctx);
  KernelGCspline(const KernelGCspline& r);
  KernelGCspline& operator=(const KernelGCspline& r);
  virtual ~KernelGCspline();

  Id hasRange() const override { return -1; }
  String getCovName() const override { return "Spline G.C."; }
  Id getMinOrder() const override { return 1; }
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