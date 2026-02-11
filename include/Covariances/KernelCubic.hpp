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
// Forward declaration
class CovContext;
class TurningBandOperate;

class GSTLEARN_EXPORT KernelCubic: public AKernel
{
public:
  KernelCubic(const CovContext& ctx);
  KernelCubic(const KernelCubic& r);
  KernelCubic& operator=(const KernelCubic& r);
  virtual ~KernelCubic();

  size_t getMaxNDim() const override { return 3; }

  String getFormula() const override;
  String getCovName() const override { return "Cubic"; }
  Id getMinOrder() const override { return -1; }
  bool getCompatibleSpaceR() const override { return true; }
  bool hasCovDerivative() const override { return true; }

  bool isValidForSimulation(const ESimuType& simuType) const override
  {
    return (simuType == ESimuType::TB);
  }
  double simulateTurningBand(double t0, TurningBandOperate& operTB) const override;

protected:
  double _evaluateCov(double h) const override;
  double _evaluateCovDerivative(Id degree, double h) const override;
  double _evaluateCovFirstDerivativeOverH(double h) const override;
};

} // namespace gstlrn