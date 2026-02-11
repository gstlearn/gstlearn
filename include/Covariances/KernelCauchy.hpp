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

class GSTLEARN_EXPORT KernelCauchy: public AKernel
{
public:
  KernelCauchy(const CovContext& ctx);
  KernelCauchy(const KernelCauchy& r);
  KernelCauchy& operator=(const KernelCauchy& r);
  virtual ~KernelCauchy();

  String getFormula() const override;
  String getCovName() const override { return "Cauchy"; }
  Id getMinOrder() const override { return -1; }
  bool getCompatibleSpaceR() const override { return true; }

  bool hasParam() const override { return true; }
  double getParMax() const override { return MAX_PARAM; }
  double getScadef() const override;

  bool isValidForSimulation(const ESimuType& simuType) const override
  {
    return (getSpaceType() == ESpaceType::RN && simuType == ESimuType::SPECTRAL);
  }
  MatrixDense simulateSpectralOmega(Id nb) const override;

protected:
  double _evaluateCov(double h) const override;
};

} // namespace gstlrn