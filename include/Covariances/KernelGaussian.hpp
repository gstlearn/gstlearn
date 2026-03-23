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
  class TurningBandOperate;
  class MatrixDense;

  class GSTLEARN_EXPORT KernelGaussian: public AKernel
  {
  public:
    KernelGaussian(const CovContext& ctx);
    KernelGaussian(const KernelGaussian& r);
    KernelGaussian& operator=(const KernelGaussian& r);
    virtual ~KernelGaussian();

    String getFormula() const override;

    String getCovName() const override { return "Gaussian"; }

    Id getMinOrder() const override { return -1; }

    double getScadef() const override;

    bool hasCovDerivative() const override { return true; }

    bool getCompatibleSpaceR() const override { return true; }

    double
      simulateTurningBand(double t0, TurningBandOperate& operTB) const override;

    bool isValidForSimulation(const ESimuType& simuType) const override
    {
      return (
        (getSpaceType() == ESpaceType::RN && simuType == ESimuType::SPECTRAL)
        || simuType == ESimuType::TB);
    }

    double evaluateSpectrum(double freq) const override;
    MatrixDense simulateSpectralOmega(Id nb) const override;

  protected:
    double _evaluateCov(double h) const override;
    double _evaluateCovDerivative(Id degree, double h) const override;
    double _evaluateCovFirstDerivativeOverH(double h) const override;
  };

} // namespace gstlrn
