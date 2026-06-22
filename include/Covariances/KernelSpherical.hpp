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
#include "Enum/ESimuType.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
  class CovContext;
  class TurningBandOperate;

  class GSTLEARN_EXPORT KernelSpherical: public AKernel
  {
  public:
    KernelSpherical(const CovContext& ctx);
    KernelSpherical(const KernelSpherical& r);
    KernelSpherical& operator=(const KernelSpherical& r);
    virtual ~KernelSpherical();

    size_t getMaxNDim() const override { return 3; }

    String getFormula() const override;

    String getCovName() const override { return "Spherical"; }

    Id getMinOrder() const override { return -1; }

    bool getCompatibleSpaceR() const override { return true; }

    bool isValidForSimulation(const ESimuType& simuType) const override
    {
      return (simuType == ESimuType::TB);
    }

    double
      simulateTurningBand(double t0, TurningBandOperate& operTB) const override;

  protected:
    double _evaluateCov(double h) const override;
    double _evaluateCovFirstDerivative(double h) const override;
  };

} // namespace gstlrn
