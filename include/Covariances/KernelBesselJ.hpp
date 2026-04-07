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

  class GSTLEARN_EXPORT KernelBesselJ: public AKernel
  {
  public:
    KernelBesselJ(const CovContext& ctx);
    KernelBesselJ(const KernelBesselJ& r);
    KernelBesselJ& operator=(const KernelBesselJ& r);
    virtual ~KernelBesselJ();

    bool hasParam() const override { return true; }

    double getParMax() const override { return 2; }

    String getFormula() const override;

    String getCovName() const override { return "J-Bessel"; }

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
  };
} // namespace gstlrn
