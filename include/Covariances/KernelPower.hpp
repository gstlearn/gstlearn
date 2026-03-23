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

/* Be careful ! This is not a real covariance */

namespace gstlrn
{
  class CovContext;
  class TurningBandOperate;

  class GSTLEARN_EXPORT KernelPower: public AKernel
  {
  public:
    KernelPower(const CovContext& ctx);
    KernelPower(const KernelPower& r);
    KernelPower& operator=(const KernelPower& r);
    virtual ~KernelPower();

    Id hasRange() const override { return -1; }

    bool hasParam() const override { return true; }

    double getParMax() const override { return 1.99; }

    Id getMinOrder() const override { return 0; }

    String getCovName() const override { return "Power"; }

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
