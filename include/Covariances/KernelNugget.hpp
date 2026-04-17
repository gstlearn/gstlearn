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

  class GSTLEARN_EXPORT KernelNugget: public AKernel
  {
  public:
    KernelNugget(const CovContext& ctx);
    KernelNugget(const KernelNugget& r);
    KernelNugget& operator=(const KernelNugget& r);
    virtual ~KernelNugget();

    String getFormula() const override;

    String getCovName() const override { return "Nugget Effect"; }

    Id getMinOrder() const override { return -1; }

    bool getCompatibleSpaceR() const override { return true; }

    Id hasRange() const override { return 0; }

    bool isValidForSimulation(const ESimuType& simuType) const override
    {
      return (simuType == ESimuType::TB || simuType == ESimuType::SPECTRAL);
    }

  protected:
    double _evaluateCov(double h) const override;
  };

} // namespace gstlrn
