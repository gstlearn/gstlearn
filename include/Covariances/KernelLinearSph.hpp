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
  // Forward declaration
  class CovContext;

  class GSTLEARN_EXPORT KernelLinearSph: public AKernel
  {
  public:
    KernelLinearSph(const CovContext& ctx);
    KernelLinearSph(const KernelLinearSph& r);
    KernelLinearSph& operator=(const KernelLinearSph& r);
    virtual ~KernelLinearSph();

    String getCovName() const override { return "LinearSph"; }

    Id getMinOrder() const override { return -1; }

    bool getCompatibleSpaceS() const override { return true; }

    bool isValidForSimulation(const ESimuType& simuType) const override
    {
      return (
        getSpaceType() == ESpaceType::SN && simuType == ESimuType::SPECTRAL);
    }

  protected:
    double _evaluateCovOnSphere(double alpha, double scale = 1., Id degree = 50)
      const override;
    VectorDouble
      _evaluateSpectrumOnSphere(Id n, double scale = 1., bool flagScale = true)
        const override;
  };

} // namespace gstlrn
