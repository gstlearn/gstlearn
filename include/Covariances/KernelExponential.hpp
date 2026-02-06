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

#include "Covariances/AKernelWithAutoDiff.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
// Forward declaration
class CovContext;
class TurningBandOperate;

class GSTLEARN_EXPORT KernelExponential:
#ifndef SWIG
  public AKernelWithAutoDiff<KernelExponential>
#else
  public AKernel
#endif
{
public:
  KernelExponential(const CovContext& ctxt)
    : AKernelWithAutoDiff<KernelExponential>(ECov::EXPONENTIAL, ctxt)
  {
  }

  KernelExponential(const KernelExponential& r)
    : AKernelWithAutoDiff<KernelExponential>(r)
  {
  }

  KernelExponential& operator=(const KernelExponential& r)
  {
    if (this != &r)
    {
      AKernelWithAutoDiff<KernelExponential>::operator=(r);
    }
    return *this;
  }
  virtual ~KernelExponential();

  String getFormula() const override;
  double getScadef() const override;
  String getCovName() const override { return "Exponential"; }
  Id getMinOrder() const override { return -1; }
  bool getCompatibleSpaceR() const override { return true; }
  bool getCompatibleSpaceS() const override { return true; }

  bool isValidForTurningBand() const override { return true; }
  double simulateTurningBand(double t0, TurningBandOperate& operTB) const override;

  bool isValidForSpectralOnRn() const override { return true; }
  bool hasCovOnSphere() const override { return true; }
  bool isValidForSpectralOnSphere() const override { return true; }
  double evaluateSpectrum(double freq) const override;

  MatrixDense simulateSpectralOmega(Id nb) const override;
  template<typename T>
  T evalImpl(T h) const
  {
    if (h < 0) return 0;
    return exp(-h);
  }

protected:
  double _evaluateCovOnSphere(double alpha, double scale = 1., Id degree = 50) const override;
  VectorDouble _evaluateSpectrumOnSphere(Id n, double scale = 1., bool flagScale = true) const override;
  // double _evaluateCovFirstDerivative(double h) const override;
};
} // namespace gstlrn