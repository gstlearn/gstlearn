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

#include "Covariances/ACovFuncWithAutoDiff.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
// Forward declaration
class CovContext;
class TurningBandOperate;

class GSTLEARN_EXPORT CovExponential:
#ifndef SWIG
  public ACovFuncWithAutoDiff<CovExponential>
#else
  public ACovFunc
#endif
{
public:
  CovExponential(const CovContext& ctxt)
    : ACovFuncWithAutoDiff<CovExponential>(ECov::EXPONENTIAL, ctxt)
  {
  }

  CovExponential(const CovExponential& r)
    : ACovFuncWithAutoDiff<CovExponential>(r)
  {
  }

  CovExponential& operator=(const CovExponential& r)
  {
    if (this != &r)
    {
      ACovFuncWithAutoDiff<CovExponential>::operator=(r);
    }
    return *this;
  }
  virtual ~CovExponential();

  String getFormula() const override;
  double getScadef() const override;
  String getCovName() const override { return "Exponential"; }
  Id getMinOrder() const override { return -1; }
  bool getCompatibleSpaceR() const override { return true; }
  bool getCompatibleSpaceS() const override { return true; }
  bool hasCovOnSphere() const override { return true; }
  bool hasSpectrumOnSphere() const override { return true; }

  bool isValidForTurningBand() const override { return true; }
  double simulateTurningBand(double t0, TurningBandOperate& operTB) const override;

  bool isValidForSpectral() const override { return true; }
  MatrixDense simulateSpectralOmega(Id nb) const override;

  template<typename T>
  T evalImpl(T h) const
  {
    if (h < 0) return 0;
    return exp(-h);
  }

protected:
  double _evaluateCovOnSphere(double alpha,
                              double scale = 1.,
                              Id degree    = 50) const override;
  VectorDouble _evaluateSpectrumOnSphere(Id n, double scale = 1.) const override;
  // double _evaluateCovFirstDerivative(double h) const override;
};
} // namespace gstlrn