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

#include "gstlearn_export.hpp"

#include "Covariances/AKernel.hpp"

namespace gstlrn
{

class CovContext;

class GSTLEARN_EXPORT KernelMarkov: public AKernel
{
public:
  KernelMarkov(const CovContext& ctx);
  KernelMarkov(const KernelMarkov& r);
  KernelMarkov& operator=(const KernelMarkov& r);
  virtual ~KernelMarkov();

  String getFormula() const override;
  String getCovName() const override { return "Markov"; }
  bool needCorrec() const override { return true; }

  bool hasParam() const override { return true; }
  double getParMax() const override { return MAX_PARAM; }
  double getScadef() const override;
  bool isValidForSimulation(const ESimuType& simuType) const override
  {
    return (simuType == ESimuType::SPECTRAL);
  }
  bool hasMarkovCoeffs() const override { return true; }
  double normalizeOnSphere(Id n = 50, double scale = 1.) const override;
  double evaluateSpectrum(double freq) const override;
  VectorDouble getMarkovCoeffs() const override { return _markovCoeffs; }
  void setMarkovCoeffs(const VectorDouble& coeffs) override { _markovCoeffs = coeffs; }
  double getCorrec() const override { return _correc; }
  void setCorrec(double val) override;

protected:
  VectorDouble _evaluateSpectrumOnSphere(Id n, double scale = 1., bool flagScale = true) const override;

private:
  VectorDouble _markovCoeffs;
  double _correc;
};
} // namespace gstlrn