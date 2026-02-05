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

class GSTLEARN_EXPORT KernelPoisson: public AKernel
{
public:
  KernelPoisson(const CovContext& ctxt);
  KernelPoisson(const KernelPoisson& r);
  KernelPoisson& operator=(const KernelPoisson& r);
  virtual ~KernelPoisson();

  String getCovName() const override { return "Poisson"; }
  bool hasParam() const override { return true; }
  double getParMax() const override { return MAX_PARAM; }
  Id getMinOrder() const override { return -1; }
  bool getCompatibleSpaceS() const override { return true; }

  bool isValidForSpectralOnRn() const override { return true; }
  bool hasCovOnSphere() const override { return true; }
  bool isValidForSpectralOnSphere() const override { return true; }

protected:
  double _evaluateCovOnSphere(double alpha, double scale = 1., Id degree = 50) const override;
  VectorDouble _evaluateSpectrumOnSphere(Id n, double scale = 1., bool flagScale = true) const override;
};

} // namespace gstlrn
