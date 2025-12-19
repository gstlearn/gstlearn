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

class GSTLEARN_EXPORT KernelCauchyGen: public AKernel
{
public:
  KernelCauchyGen(const CovContext& ctx);
  KernelCauchyGen(const KernelCauchyGen& r);
  KernelCauchyGen& operator=(const KernelCauchyGen& r);
  virtual ~KernelCauchyGen();

  String getFormula() const override;
  String getCovName() const override { return "Generalized Cauchy"; };
  Id getMinOrder() const override { return -1; };
  bool getCompatibleSpaceR() const override { return true; };

  bool hasParam() const override { return true; };
  double getParMax() const override { return MAX_PARAM; };
  double getScadef() const override;
  //  double getAlpha() const { return _alpha; };
  //  bool setAlpha(double alpha);

protected:
  double _evaluateCov(double h) const override;

private:
  //  double _alpha;
};

} // namespace gstlrn