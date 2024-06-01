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
#include "Covariances/ACovFunc.hpp"

class CovContext;
class TurningBandOperate;

class GSTLEARN_EXPORT CovGaussian : public ACovFunc
{
public:
  CovGaussian(const CovContext& ctx);
  CovGaussian(const CovGaussian &r);
  CovGaussian& operator= (const CovGaussian &r);
  virtual ~CovGaussian();

  virtual String getFormula() const override;
  String         getCovName() const override { return "Gaussian"; }
  int            getMinOrder() const override { return -1; }
  double         getScadef() const override;
  virtual bool   hasCovDerivative() const override { return true; }

  bool isValidForTurningBand() const override { return true; }
  double simulateTurningBand(double t0,
                             const VectorDouble &t,
                             TurningBandOperate &operTB) const override;

protected:
  double _evaluateCov(double h)  const override;
  double _evaluateCovDerivative(int degree, double h) const override;
};

