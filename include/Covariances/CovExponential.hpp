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

class GSTLEARN_EXPORT CovExponential : public ACovFunc
{
public:
  CovExponential(const CovContext& ctx);
  CovExponential(const CovExponential &r);
  CovExponential& operator= (const CovExponential &r);
  virtual ~CovExponential();

  double getScadef() const override;
  virtual String getFormula() const override;
  String         getCovName() const override { return "Exponential"; }
  int            getMinOrder() const override { return -1; }

  bool isValidForTurningBand() const override { return true; }
  double simulateTurningBand(double t0, TurningBandOperate &operTB) const override;

  bool isValidForSpectral() const override { return true; }
  MatrixRectangular simulateSpectralOmega(int nb) const override;

protected:
  double _evaluateCov(double h) const override;
};

