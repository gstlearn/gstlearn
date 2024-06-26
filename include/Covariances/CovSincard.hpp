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

class GSTLEARN_EXPORT CovSincard : public ACovFunc
{
public:
  CovSincard(const CovContext& ctx);
  CovSincard(const CovSincard &r);
  CovSincard& operator= (const CovSincard &r);
  virtual ~CovSincard();

  double         getScadef() const override;
  virtual String getFormula() const override;
  String         getCovName() const override { return "Cardinal Sine"; }
  int            getMinOrder() const override { return -1; }
  bool           getCompatibleSpaceR() const override { return true; }

  bool isValidForTurningBand() const override { return true; }
  double simulateTurningBand(double t0, TurningBandOperate &operTB) const override;

protected:
  double _evaluateCov(double h) const override;
};

