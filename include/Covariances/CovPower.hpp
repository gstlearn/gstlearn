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

/* Be careful ! This is not a real covariance */

namespace gstlrn
{
class CovContext;
class TurningBandOperate;

class GSTLEARN_EXPORT CovPower : public ACovFunc
{
public:
  CovPower(const CovContext& ctx);
  CovPower(const CovPower &r);
  CovPower& operator= (const CovPower &r);
  virtual ~CovPower();

  Id            hasRange()    const override { return -1; }
  bool           hasParam()    const override { return true; }
  double         getParMax()   const override { return 1.99; }
  Id            getMinOrder() const override { return 0; }
  String         getCovName()  const override { return "Power"; }
  bool           getCompatibleSpaceR() const override { return true; }

  bool isValidForTurningBand() const override { return true; }
  double simulateTurningBand(double t0, TurningBandOperate &operTB) const override;

protected:
  double _evaluateCov(double h) const override;
};

}