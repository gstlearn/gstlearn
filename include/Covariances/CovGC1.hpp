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
#include "Basic/Law.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/ACovFunc.hpp"

/* Be careful ! This is not a real covariance */

class CovContext;
class TurningBandOperate;

class GSTLEARN_EXPORT CovGC1 : public ACovFunc
{
public:
  CovGC1(const CovContext& ctx);
  CovGC1(const CovGC1 &r);
  CovGC1& operator= (const CovGC1 &r);
  virtual ~CovGC1();

  int    hasRange()    const override { return -1; }
  int    getMinOrder() const override { return 0; }
  String getCovName()  const override { return "Order-1 G.C."; }
  bool   getCompatibleSpaceR() const override { return true; }

  bool isValidForTurningBand() const override { return true; }
  double simulateTurningBand(double t0,TurningBandOperate &operTB) const override;

protected:
  double _evaluateCov(double h) const override;
};

