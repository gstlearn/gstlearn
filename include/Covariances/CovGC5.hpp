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
#include "Basic/Law.hpp"
#include "Covariances/CovContext.hpp"

/* Be careful ! This is not a real covariance */

class CovContext;
class TurningBandOperate;

class GSTLEARN_EXPORT CovGC5 : public ACovFunc
{
public:
  CovGC5(const CovContext& ctx);
  CovGC5(const CovGC5 &r);
  CovGC5& operator= (const CovGC5 &r);
  virtual ~CovGC5();

  int    hasRange() const override { return -1; }
  int    getMinOrder()  const override { return 2; }
  String getCovName() const override { return "Order-5 G.C."; }

  bool isValidForTurningBand() const override { return true; }
  double simulateTurningBand(double t0,
                             const VectorDouble &t,
                             TurningBandOperate &operTB) const override;

protected:
  double _evaluateCov(double h) const override;
};

