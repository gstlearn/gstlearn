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

class CovContext;
class TurningBandOperate;

class GSTLEARN_EXPORT CovGCspline : public ACovFunc
{
public:
  CovGCspline(const CovContext& ctx);
  CovGCspline(const CovGCspline &r);
  CovGCspline& operator= (const CovGCspline &r);
  virtual ~CovGCspline();

  int            hasRange() const override { return -1; }
  String         getCovName() const override { return "Spline G.C."; }
  int            getMinOrder() const override { return 1; }
  bool           getCompatibleSpaceR() const override { return true; }

  bool isValidForTurningBand() const override { return true; }
  double simulateTurningBand(double t0, TurningBandOperate &operTB) const override;

protected:
  double _evaluateCov(double h) const override;
};

