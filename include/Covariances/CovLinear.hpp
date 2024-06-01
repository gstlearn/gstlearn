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

class GSTLEARN_EXPORT CovLinear : public ACovFunc
{
public:
  CovLinear(const CovContext& ctx);
  CovLinear(const CovLinear &r);
  CovLinear& operator= (const CovLinear &r);
  virtual ~CovLinear();

  int    hasRange() const override { return -1; }
  int    getMinOrder()  const override { return 0; }
  String getCovName() const override { return "Linear"; }

  bool isValidForTurningBand() const override { return true; }
  double simulateTurningBand(double t0,
                             const VectorDouble &t,
                             TurningBandOperate &operTB) const override;

protected:
  double _evaluateCov(double h) const override;
};

