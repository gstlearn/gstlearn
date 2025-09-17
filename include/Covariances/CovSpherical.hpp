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

namespace gstlrn
{
class CovContext;
class TurningBandOperate;

class GSTLEARN_EXPORT CovSpherical : public ACovFunc
{
public:
  CovSpherical(const CovContext& ctx);
  CovSpherical(const CovSpherical &r);
  CovSpherical& operator= (const CovSpherical &r);
  virtual ~CovSpherical();

  size_t getMaxNDim() const override { return 3; }
  String         getFormula() const override;
  String         getCovName() const override { return "Spherical"; }
  Id            getMinOrder() const override { return -1; }
  bool           getCompatibleSpaceR() const override { return true; }
  bool isValidForTurningBand() const override { return true; }
  double simulateTurningBand(double t0, TurningBandOperate &operTB) const override;

protected:
  double _evaluateCov(double h) const override;
  double _evaluateCovFirstDerivative(double h) const override;
};

}