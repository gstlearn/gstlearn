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

#include "Covariances/ACovFunc.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
// Forward declaration
class CovContext;
class TurningBandOperate;

class GSTLEARN_EXPORT CovCubic: public ACovFunc
{
public:
  CovCubic(const CovContext& ctx);
  CovCubic(const CovCubic& r);
  CovCubic& operator=(const CovCubic& r);
  virtual ~CovCubic();

  size_t getMaxNDim() const override { return 3; }

  String getFormula() const override;
  String getCovName() const override { return "Cubic"; }
  Id getMinOrder() const override { return -1; }
  bool getCompatibleSpaceR() const override { return true; }
  bool hasCovDerivative() const override { return true; }

  bool isValidForTurningBand() const override { return true; }
  double simulateTurningBand(double t0, TurningBandOperate& operTB) const override;

protected:
  double _evaluateCov(double h) const override;
  double _evaluateCovDerivative(Id degree, double h) const override;
  double _evaluateCovFirstDerivativeOverH(double h) const override;
};

} // namespace gstlrn