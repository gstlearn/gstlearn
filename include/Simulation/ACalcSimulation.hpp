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

#include "Calculators/ACalcInterpolator.hpp"

namespace gstlrn
{
class GSTLEARN_EXPORT ACalcSimulation: public ACalcInterpolator
{
public:
  ACalcSimulation(Id nbsimu, Id seed = 4324324);
  ACalcSimulation(const ACalcSimulation& r)            = delete;
  ACalcSimulation& operator=(const ACalcSimulation& r) = delete;
  virtual ~ACalcSimulation();

  Id getSeed() const { return _seed; }
  Id getNbSimu() const { return _nbsimu; }
  void setSeed(Id seed) { _seed = seed; }
  void setNbSimu(Id nbsimu) { _nbsimu = nbsimu; }

protected:
  bool _check() override;
  bool _preprocess() override;

private:
  Id _nbsimu;
  Id _seed;
};

} // namespace gstlrn