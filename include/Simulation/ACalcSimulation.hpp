/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Calculators/ACalcInterpolator.hpp"

class GSTLEARN_EXPORT ACalcSimulation: public ACalcInterpolator
{
public:
  ACalcSimulation(int nbimu, int seed = 4324324);
  ACalcSimulation(const ACalcSimulation& r) = delete;
  ACalcSimulation& operator=(const ACalcSimulation& r) = delete;
  virtual ~ACalcSimulation();

  int getSeed() const { return _seed; }
  int getNbSimu() const { return _nbsimu; }
  void setSeed(int seed) { _seed = seed; }
  void setNbSimu(int nbsimu) { _nbsimu = nbsimu; }

private:
  int _nbsimu;
  int _seed;
};
