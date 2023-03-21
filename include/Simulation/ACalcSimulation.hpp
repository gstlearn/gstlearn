/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
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

protected:
  virtual bool _check() override;
  virtual bool _preprocess() override;

private:
  int _nbsimu;
  int _seed;
};
