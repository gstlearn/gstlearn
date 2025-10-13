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
#include "Simulation/ACalcSimulation.hpp"
#include "Calculators/ACalcInterpolator.hpp"

namespace gstlrn
{

ACalcSimulation::ACalcSimulation(Id nbsimu, Id seed)
  : ACalcInterpolator()
  , _nbsimu(nbsimu)
  , _seed(seed)
{
}

ACalcSimulation::~ACalcSimulation()
{
}

bool ACalcSimulation::_check()
{
  if (! ACalcInterpolator::_check()) return false;

  if (getNbSimu() <= 0)
  {
    messerr("You must define 'nbsimu' and 'nbtuba'");
    return false;
  }
  return true;
}

bool ACalcSimulation::_preprocess()
{
  return ACalcInterpolator::_preprocess();
}
}