/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Simulation/ACalcSimulation.hpp"
#include "Calculators/ACalcInterpolator.hpp"

ACalcSimulation::ACalcSimulation(int nbsimu, int seed)
    : ACalcInterpolator(),
      _nbsimu(nbsimu),
      _seed(seed)
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
  if (!ACalcInterpolator::_preprocess()) return false;

  return true;
}
