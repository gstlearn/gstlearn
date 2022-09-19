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
