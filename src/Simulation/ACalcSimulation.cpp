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
#include "Basic/Law.hpp"
#include "Calculators/ACalcInterpolator.hpp"
#include "Estimation/KrigingSystem.hpp"

namespace gstlrn
{

  ACalcSimulation::ACalcSimulation(Id nbsimu, Id seed, bool verbose)
    : ACalcInterpolator(verbose)
    , _nbsimu(nbsimu)
    , _seed(seed)
    , _shift(0)
    , _seedPerSimu()
  {
  }

  ACalcSimulation::~ACalcSimulation() {}

  bool ACalcSimulation::_check()
  {
    if (!ACalcInterpolator::_check()) return false;

    if (getNbSimu() <= 0)
    {
      messerr("You must define 'nbsimu' and 'nbtuba'");
      return false;
    }

    // Memorize the current seed
    Id mem_seed = law_get_random_seed();
    law_set_random_seed(getSeed());

    _seedPerSimu.resize(getNbSimu());
    for (Id isimu = 0; isimu < getNbSimu(); isimu++)
    {
      // Next statement is simply meant to move the seed for random number generator
      law_uniform();
      _seedPerSimu[isimu] = law_get_random_seed();
    }

    // Restore the seed
    law_set_random_seed(mem_seed);
    return true;
  }

  Id ACalcSimulation::getNVar() const
  {
    if (getModelGeneric() == nullptr) return 0;
    return getModelGeneric()->getNVar();
  }

} // namespace gstlrn
