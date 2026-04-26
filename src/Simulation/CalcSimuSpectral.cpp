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
#include "Simulation/CalcSimuSpectral.hpp"
#include "Basic/Law.hpp"
#include "Basic/Message.hpp"
#include "Basic/NamingConvention.hpp"
#include "Covariances/ACov.hpp"
#include "Db/Db.hpp"
#include "Enum/ESimuType.hpp"
#include "Model/Model.hpp"
#include "Model/ModelGeneric.hpp"
#include "Simulation/ACalcSimuGaussian.hpp"
#include "Simulation/SimuSpectralRN.hpp"
#include "Simulation/SimuSpectralS2.hpp"
#include "Stats/Classical.hpp"
#include "geoslib_define.h"

#include <cmath>

namespace gstlrn
{
  CalcSimuSpectral::CalcSimuSpectral(
    Id nbsimu,
    Id ns,
    Id nd,
    Id seed,
    bool verbose)
    : ACalcSimuGaussian(nbsimu, seed, verbose)
    , _ns(ns)
    , _nd(nd)
  {
  }

  CalcSimuSpectral::~CalcSimuSpectral() {}

  bool CalcSimuSpectral::_check()
  {
    if (!ACalcSimuGaussian::_check()) return false;

    if (!hasModelGeneric()) return false;

    // Check that the Model is compatible with Spectral Simulation
    if (!getModelGeneric()->isValidForSimulation(ESimuType::SPECTRAL))
      return false;

    if (_getNs() <= 0)
    {
      messerr("The number of simulated harmonic components should be positive");
      return false;
    }
    return true;
  }

  /**
   * Simulate the spectrum components for Rn or S2 for one simulation
   */
  Id CalcSimuSpectral::simulateSpectral()
  {
    // When called outside the calculator framework, we need to set the random seed here
    law_set_random_seed(getSeed());
    if (!_check()) return 1;

    return _simulate(0);
  }

  /**
   * @brief Compute one non-conditional simulation on the samples of Dbout using Spectral Method
   */
  Id CalcSimuSpectral::computeSpectral(Db* dbout, Id isimu)
  {
    VectorVectorDouble tab;
    VectorBool activeArray;
    _allocateForOneSimulation(dbout, _getNVar(), activeArray, tab);

    // The next line has been added to allow using method 'computeSpectral' independently.
    if (getDbout() == nullptr) setDbout(dbout);

    // Initialize the output variable
    _initializeOutputAttribute();

    // Compute one simulation
    if (_compute(dbout, activeArray, tab)) return 1;

    // Save the resulting array
    saveResults(dbout, isimu, activeArray, tab);

    return 0;
  }

  bool CalcSimuSpectral::_preprocess()
  {
    if (!ACalcSimuGaussian::_preprocess()) return false;

    // Factorize the matrix of sills
    auto* modelLocal = dynamic_cast<Model*>(getModelGeneric());
    if (modelLocal != nullptr) modelLocal->computeAic();
    return true;
  }

} // namespace gstlrn
