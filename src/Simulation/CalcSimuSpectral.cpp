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
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "Enum/ESimuType.hpp"
#include "Model/Model.hpp"
#include "Model/ModelGeneric.hpp"
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
    : ACalcSimulation(nbsimu, seed, verbose)
    , _iattOut(-1)
    , _ns(ns)
    , _nd(nd)
  {
  }

  CalcSimuSpectral::~CalcSimuSpectral() {}

  Id CalcSimuSpectral::_getNDim() const
  {
    if (getModelGeneric() == nullptr) return 0;
    return getModelGeneric()->getNDim();
  }

  Id CalcSimuSpectral::_getNVar() const
  {
    if (getModelGeneric() == nullptr) return 0;
    return getModelGeneric()->getNVar();
  }

  bool CalcSimuSpectral::_check()
  {
    if (!ACalcSimulation::_check()) return false;

    if (!hasDbout()) return false;
    if (!hasModelGeneric()) return false;
    if (hasDbin(false))
    {
      if (!hasNeigh()) return false;
    }
    if (getDbout()->getNDim() != _getNDim())
    {
      messerr(
        "The Space dimension of 'dbout'(%d) should match the one of Model(%d)",
        getDbout()->getNDim(), _getNDim());
      return 1;
    }

    // Check that the Model is compatible with Spectral Simulation
    if (!isValidForSimulation(ESimuType::SPECTRAL)) return false;

    if (_getNs() <= 0)
    {
      messerr("The number of simulated harmonic components should be positive");
      return false;
    }
    return true;
  }

  /****************************************************************************/
  /*!
   **  Check if the Model can be simulated using Spectral Method
   **
   ** \return  True if the Model is valid; 0 otherwise
   **
   *****************************************************************************/
  bool CalcSimuSpectral::isValidForSimulation(const ESimuType& simuType) const
  {
    auto ncova = _getNCov();
    const auto* modellist =
      dynamic_cast<const ModelCovList*>(getModelGeneric());

    // Loop on the simulations
    for (Id is = 0, ns = MAX(ncova, 1); is < ns; is++)
    {
      if (ncova <= 0)
      {
        const auto* cova = getModelGeneric()->getCov();
        if (!cova->isValidForSimulation(simuType))
        {
          messerr(
            "The covariance component %d of the Model is not valid for "
            "%s simulation",
            is + 1, simuType.getKey());
          return false;
        }
      }
      else
      {
        const auto* covbase =
          dynamic_cast<const CovAniso*>(modellist->getCovBase(is));
        if (!covbase->isValidForSimulation(simuType))
        {
          messerr(
            "The covariance component %d of the Model is not valid for "
            "%s Simulation",
            is + 1, simuType.getKey());
          return false;
        }
      }
    }
    return true;
  }

  /**
   * Simulate the spectrum components for Rn or S2 for one simulation
   */
  Id CalcSimuSpectral::simulateTest()
  {
    // When called outside the calculator framework, we need to set the random seed here
    law_set_random_seed(getSeed());

    return _simulate();
  }

  /**
   * @brief Compute one non-conditional simulation on the samples of Dbout using Spectral Method
   */
  Id CalcSimuSpectral::computeTest(Db* dbout, Id isimu)
  {
    VectorVectorDouble tab;
    VectorBool activeArray;
    _allocateForOneSimulation(dbout, _getNVar(), activeArray, tab);

    // The next line has been added to allow using method 'compute' independently
    // In the calculator framework, it is simply bypassed
    if (getDbout() == nullptr) setDbout(dbout);
    if (_iattOut < 0)
      _iattOut = dbout->addColumnsByConstant(
        _getNVar() * getNbSimu(), 0., "Simu", ELoc::SIMU);

    // Compute one simulation
    if (_compute(dbout, isimu, activeArray, tab)) return 1;

    // Save the resulting array
    saveResults(dbout, isimu, activeArray, tab);

    return 0;
  }

  bool CalcSimuSpectral::_preprocess()
  {
    if (!ACalcSimulation::_preprocess()) return false;

    auto nvar = _getNVar();
    auto nbsimu = getNbSimu();

    // Add the attributes for storing the results
    if (getDbin() != nullptr)
    {
      Id iptr_in = _addVariableDb(1, 2, ELoc::SIMU, 0, nvar * nbsimu);
      if (iptr_in < 0) return false;
    }

    // Factorize the matrix of sills
    auto* modelLocal = dynamic_cast<Model*>(getModelGeneric());
    if (modelLocal != nullptr) modelLocal->computeAic();

    _iattOut = _addVariableDb(2, 1, ELoc::SIMU, 0, nvar * nbsimu);
    return _iattOut >= 0;
  }

  bool CalcSimuSpectral::_postprocess()
  {
    // Free the temporary variables
    _cleanVariableDb(2);

    // _renameVariable(2, VectorString(), ELoc::Z, _getNVar(), _iattOut, String(), getNbSimu());

    NamingConvention namconv = getNamingConvention();
    String prefix(namconv.getPrefix());
    String delim(namconv.getDelim());
    NamingConvention namconvS(namconv);

    // Loop on the simulations
    Id nbsimu = getNbSimu();
    Id nvar = _getNVar();
    for (Id isimu = 0; isimu < nbsimu; isimu++)
      for (Id ivar = 0; ivar < nvar; ivar++)
      {
        String ps(prefix);
        ps.append(delim + "V" + std::to_string(ivar + 1));
        ps.append(delim + "S" + std::to_string(isimu + 1));
        namconvS.setPrefix(ps);
        namconvS.setNamesAndLocators(
          nullptr, VectorString(), ELoc::Z, 1, getDbout(),
          //                                   _iattOut + isimu * nvar + ivar, "", 1);
          _iattOut + ivar * nbsimu + isimu, "", 1);
      }
    return true;
  }

  bool CalcSimuSpectral::_run()
  {
    auto* db = getDbout();

    law_set_random_seed(getSeed());

    VectorVectorDouble tab;
    VectorBool activeArray;
    _allocateForOneSimulation(db, _getNVar(), activeArray, tab);

    // Loop on the simulations
    for (Id isimu = 0, nbsimu = getNbSimu(); isimu < nbsimu; isimu++)
    {
      if (getVerbose()) message(">>> computing simulation %d\n", isimu + 1);

      // simulate the spectrum
      _simulate();

      // Compute one simulation
      _compute(db, isimu, activeArray, tab);

      // Save the resulting array
      saveResults(db, isimu, activeArray, tab);
    }

    return true;
  }

} // namespace gstlrn
