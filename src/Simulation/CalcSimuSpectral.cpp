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
#include "Basic/NamingConvention.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovList.hpp"
#include "Db/Db.hpp"
#include "Enum/ESpaceType.hpp"
#include "Model/Model.hpp"
#include "Model/ModelGeneric.hpp"
#include "Simulation/SimuSpectralRN.hpp"
#include "Simulation/SimuSpectralS2.hpp"
#include "Stats/Classical.hpp"
#include "geoslib_define.h"

#include <cmath>

namespace gstlrn
{
CalcSimuSpectral::CalcSimuSpectral(Id nbsimu, Id ns, Id nd, Id seed, bool verbose, bool performedOnRN)
  : ACalcSimulation(nbsimu, seed)
  , _verbose(verbose)
  , _performedOnRN(performedOnRN)
  , _iattOut(-1)
  , _ns(ns)
  , _nd(nd)
  , _phi()
{
}

CalcSimuSpectral::~CalcSimuSpectral()
{
}

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
    messerr("The Space dimension of 'dbout'(%d) should match the one of Model(%d)",
            getDbout()->getNDim(), _getNDim());
    return 1;
  }

  // Check that the Model is compatible with Spectral Simulation
  if (!isValidForSpectral()) return false;

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
bool CalcSimuSpectral::isValidForSpectral() const
{
  auto ncova            = _getNCov();
  const auto* modellist = dynamic_cast<const ModelCovList*>(getModelGeneric());

  // Loop on the simulations
  for (Id is = 0, ns = MAX(ncova, 1); is < ns; is++)
  {
    if (ncova <= 0)
    {
      const auto* cova = getModelGeneric()->getCov();
      if (_performedOnRN)
      {
        if (!cova->isValidForSpectralOnRn())
        {
          messerr("The covariance component %s of the Model is not valid for Spectral Simulation on Rn", is + 1);
          return false;
        }
      }
      else
      {
        if (!cova->isValidForSpectralOnSphere())
        {
          messerr("The covariance component %d of the Model is not valid for Spectral Simulation on the Sphere", is + 1);
          return false;
        }
      }
    }
    else
    {
      const auto* covbase = dynamic_cast<const CovAniso*>(modellist->getCovBase(is));
      if (_performedOnRN)
      {
        if (!covbase->isValidForSpectralOnRn())
        {
          messerr("The covariance component %d of the Model is not valid for Spectral Simulation on Rn", is + 1);
          return false;
        }
      }
      else
      {
        if (!covbase->isValidForSpectralOnSphere())
        {
          messerr("The covariance component %d of the Model is not valid for Spectral Simulation on the Sphere", is + 1);
          return false;
        }
      }
    }
  }
  return true;
}

/**
 * Simulate the spectrum components for Rn or S2 for one simulation
 *
 * @param cova the covariance function to be simulated
 */
Id CalcSimuSpectral::simulate(const ACov* cova)
{
  // simulation of random phases, uniform on [0,2 pi]
  auto ns    = _getNs();
  _phi       = VectorDouble(ns);
  double pi2 = 2. * GV_PI;
  for (Id is = 0; is < ns; is++)
    _phi[is] = pi2 * law_uniform();

  return _simulate(cova);
}

/**
 * @brief Compute one non-conditional simulation on the samples of Dbout using Spectral Method
 */
Id CalcSimuSpectral::compute(Db* dbout, Id isimu)
{
  auto nvar   = _getNVar();
  auto nbsimu = getNbSimu();
  auto nech   = dbout->getNSample();

  VectorVectorDouble tab(nvar);
  for (Id ivar = 0; ivar < nvar; ivar++) tab[ivar].resize(nech);
  VectorBool activeArray = dbout->getActiveArray();

  // The next line has been added to allow using method 'compute' independently
  // In the calculator framework, it is simply bypassed
  if (getDbout() == nullptr) setDbout(dbout);
  if (_iattOut < 0)
    _iattOut = dbout->addColumnsByConstant(nvar * nbsimu, 0., "Simu", ELoc::SIMU);

  // Compute one simulation
  if (_compute(dbout, activeArray, tab)) return 1;

  // Save the resulting array
  saveResults(dbout, 0, 0, isimu, activeArray, tab);

  return 0;
}

bool CalcSimuSpectral::_preprocess()
{
  if (!ACalcSimulation::_preprocess()) return false;

  auto nvar   = _getNVar();
  auto nbsimu = getNbSimu();

  // Add the attributes for storing the results
  if (getDbin() != nullptr)
  {
    Id iptr_in = _addVariableDb(1, 2, ELoc::SIMU, 0, nvar * nbsimu);
    if (iptr_in < 0) return false;
  }

  // Factorize the matrix of sills
  auto* modelLocal = dynamic_cast<Model*>(getModelGeneric());
  if (modelLocal != nullptr)
    modelLocal->computeAic();

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
  Id nvar   = _getNVar();
  for (Id isimu = 0; isimu < nbsimu; isimu++)
    for (Id ivar = 0; ivar < nvar; ivar++)
    {
      String ps(prefix);
      ps.append(delim + "V" + std::to_string(ivar + 1));
      ps.append(delim + "S" + std::to_string(isimu + 1));
      namconvS.setPrefix(ps);
      namconvS.setNamesAndLocators(nullptr, VectorString(), ELoc::Z, 1, getDbout(),
                                   _iattOut + isimu * nvar + ivar, "", 1);
    }
  return true;
}

bool CalcSimuSpectral::_run()
{
  auto nbsimu = getNbSimu();
  auto nvar   = _getNVar();
  auto* db    = getDbout();
  auto nech   = db->getNSample();

  // Set the random seed
  Id mem_seed = law_get_random_seed();
  law_set_random_seed(getSeed());

  auto ncova = _getNCov();
  const ACov* cova;
  const CovAniso* covbase;
  const auto* modellist = dynamic_cast<const ModelCovList*>(getModelGeneric());

  VectorVectorDouble tab(nvar);
  for (Id ivar = 0; ivar < nvar; ivar++) tab[ivar].resize(nech);
  VectorBool activeArray = db->getActiveArray();

  // Loop on the simulations
  for (Id isimu = 0; isimu < nbsimu; isimu++)
  {
    for (Id is = 0, ns = MAX(ncova, 1); is < ns; is++)
    {
      if (getVerbose())
        messerr(">>> computing simulation %d for covariance %d", isimu + 1, is + 1);
      if (ncova <= 0)
        cova = getModelGeneric()->getCov();
      else
      {
        covbase = dynamic_cast<const CovAniso*>(modellist->getCovBase(is));
        if (covbase->getType() == ECov::NUGGET) continue;
        cova = dynamic_cast<const ACov*>(covbase);
      }

      // Blank out the array 'tab'
      for (Id ivar = 0; ivar < nvar; ivar++)
        tab[ivar].fill(0.);

      if (simulate(cova)) return false;

      // Compute one simulation
      if (_compute(db, activeArray, tab)) return 1;

      // Save the resulting array
      if (ncova <= 0)
        saveResults(db, 0, 0, isimu, activeArray, tab);
      else
        scaleAndSaveResults(db, covbase, 0, 0, isimu, activeArray, tab);
    }
  }

  // Set the initial seed back
  law_set_random_seed(mem_seed);
  return true;
}

/**
 * Perform a series of simulations (on Rn or on the Sphere) using Spectral Method
 *
 * @param dbin Input Db where the conditioning data are read
 * @param dbout Output Db where the results are stored
 * @param model ModelGeneric structure
 * @param neigh Neighborhood structure
 * @param nbsimu Number of simulations processed simultaneously
 * @param seed Seed used for the Random number generator
 * @param ns Number of spectral components
 * @param nd Maximum number of spectral orders on S2
 * @param cov0 Auxiliary covariance used for importance sampling
 * @param verbose Verbose flag
 * @param namconv Naming Convention
 *
 * @note The conditional version is not yet available
 */
Id simuSpectral(Db* dbin,
                Db* dbout,
                ModelGeneric* model,
                ANeigh* neigh,
                Id nbsimu,
                Id seed,
                Id ns,
                Id nd,
                const ACov* cov0,
                bool verbose,
                const NamingConvention& namconv)
{
  // Check the space type
  bool isSimuRN;
  const auto space = model->getContext()->getSpace();
  if (space->getType() == ESpaceType::COMPOSITE)
  {
    // The RN x time model is simulated as a R(N+1) model (see CorGneiting)
    isSimuRN = (space->getComponent(0)->getType() == ESpaceType::RN);
  }
  else
  {
    isSimuRN = (space->getType() == ESpaceType::RN);
  }

  // Instantiate the Calculator
  std::unique_ptr<CalcSimuSpectral> spectral;

  // Instantiate
  if (isSimuRN)
  {
    spectral = std::make_unique<SimuSpectralRN>(nbsimu, ns, nd, seed, cov0, verbose);
  }
  else
  {
    spectral = std::make_unique<SimuSpectralS2>(nbsimu, ns, nd, seed, verbose);
  }

  // Set the members of the Calculator
  spectral->setDbin(dbin);
  spectral->setDbout(dbout);
  spectral->setModelGeneric(model);
  spectral->setNeigh(neigh);
  spectral->setNamingConvention(namconv);

  // Run the calculator
  Id error = (spectral->run()) ? 0 : 1;
  return error;
}
} // namespace gstlrn