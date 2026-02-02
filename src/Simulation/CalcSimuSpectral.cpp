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
#include "Db/Db.hpp"
#include "Enum/ESpaceType.hpp"
#include "Model/Model.hpp"
#include "Model/ModelGeneric.hpp"
#include "Simulation/SimuSpectralRN.hpp"
#include "Simulation/SimuSpectralS2.hpp"
#include "Space/ASpaceObject.hpp"
#include "Stats/Classical.hpp"
#include "geoslib_define.h"

#include <cmath>

namespace gstlrn
{
CalcSimuSpectral::CalcSimuSpectral(Id nbsimu, Id ns, Id nd, Id seed, bool verbose)
  : ACalcSimulation(nbsimu, seed)
  , _isPrepared(false)
  , _verbose(verbose)
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

  // Check that the Model is compatible with Spectral Simulation
  if (!isValidForSpectral(getModelGeneric())) return false;

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
 ** \param[in]  model    Model structure
 **
 *****************************************************************************/
bool CalcSimuSpectral::isValidForSpectral(const ModelGeneric* model)
{
  /* Loop on the structures */

  const auto* cova = model->getCov();
  if (cova == nullptr) return false;
  return cova->isValidForSpectral();
}

/**
 * Simulate the spectrum components for Rn or S2 for one simulation
 *
 * @param verbose Verbose flag
 * @param cov0 the auxiliary covariance function used for importance sampling
 */
Id CalcSimuSpectral::simulate()
{
  // simulation of random phases, uniform on [0,2 pi]
  auto ns    = _getNs();
  _phi       = VectorDouble(ns);
  double pi2 = 2. * GV_PI;
  for (Id is = 0; is < ns; is++)
    _phi[is] = pi2 * law_uniform();

  // simulation of the random frequencies
  if (_simulate() != 0) return 1;
  _setIsPrepared(true);
  return 0;
}

/**
 * @brief Compute one non-conditional simulation on the samples of Dbout using Spectral Method
 *
 * @param dbout
 * @param isimu
 * @return Id
 */
Id CalcSimuSpectral::compute(Db* dbout, Id isimu)
{
  if (!_getIsPrepared())
  {
    messerr("You should run 'simulate' beforehand");
    return 1;
  }
  Id ndim = dbout->getNDim();
  Id nech = dbout->getNSample(true);

  if (ndim != _getNDim())
  {
    messerr("The Space dimension of 'dbout'(%d) should match the one of Model(%d)", ndim);
    return 1;
  }
  if (nech <= 0)
  {
    messerr("'dbout' must have a positive number of active samples");
    return 1;
  }

  // computing the values
  if (_compute(dbout, isimu) > 0) return 1;

  return 0;
}

bool CalcSimuSpectral::_preprocess()
{
  if (!ACalcSimulation::_preprocess()) return false;

  auto nvar   = _getNVar();
  auto nbsimu = getNbSimu();

  /* Add the attributes for storing the results */

  if (getDbin() != nullptr)
  {
    Id iptr_in = _addVariableDb(1, 2, ELoc::SIMU, 0, nvar * nbsimu);
    if (iptr_in < 0) return false;
  }

  _iattOut = _addVariableDb(2, 1, ELoc::SIMU, 0, nvar * nbsimu);
  return _iattOut >= 0;
}

bool CalcSimuSpectral::_postprocess()
{
  // Free the temporary variables
  _cleanVariableDb(2);

  // Clean variables created for Expansion
  if (_expandInformation(-1, ELoc::F)) return false;
  if (_expandInformation(-1, ELoc::NOSTAT)) return false;

  // Set the error return flag

  // _renameVariable(2, VectorString(), ELoc::Z, _getNVar(), _iattOut, String(), getNbSimu());

  NamingConvention namconv = getNamingConvention();
  String prefix(namconv.getPrefix());
  String delim(namconv.getDelim());
  NamingConvention namconvS(namconv);

  // Loop on the simulations
  for (Id isimu = 0; isimu < getNbSimu(); isimu++)
    for (Id ivar = 0; ivar < _getNVar(); ivar++)
    {
      String ps(prefix);
      ps.append(delim + "V" + std::to_string(ivar + 1));
      ps.append(delim + "S" + std::to_string(isimu + 1));
      namconvS.setPrefix(ps);
      namconvS.setNamesAndLocators(nullptr, VectorString(), ELoc::Z, 1, getDbout(),
                                   _iattOut + isimu * _getNVar() + ivar, "", 1);
    }
  return true;
}

bool CalcSimuSpectral::_run()
{
  law_set_random_seed(getSeed());
  auto nbsimu = getNbSimu();

  // Loop on the simulations
  for (Id isimu = 0; isimu < nbsimu; isimu++)
  {
    if (getVerbose())
      messerr(">>> computing simulation %d", isimu + 1);
    if (simulate()) return false;
    if (compute(getDbout(), isimu)) return false;
  }
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
  if (getDefaultSpaceType() == ESpaceType::COMPOSITE)
  {
    // The RN x time model is simulated as a R(N+1) model (see CorGneiting)
    isSimuRN = (getDefaultSpace()->getComponent(0)->getType() == ESpaceType::RN);
  }
  else
  {
    isSimuRN = (getDefaultSpaceType() == ESpaceType::RN);
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