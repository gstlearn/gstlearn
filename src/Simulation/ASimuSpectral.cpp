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
#include "Simulation/ASimuSpectral.hpp"
#include "Basic/Law.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/VectorNumT.hpp"
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
#include <string>

namespace gstlrn
{
ASimuSpectral::ASimuSpectral(const ACov* cova)
  : _isPrepared(false)
  , _phi()
  , _cova(cova)
{
}

ASimuSpectral::ASimuSpectral(const ASimuSpectral& r)
  : _isPrepared(r._isPrepared)
  , _phi(r._phi)
  , _cova(r._cova)
{
}

ASimuSpectral& ASimuSpectral::operator=(const ASimuSpectral& r)
{
  if (this != &r)
  {
    _isPrepared = r._isPrepared;
    _phi        = r._phi;
    _cova       = r._cova;
  }
  return *this;
}

ASimuSpectral::~ASimuSpectral()
{
}

Id ASimuSpectral::getNs() const
{
  return _phi.length();
}
Id ASimuSpectral::getNDim() const
{
  return _cova->getNDim();
}
Id ASimuSpectral::getNVar() const
{
  return _cova->getNVar();
}

/**
 * Simulate the spectrum components for Rn or S2
 *
 * @param ns Number of components
 * @param seed Seed for random number generation: avoid setting the seed)
 * @param verbose Verbose flag
 * @param cov0 the auxiliary covariance function used for importance sampling
 * @param nd Maximum order of the spectrum on S2
 */
Id ASimuSpectral::simulate(Id ns, Id seed, bool verbose, const ACov* cov0, Id nd)
{
  if (_cova == nullptr)
  {
    messerr("A Covariance should be attached beforehand");
    return 1;
  }
  if (!_cova->isValidForSpectral())
  {
    messerr("Simulation of the harmonic components is not implemented for this covariance");
    return 1;
  }
  Id ndim = _cova->getNDim();
  Id nvar = _cova->getNVar();
  if (ns <= 0)
  {
    messerr("The number of simulated harmonic components should be positive");
    return 1;
  }
  if (cov0 != nullptr)
  {
    if (!cov0->isValidForSpectral())
    {
      messerr("Simulation of the harmonic components is not implemented for the auxiliary covariance");
      return 1;
    }
    if (cov0->getNVar() > 1)
    {
      messerr("The auxiliary covariance should be scalar");
      return 1;
    }
  }

  // Optional printout
  if (verbose)
  {
    message("Simulation of the spectrum\n");
    message("- Space dimension   = %d\n", ndim);
    message("- Number of variables  = %d\n", nvar);
    message("- Number of spectral components = %d\n", ns);
    if (cov0 != nullptr)
    {
      message("Simulation using importance sampling\n");
    }
  }

  // perform the simulation
  if (seed > 0) law_set_random_seed(seed);

  // simulation of random phases, uniform on [0,2 pi]
  _phi       = VectorDouble(ns);
  double pi2 = 2. * GV_PI;
  for (Id is = 0; is < ns; is++)
    _phi[is] = pi2 * law_uniform();

  // simulation of the random frequencies
  if (_simulate(ns, nd, cov0, verbose) != 0) return 1;
  _isPrepared = true;
  return 0;
}

Id ASimuSpectral::compute(Db* dbout, Id iuid, bool verbose, const NamingConvention& namconv, const String& qualifier)
{
  Id ndim              = dbout->getNDim();
  Id nech              = dbout->getNSample(true);
  Id nvar              = _cova->getNVar();
  bool flagNewVariable = (iuid <= 0);

  if (ndim != _cova->getNDim())
  {
    messerr("The Space dimension of 'dbout'(%d) should match the one of Model(%d)", ndim);
    return 1;
  }
  if (nech <= 0)
  {
    messerr("'dbout' must have a positive number of active samples");
    return 1;
  }
  if (!_isPrepared)
  {
    messerr("You should run 'simulate' beforehand");
    return 1;
  }

  // Optional printout
  if (verbose)
  {
    message("Spectral Simulation on a set of Isolated Points\n");
    message("- Number of samples = %d\n", nech);
    message("- Space dimension   = %d\n", ndim);
    message("- Number of variables = %d\n", nvar);
  }

  // Create the variables
  if (flagNewVariable)
  {
    iuid = dbout->addColumnsByConstant(nvar, 0., String(), ELoc::Z);
    if (iuid < 0) return 1;
  }

  // computing the values
  if (_compute(dbout, iuid, verbose) > 0) return 1;

  // Modify the name of the output variables
  if (flagNewVariable)
  {
    for (Id ivar = 0; ivar < nvar; ivar++)
    {
      NamingConvention nmSim(namconv);
      String nm(namconv.getPrefix() + namconv.getDelim() + "V" + std::to_string(ivar + 1));
      nmSim.setPrefix(nm);
      nmSim.setNamesAndLocators(dbout, iuid + ivar, qualifier, 1);
    }
  }
  return 0;
}

/**
 * Perform a series of simulations (on Rn or on the Sphere) using Spectral Method
 *
 * @param dbin Input Db where the conditioning data are read
 * @param dbout Output Db where the results are stored
 * @param cova a covariance that can cope with spectral method
 * @param ns Number of spectral components
 * @param nd Maximum number of spectral orders on S2
 * @param nbsimu Number of simulations processed simultaneously
 * @param seed Seed used for the Random number generator
 * @param cov0 Auxiliary covariance used for importance sampling
 * @param verbose Verbose flag
 * @param namconv Naming Convention
 *

   if (getDefaultSpaceType() == ESpaceType::RN)
    SimuSpectralRN is used
  else
    SimuSpectralS2 is used

 * @note The conditional version is not yet available
 */
Id simuSpectral(Db* dbin,
                Db* dbout,
                const ACov* cova,
                Id nbsimu,
                Id seed,
                Id ns,
                Id nd,
                const ACov* cov0,
                bool verbose,
                const NamingConvention& namconv)
{
  if (dbin != nullptr)
  {
    messerr("The current version does not allow Conditional Simulations");
    return 1;
  }
  if (nbsimu <= 0)
  {
    messerr("You must provide a positive number of simulations");
    return 1;
  }
  if (dbout->getNDim() != cova->getNDim())
  {
    messerr("The Space dimension of 'dbout'(%d) should match the one of Model(%d)",
            dbout->getNDim(), cova->getNDim());
    return 1;
  }

  // Set the seed for Random number generator (for all simulations)
  if (seed > 0)
  {
    law_set_random_seed(seed);
    seed = 0;
  }

  // Creating the output variables
  Id nvar = cova->getNVar();
  Id iuid = dbout->addColumnsByConstant(nbsimu * nvar, 0., String(), ELoc::Z);
  if (iuid < 0) return 1;

  // Instantiate the SimuSpectral class
  std::unique_ptr<ASimuSpectral> simu = nullptr;
  if (getDefaultSpaceType() == ESpaceType::COMPOSITE)
  {
    if (getDefaultSpace()->getComponent(0)->getType() == ESpaceType::RN)
    { // The RN x time model is simulated as a R(N+1) model (see CorGneiting)
      simu = std::make_unique<SimuSpectralRN>(cova);
    }
    else
    {
      messerr("Space time model on S2 not yet implemented");
      // simu = std::make_unique<SimuSpectralS2>(cova);
    }
  }
  else if (getDefaultSpaceType() == ESpaceType::RN)
  {
    simu = std::make_unique<SimuSpectralRN>(cova);
  }
  else
  {
    simu = std::make_unique<SimuSpectralS2>(cova);
  }

  String prefix(namconv.getPrefix());
  String delim(namconv.getDelim());
  NamingConvention namconvS(namconv);

  // Loop on the simulations
  for (Id isimu = 0; isimu < nbsimu; isimu++)
  {
    if (verbose)
      messerr(">>> computing simulation %d", isimu + 1);
    if (simu->simulate(ns, 0, verbose, cov0, nd)) return 1;
    if (simu->compute(dbout, iuid + isimu * nvar, verbose)) return 1;

    // Modify the name of the output
    for (Id ivar = 0; ivar < nvar; ivar++)
    {
      String ps(prefix);
      ps.append(delim + "V" + std::to_string(ivar + 1));
      ps.append(delim + "S" + std::to_string(isimu + 1));
      namconvS.setPrefix(ps);
      namconvS.setNamesAndLocators(dbin, VectorString(), ELoc::Z, 1, dbout, iuid + isimu * nvar + ivar, "", 1);
    }
  }
  return 0;
}
} // namespace gstlrn