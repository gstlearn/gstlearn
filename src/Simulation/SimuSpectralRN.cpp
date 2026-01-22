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
#include "Basic/VectorNumT.hpp"
#include "Db/Db.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Model/Model.hpp"
#include "Simulation/ASimuSpectral.hpp"
#include "Simulation/SimuSpectralRN.hpp"
#include "Simulation/SpectrumRN.hpp"
#include "Stats/Classical.hpp"

#include <cmath>

namespace gstlrn
{

/**
 * ---------------------------------
 * Spectral simulation on Rn
 * ---------------------------------
 */
SimuSpectralRN::SimuSpectralRN(const ACov* cova)
  : ASimuSpectral(cova)
  , _gamma()
  , _omega()
{
}

SimuSpectralRN::SimuSpectralRN(const SimuSpectralRN& r)
  : ASimuSpectral(r)
  , _gamma(r._gamma)
  , _omega(r._omega)
{
}

SimuSpectralRN& SimuSpectralRN::operator=(const SimuSpectralRN& r)
{
  if (this != &r)
  {
    _cova       = r._cova;
    _isPrepared = r._isPrepared;
    _phi        = r._phi;
    _gamma      = r._gamma;
    _omega      = r._omega;
  }
  return *this;
}

SimuSpectralRN::~SimuSpectralRN()
{
}

/**
 * Simulate the spectrum components for Rn
 *
 * @param ns Number of components
 * @param nd not used (for S2 only)
 * @param cov0 the auxiliary covariance function used for importance sampling
 * @param verbose Verbose flag
 */
Id SimuSpectralRN::_simulate(Id ns,
                             Id nd,
                             const ACov* cov0,
                             bool verbose)
{
  DECLARE_UNUSED(nd)
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
  Id ndim = _cova->getNDim();
  Id nvar = _cova->getNVar();
  Id ierr = 0;
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

  // Cleaning any previously allocated memory
  _gamma.reset(0, 0);
  _omega.reset(0, 0);

  // Simulation of the spectrum by the covariance
  SpectrumRN sp = _cova->simulateSpectrumRN(ns, cov0);
  if (sp.getNs() > 0)
  {
    _gamma      = sp.getGamma();
    _omega      = sp.getOmega();
    _isPrepared = true;
    ierr        = 0;
  }
  else
  {
    _isPrepared = false;
    ierr        = 1;
  }
  return ierr;
}

/**
 * Compute the simulation on Dbout using Spectral Method
 *
 * @param dbout Db containing the results
 * @param iuid  Address for storage (or 0 if the variable must be created locally)
 * @param verbose Verbose flag
 */
Id SimuSpectralRN::_compute(Db* dbout, Id iuid, bool verbose)
{
  if (!_isPrepared)
  {
    messerr("You should run 'simulate' beforehand");
    return 1;
  }
  Id nvar = _cova->getNVar();
  Id ndim = dbout->getNDim();
  VectorDouble coor(ndim);
  VectorInt ranks;
  dbout->getSampleRanksPerVariable(ranks);
  Id ns   = getNs();
  Id nech = ranks.length();
  if (nech <= 0)
  {
    messerr("'dbout' must have a positive number of active samples");
    return 1;
  }
  // Optional printout
  if (verbose)
  {
    message("Spectral Simulation on a set of Isolated Points\n");
    message("- Number of samples = %d\n", nech);
    message("- Space dimension   = %d\n", ndim);
    message("- Number of spectral components = %d\n", ns);
    message("- Number of variables = %d\n", nvar);
  }
  // Loop on the active samples
  for (Id jech = 0; jech < nech; jech++)
  {
    Id iech = ranks[jech];
    dbout->getCoordinatesInPlace(coor, iech);
    VectorDouble u = _omega.prodMatVec(coor);

    for (Id ib = 0; ib < ns; ib++)
      u[ib] = cos(u[ib] + _phi[ib]);
    VectorDouble values = _gamma.prodMatVec(u, true);
    for (Id iv = 0; iv < nvar; iv++)
      dbout->setArray(iech, iuid + iv, values[iv]);
  }
  return 0;
}

} // namespace gstlrn