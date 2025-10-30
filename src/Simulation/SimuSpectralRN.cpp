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
#include "Simulation/SimuSpectralRN.hpp"
#include "Basic/Law.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Simulation/ASimuSpectral.hpp"
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

  // simulation of the spatial frequencies
  if (cov0 == nullptr)
  {
    _omega = _cova->simulateSpectralOmega(ns);
  }
  else
  {
    _omega = cov0->simulateSpectralOmega(ns);
  }
  // simulation of the normalizing coefficients
  _gamma.reset(ns, nvar);
  for (Id ib = 0; ib < ns; ib++)
  {
    VectorDouble values(nvar);
    double gamma = sqrt(-log(law_uniform()) * 2 * nvar / ns);
    // importance sampling or multivariate model
    if ((nvar > 1) | (cov0 != nullptr))
    {
      VectorDouble freq = _omega.getRow(ib);
      MatrixSymmetric H(nvar);
      for (Id ivar = 0; ivar < nvar; ivar++)
      {
        for (Id jvar = 0; jvar <= ivar; jvar++)
        {
          double ratioIS = _cova->evalSpectrumRatio(freq, ivar, jvar, cov0);
          H.setValue(ivar, jvar, ratioIS);
        }
      }
      // square root of the symmetric matrix
      if (H.computeSquareRoot(H) != 0)
      {
        message("Error in computing square root matrix\n");
        return 1;
      }
      Id icol        = law_int_uniform(0, nvar - 1);
      VectorDouble v = H.getColumn(icol);
      for (Id ivar = 0; ivar < nvar; ivar++)
        values[ivar] = gamma * v[ivar];
    }
    else
    {
      for (Id ivar = 0; ivar < nvar; ivar++)
        values[ivar] = gamma;
    }
    _gamma.setRow(ib, values);
  } // loop over the spectral components
  _isPrepared = true;
  return 0;
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