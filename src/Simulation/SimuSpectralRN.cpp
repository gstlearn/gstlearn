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
#include "Basic/VectorNumT.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Model/Model.hpp"
#include "Simulation/CalcSimuSpectral.hpp"
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
SimuSpectralRN::SimuSpectralRN(Id nbsimu, Id ns, Id nd, Id seed, const ACov* cov0, bool verbose)
  : CalcSimuSpectral(nbsimu, ns, nd, seed, verbose)
  , _gamma()
  , _omega()
  , _cov0(cov0)
{
}

SimuSpectralRN::~SimuSpectralRN()
{
}

bool SimuSpectralRN::_check()
{
  if (!CalcSimuSpectral::_check()) return false;

  bool hasCov0 = (_cov0 != nullptr);
  if (hasCov0)
  {
    if (!_cov0->isValidForSimulation(ESimuType::SPECTRAL))
    {
      messerr("Simulation of the harmonic components is not implemented for the auxiliary covariance");
      return false;
    }
    if (_cov0->getNVar() > 1)
    {
      messerr("The auxiliary covariance should be scalar");
      return false;
    }
  }

  return true;
}

/**
 * Simulate the spectrum components for Rn
 */
Id SimuSpectralRN::_simulate(const ACov* cova)
{
  Id ns   = _getNs();
  Id ndim = _getNDim();
  Id nvar = _getNVar();

  // Optional printout
  if (getVerbose())
  {
    message("Simulation of the spectrum\n");
    message("- Space dimension   = R%d\n", ndim);
    message("- Number of variables  = %d\n", nvar);
    message("- Number of spectral components = %d\n", ns);
    if (_cov0 != nullptr)
      message("Simulation using importance sampling\n");
  }

  // Cleaning any previously allocated memory
  _gamma.reset(0, 0);
  _omega.reset(0, 0);

  // Simulation of the spectrum by the covariance
  SpectrumRN sp = cova->simulateSpectrumRN(ns, _cov0);

  if (sp.getNs() > 0)
  {
    _gamma = sp.getGamma();
    _omega = sp.getOmega();
    return 0;
  }
  return 1;
}

/**
 * Compute the simulation on Dbout using Spectral Method
 *
 * @param dbout Db containing the results
 * @param activeArray Array of booleans indicating the active samples in dbout
 * @param tab Array for storing one (multivariate) simulation on 'dbout'
 */
Id SimuSpectralRN::_compute(Db* dbout, const VectorBool& activeArray, VectorVectorDouble& tab)
{
  auto nvar = _getNVar();
  auto ns   = _getNs();
  auto ndim = dbout->getNDim();
  auto nech = dbout->getNSample();
  VectorDouble coor(ndim);

  // Optional printout
  if (getVerbose())
  {
    message("Spectral Simulation on a set of Isolated Points\n");
    message("- Number of samples = %d\n", nech);
  }

  // Loop on the active samples
  VectorDouble u(ns);
  VectorDouble values(nvar);
  for (Id iech = 0; iech < nech; iech++)
  {
    if (!activeArray[iech]) continue;
    dbout->getCoordinatesInPlace(coor, iech);
    AMatrix::productInPlace(u, _omega, coor);
    for (Id ib = 0; ib < ns; ib++)
      u[ib] = cos(u[ib] + getPhi(ib));
    AMatrix::productInPlace(values, _gamma, u, true);
    for (Id ivar = 0; ivar < nvar; ivar++)
      tab[ivar][iech] = values[ivar];
  }
  return 0;
}

} // namespace gstlrn