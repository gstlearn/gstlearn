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
#include "Basic/Message.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "Enum/ESimuType.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Model/Model.hpp"
#include "Simulation/CalcSimuSpectral.hpp"
#include "Simulation/SpectrumOnRN.hpp"
#include "Stats/Classical.hpp"
#include "geoslib_define.h"
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
  , _sp()
  , _cov0(cov0)
{
}

SimuSpectralRN::~SimuSpectralRN()
{
  delete _sp;
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
  DECLARE_UNUSED(cova)
  const ACov* cov = getModelGeneric()->getCov();
  if (cov == nullptr)
  {
    messerr("Covariance model not defined.");
    return -1;
  }
  if (!cov->isValidForSimulation(ESimuType::SPECTRAL))
  {
    messerr("Covariance not valid for spectral simulation.");
    return -2;
  }

  // Optional printout
  if (getVerbose())
  {
    message("Simulation of the spectrum\n");
    message("- Space dimension   = R%d\n", _getNDim());
    message("- Number of variables  = %d\n", _getNVar());
    message("- Number of spectral components = %d\n", _getNs());
    if (_cov0 != nullptr)
      message("Simulation using importance sampling\n");
  }
  delete _sp;
  _sp = cov->simulateOnRN(_getNs());
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
  auto nech = dbout->getNSample();
  if (_sp == nullptr)
  {
    messerr("SpectrumOnRN not initialized.\n");
    return 1;
  }
  // Optional printout
  if (getVerbose())
  {
    message("Spectral Simulation on a set of Isolated Points\n");
    message("- Number of samples = %d\n", nech);
  }
  _sp->compute(dbout, activeArray, tab);
  return 0;
}

} // namespace gstlrn