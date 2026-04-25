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
#include "Basic/Message.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
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
  SimuSpectralRN::SimuSpectralRN(Id nbsimu, Id ns, Id nd, Id seed, bool verbose)
    : CalcSimuSpectral(nbsimu, ns, nd, seed, verbose)
    , _gamma()
    , _omega()
    , _sp()
  {
  }

  SimuSpectralRN::~SimuSpectralRN()
  {
    delete _sp;
  }

  /**
   * Simulate the spectrum components for Rn
   */
  Id SimuSpectralRN::_simulate(Id isimu)
  {
    law_set_random_seed(getSeedPerSimu(isimu));

    const ACov* cov = getModelGeneric()->getCov();
    if (cov == nullptr) return -1;

    // Optional printout
    if (getVerbose())
    {
      message("Simulation of the spectrum\n");
      message("- Space dimension   = R%d\n", _getNDim());
      message("- Number of variables  = %d\n", _getNVar());
      message("- Number of spectral components = %d\n", _getNs());
    }
    delete _sp;
    _sp = cov->simulateOnRN(_getNs());
    return 1;
  }

  /**
   * Compute the simulation on Dbout using Spectral Method
   *
   * @param db Db containing the results
   * @param activeArray Array of booleans indicating the active samples in dbout
   * @param tab Array for storing one (multivariate) simulation on 'dbout'
   */
  Id SimuSpectralRN::_compute(
    Db* db,
    const VectorBool& activeArray,
    VectorVectorDouble& tab)
  {
    if (_sp == nullptr)
    {
      messerr("SpectrumOnRN not initialized.\n");
      return 1;
    }
    // Optional printout
    if (getVerbose())
    {
      message("Spectral Simulation on a set of Isolated Points\n");
      message("- Number of samples = %d\n", db->getNSample());
    }
    _sp->compute(db, activeArray, tab);
    return 0;
  }

} // namespace gstlrn
