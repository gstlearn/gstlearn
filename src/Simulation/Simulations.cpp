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
#include "Simulation/Simulations.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Neigh/ANeigh.hpp"
#include "Simulation/CalcSimuEden.hpp"
#include "Simulation/CalcSimuFFT.hpp"
#include "Simulation/CalcSimuPartition.hpp"
#include "Simulation/CalcSimuRefine.hpp"
#include "Simulation/CalcSimuSpectral.hpp"
#include "Simulation/CalcSimuSubstitution.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Simulation/SimuSpectralRN.hpp"
#include "Simulation/SimuSpectralS2.hpp"

namespace gstlrn
{

  /*****************************************************************************/
  /*!
  **  Multivariate multiphase propagation into a set of components
  **  constrained by initial conditions and fluid densities
  **
  ** \return  Error return code : 1 no fluid to propagate
  **
  ** \param[in]  dbgrid        Db grid structure
  ** \param[in]  name_facies   Name of the variable containing the Facies
  ** \param[in]  name_fluid    Name of the variable containing the Fluid
  ** \param[in]  name_perm     Name of the variable containing the Permeability
  ** \param[in]  name_poro     Name of the variable containing the Porosity
  ** \param[in]  nfacies       number of facies (facies 0 excluded)
  ** \param[in]  nfluids       number of fluids
  ** \param[in]  niter         Number of iterations
  ** \param[in]  speeds        Array containing the travel speeds
  ** \param[in]  show_fluid    1 for modifying the value of the cells to show
  ** \li                       the initial valid fluid information
  ** \li                       the cork (different from shale)
  ** \param[in]  number_max    Maximum count of cells invaded (or TEST)
  ** \param[in]  volume_max    Maximum volume invaded (or TEST)
  ** \param[in]  seed          Seed for random number generator (or 0)
  ** \param[in]  verbose       1 for a verbose option
  ** \param[in]  namconv       Naming convention
  **
  ** \remark  Directions are ordered as follows :
  ** \remark  0: +X; 1: -X; 2: +Y; 3: -Y; 4: +Z(up); 5: -Z(down)
  ** \remark  The coding of the matrix is:
  ** \remark              facies + nfacies * fluid
  ** \remark  Facies: 0 (Shale), 1 to nfacies, -1 (Cork)
  ** \remark  Fluids: 0 (undefined), 1 to nfluids, -1 (No Fluid)
  ** \remark  Fluids should be ordered by increasing weight
  ** \remark  A Permeability variable is a value (>=1) which divides
  ** \remark  the velocities. This variable is optional.
  ** \remark  A Porosity variable is a value (in [0,1]) which multiplies
  ** \remark  the velocities. This variable is optional.
  **
  ** \remark  the volumes. This variable is optional.
  ** \remark  Volume_max represents the volumic part of the invaded area:
  ** \remark  it is always <= number of cells invaded.
  **
  *****************************************************************************/
  Id fluidPropagation(
    DbGrid* dbgrid,
    const String& name_facies,
    const String& name_fluid,
    const String& name_perm,
    const String& name_poro,
    Id nfacies,
    Id nfluids,
    Id niter,
    const VectorInt& speeds,
    bool show_fluid,
    double number_max,
    double volume_max,
    Id seed,
    bool verbose,
    const NamingConvention& namconv)
  {
    CalcSimuEden seden(nfacies, nfluids, niter, 1, seed, verbose);

    seden.setDbout(dbgrid);
    seden.setNamingConvention(namconv);

    seden.setIndFacies(dbgrid->getUID(name_facies));
    seden.setIndFluid(dbgrid->getUID(name_fluid));
    if (!name_poro.empty()) seden.setIndPoro(dbgrid->getUID(name_poro));
    if (!name_perm.empty()) seden.setIndPerm(dbgrid->getUID(name_perm));

    seden.setSpeeds(speeds);
    seden.setShowFluid(show_fluid);
    seden.setNMax(number_max);
    seden.setVolumeMax(volume_max);

    // Run the calculator
    Id error = (seden.run()) ? 0 : 1;
    return error;
  }

  /****************************************************************************/
  /*!
   **  Perform the non-conditional simulation by FFT method on a grid
   **
   ** \return  Error return code
   **
   ** \param[in]  dbin    Input Db structure
   ** \param[in]  dbout   Output Db Grid structure
   ** \param[in]  model   ModelGeneric structure
   ** \param[in]  param   SimuFFTParam structure
   ** \param[in]  nbsimu  Number of simulations
   ** \param[in]  seed    Value of the seed
   ** \param[in]  verbose Verbose flag
   ** \param[in]  namconv Naming Convention
   **
   *****************************************************************************/
  Id simuFFT(
    Db* dbin,
    DbGrid* dbout,
    ModelGeneric* model,
    SimuFFTParam& param,
    Id nbsimu,
    Id seed,
    bool verbose,
    const NamingConvention& namconv)
  {
    CalcSimuFFT simufft(nbsimu, verbose, seed);
    simufft.setDbin(dbin);
    simufft.setDbout(dbout);
    simufft.setModelGeneric(model);
    simufft.setNamingConvention(namconv);
    simufft.setParam(param);

    Id error = (simufft.run()) ? 0 : 1;
    return error;
  }

  /****************************************************************************/
  /*!
   **  Calculate the change of support coefficients by FFT method
   **  in the lognormal case on a grid
   **
   ** \return  r^2 coefficients for the different logarithmic variances
   **
   ** \param[in]  db      Db structure
   ** \param[in]  model   ModelGeneric structure
   ** \param[in]  param   SimuFFTParam structure
   ** \param[in]  sigma   Array of logarithmic variances
   ** \param[in]  seed    Seed for random number generator
   ** \param[in]  verbose Verbose flag
   **
   *****************************************************************************/
  VectorDouble getChangeSupport(
    DbGrid* db,
    ModelGeneric* model,
    const SimuFFTParam& param,
    const VectorDouble& sigma,
    Id seed,
    bool verbose)
  {
    CalcSimuFFT simufft(1, verbose, seed);
    simufft.setDbout(db);
    simufft.setModelGeneric(model);
    simufft.setParam(param);
    return simufft.changeSupport(sigma);
  }

  /*****************************************************************************
   **
   ** Generate a simulation on a regular 3D grid using Poisson Polyhedra Model
   **
   ** \returns Error return code
   **
   ** \param[in]  dbgrid      Db structure (should be a grid)
   ** \param[in]  model       Model used for the valuation of tesselation
   ** \param[in]  parparam    SimuPartitionParam structure
   ** \param[in]  seed        Seed
   ** \param[in]  verbose     Verbose option
   ** \param[in]  namconv     Naming Convention
   **
   *****************************************************************************/
  Id tessellationPoisson(
    DbGrid* dbgrid,
    Model* model,
    const SimuPartitionParam& parparam,
    Id seed,
    bool verbose,
    const NamingConvention& namconv)
  {
    CalcSimuPartition simpart(2, 1, seed, verbose);
    simpart.setDbout(dbgrid);
    simpart.setModelGeneric(model);
    simpart.setNamingConvention(namconv);
    simpart.setParparam(parparam);

    Id error = (simpart.run()) ? 0 : 1;
    return error;
  }

  /*****************************************************************************
   **
   ** Generate a simulation on a regular 3D grid using Voronoi Mosaic Model
   **
   ** \returns Error return code
   **
   ** \param[in]  dbgrid      Db structure (should be a grid)
   ** \param[in]  model       Model used for the valuation of tesselation
   ** \param[in]  parparam    SimuPartitionParam structure
   ** \param[in]  seed        Seed
   ** \param[in]  verbose     Verbose option
   ** \param[in]  namconv     Naming Convention
   **
   *****************************************************************************/
  Id tessellationVoronoi(
    DbGrid* dbgrid,
    Model* model,
    const SimuPartitionParam& parparam,
    Id seed,
    bool verbose,
    const NamingConvention& namconv)
  {
    CalcSimuPartition simpart(1, 1, seed, verbose);
    simpart.setDbout(dbgrid);
    simpart.setModelGeneric(model);
    simpart.setNamingConvention(namconv);
    simpart.setParparam(parparam);

    Id error = (simpart.run()) ? 0 : 1;
    return error;
  }

  /****************************************************************************/
  /*!
   **  Refine the simulation
   **
   ** \return  Newly refined Grid.
   **
   ** \param[in]  dbin       Input grid Db structure
   ** \param[in]  model      Model structure
   ** \param[in]  param      SimuRefineParam structure
   ** \param[in]  seed       Seed for the random number generator
   ** \param[in]  namconv    Naming convention
   **
   ** \remark For each dimension of the space, if N stands for the number of
   ** \remark nodes in the input grid, the number of nodes of the output grid
   ** \remark will be (N-1) * 2^p + 1 where p is the param.getNmult()
   **
   *****************************************************************************/
  DbGrid* simuRefine(
    DbGrid* dbin,
    Model* model,
    const SimuRefineParam& param,
    Id seed,
    const NamingConvention& namconv)
  {
    CalcSimuRefine simfine(1, seed);
    simfine.setDbin(dbin);
    simfine.setModelGeneric(model);
    simfine.setNamingConvention(namconv);
    simfine.setParam(param);

    return (simfine.run()) ? simfine.getResultingGrid() : nullptr;
  }

  /**
   * Perform a series of simulations (on Rn or on the Sphere) using Spectral Method
   *
   * @param dbin Input Db where the conditioning data are read
   * @param dbout Output Db where the results are stored
   * @param model ModelGeneric structure
   * @param neigh Neighborhood structure (optional)
   * @param nbsimu Number of simulations processed simultaneously
   * @param seed Seed used for the Random number generator
   * @param ns Number of spectral components
   * @param nd Maximum number of spectral orders on S2
   * @param verbose Verbose flag
   * @param namconv Naming Convention
   *
   * @note The conditional version is not yet available
   */
  Id simuSpectral(
    Db* dbin,
    Db* dbout,
    ModelGeneric* model,
    ANeigh* neigh,
    Id nbsimu,
    Id seed,
    Id ns,
    Id nd,
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

    auto* neighLocal = ANeigh::createDefaultNeighborhood(neigh, dbin, dbout);

    // Instantiate the Calculator
    std::unique_ptr<CalcSimuSpectral> spectral;

    // Instantiate
    if (isSimuRN)
    {
      spectral =
        std::make_unique<SimuSpectralRN>(nbsimu, ns, nd, seed, verbose);
    }
    else
    {
      spectral =
        std::make_unique<SimuSpectralS2>(nbsimu, ns, nd, seed, verbose);
    }

    // Set the members of the Calculator
    spectral->setDbin(dbin);
    spectral->setDbout(dbout);
    spectral->setModelGeneric(model);
    spectral->setNeigh(neighLocal);
    spectral->setNamingConvention(namconv);

    // Run the calculator
    Id error = (spectral->run()) ? 0 : 1;
    if (neigh != neighLocal) delete neighLocal;
    return error;
  }

  /*****************************************************************************
   **
   ** Generate a simulation on a regular 3D grid using substitution method
   **
   ** \returns Error return code
   **
   ** \param[in]  dbgrid      Db structure (should be a grid)
   ** \param[in]  subparam    SimuSubstitutionParam structure
   ** \param[in]  seed        Seed
   ** \param[in]  verbose     Verbose option
   ** \param[in]  namconv     Naming convention
   **
   *****************************************************************************/
  Id simuSubstitution(
    DbGrid* dbgrid,
    SimuSubstitutionParam& subparam,
    Id seed,
    Id verbose,
    const NamingConvention& namconv)
  {
    CalcSimuSubstitution simsub(1, seed, verbose);
    simsub.setDbout(dbgrid);
    simsub.setNamingConvention(namconv);
    simsub.setSubparam(subparam);

    // Run the calculator
    Id error = (simsub.run()) ? 0 : 1;
    return error;
  }

  /****************************************************************************/
  /*!
   **  Perform the conditional or non-conditional simulation
   **
   ** \return  Error return code
   **
   ** \param[in]  dbin       Input Db structure (optional)
   ** \param[in]  dbout      Output Db structure
   ** \param[in]  model      Model structure
   ** \param[in]  neigh      ANeigh structure (optional)
   ** \param[in]  nbsimu     Number of simulations
   ** \param[in]  seed       Seed for random number generator
   ** \param[in]  nbtuba     Number of turning bands
   ** \param[in]  flag_dgm   1 for Direct Block Simulation
   ** \param[in]  box        Surrounding Box provided by the user (see remarks)
   ** \param[in]  namconv    Naming convention
   **
   ** \remark  The arguments 'dbin' and 'neigh' are optional: they must
   ** \remark  be defined only for conditional simulations
   **
   ** \remark The argument 'box' is optional: it must be defined only if the user
   ** \remark wants to provide a surrounding box for the Turning Bands simulation.
   ** \remark Otherwise it is calculated automatically as the box containing
   ** \remark both the input and output Db (when provided), parallel to main axes.
   **
   *****************************************************************************/
  Id simtub(
    Db* dbin,
    Db* dbout,
    Model* model,
    ANeigh* neigh,
    Id nbsimu,
    Id seed,
    Id nbtuba,
    bool flag_dgm,
    const VectorVectorDouble& box,
    const NamingConvention& namconv)
  {
    // Instantiate the Calculator
    CalcSimuTurningBands situba(nbsimu, nbtuba, seed);
    ANeigh* neighLocal = ANeigh::createDefaultNeighborhood(neigh, dbin, dbout);

    // Set the members of the Calculator
    situba.setDbin(dbin);
    situba.setDbout(dbout);
    situba.setModelGeneric(model);
    situba.setNeigh(neighLocal);
    situba.setNamingConvention(namconv);
    situba.setFlagDGM(flag_dgm);
    situba.setBox(box);

    // Run the calculator
    Id error = (situba.run()) ? 0 : 1;
    if (neigh != neighLocal) delete neighLocal;
    return error;
  }

  /****************************************************************************/
  /*!
   **  Perform the conditional or non-conditional simulation
   **  with Bayesian Drift
   **
   ** \return  Error return code
   **
   ** \param[in]  dbin       Input Db structure (optional)
   ** \param[in]  dbout      Output Db structure
   ** \param[in]  model      Model structure
   ** \param[in]  neigh      ANeigh structure (optional)
   ** \param[in]  nbsimu     Number of simulations
   ** \param[in]  seed       Seed for random number generator
   ** \param[in]  nbtuba     Number of turning bands
   ** \param[in]  namconv    Naming convention
   **
   ** \remark  The arguments 'dbin' and 'neigh' are optional: they must
   ** \remark  be defined only for conditional simulations
   **
   *****************************************************************************/
  Id simbayes(
    Db* dbin,
    Db* dbout,
    Model* model,
    ANeigh* neigh,
    Id nbsimu,
    Id seed,
    Id nbtuba,
    const NamingConvention& namconv)
  {
    // Instantiate the Calculator
    CalcSimuTurningBands situba(nbsimu, nbtuba, seed);

    // Set the members of the Calculator
    situba.setDbin(dbin);
    situba.setDbout(dbout);
    situba.setModelGeneric(model);
    situba.setNeigh(neigh);
    situba.setNamingConvention(namconv);
    situba.setFlagBayes(true);

    // Run the calculator
    Id error = (situba.run()) ? 0 : 1;
    return error;
  }

} // namespace gstlrn
