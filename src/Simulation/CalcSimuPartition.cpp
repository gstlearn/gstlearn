/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "Boolean/AShape.hpp"
#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"
#include "Simulation/ACalcSimulation.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Simulation/CalcSimuPartition.hpp"
#include "Simulation/SimuPartitionParam.hpp"
#include "Basic/Law.hpp"

#include <math.h>

#define COOR(iech,idim)    (coor[(iech) * ndim + (idim)])

CalcSimuPartition::CalcSimuPartition(int mode,
                                     int nbsimu,
                                     int seed,
                                     bool verbose)
    : ACalcSimulation(nbsimu, seed),
      _mode(mode),
      _verbose(verbose),
      _iattOut(-1),
      _parparam()
{
}

CalcSimuPartition::~CalcSimuPartition()
{
}

/*****************************************************************************
 **
 ** Generate a simulation on a regular 3D grid using Voronoi Mosaic Model
 **
 ** \returns Error return code
 **
 *****************************************************************************/
bool CalcSimuPartition::_voronoi()
{
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
  int ndim = _getNDim();
  VectorDouble simgrid(dbgrid->getSampleNumber());

  /************************************/
  /* Simulation of the Gaussian field */
  /************************************/

  double volume = 1.;
  VectorDouble field(ndim);
  VectorDouble origin(ndim);
  for (int i = 0; i < ndim; i++)
  {
    double dil = _parparam.getDilate(i);
    field[i] = dbgrid->getDX(i) * dbgrid->getNX(i);
    origin[i] = dbgrid->getX0(i) - dbgrid->getDX(i) / 2. - dil * field[i] / 2.;
    field[i] = field[i] * (1. + dil);
    volume *= field[i];
  }

  /* Derive the number of points */

  int nbpoints = (int) (volume * _parparam.getIntensity());
  if (_verbose)
    message("Boolean simulation. Intensity = %lf - Nb. seeds = %d\n",
            _parparam.getIntensity(), nbpoints);

  /* Simulate the uniform points */

  VectorDouble coor(nbpoints * ndim, 0.);
  for (int idim = 0; idim < ndim; idim++)
    for (int ip = 0; ip < nbpoints; ip++)
      COOR(ip,idim) = origin[idim] + field[idim] * law_uniform(0., 1.);

  /* Create the Point Data Base */

  Db* dbpoint = Db::createFromSamples(nbpoints, ELoadBy::SAMPLE, coor);
  dbpoint->setLocatorsByUID(ndim, 0, ELoc::X);
  VectorDouble simpoint(dbpoint->getSampleNumber());

  /* Perform the simulation at the seed points */
  if (simtub(NULL, dbpoint, getModel(), NULL, 1,
             getSeed(), _parparam.getNbtuba())) return false;

  /* Expand the data values over the grid nodes */

  int iattp = dbpoint->getColumnNumber() - 1;
  if (expand_point_to_grid(dbpoint, dbgrid,
                           iattp, -1, 0, -1, -1, -1, -1, 0,
                           VectorDouble(), simgrid)) return 1;

  /* Save the grid in dbgrid */

  dbgrid->setColumnByUID(simgrid, _iattOut);

  delete dbpoint;
  return true;
}

/*****************************************************************************
 **
 ** Generate a simulation on a regular 3D grid using Poisson Polyhedra Model
 **
 ** \returns Error return code
 **
 ** \param[in]  dbgrid      Db structure (should be a grid)
 ** \param[in]  model       Model used for the valuation of tesselation
 ** \param[in]  seed        Seed
 ** \param[in]  intensity   Intensity of the Poisson Process
 ** \param[in]  nbtuba      Number of bands (for the gaussian field simulation)
 ** \param[in]  verbose     Verbose option
 **
 *****************************************************************************/
bool CalcSimuPartition::_poisson()
{
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());

  std::vector<Stack> stacks;
  std::vector<Plane> planes;

  law_set_random_seed(getSeed());
  int np = 0;
  int iattg = -1;

  /* Preliminary checks */

  int ndim = dbgrid->getNDim();

  /************************************/
  /* Simulation of the Gaussian field */
  /************************************/

  if (simtub(NULL, dbgrid, getModel(), NULL, 1, getSeed(), _parparam.getNbtuba()))
    return false;
  iattg = dbgrid->getColumnNumber() - 1;

  /***********************/
  /* Information process */
  /***********************/

  /* Calculate the number of planes */

  double diagonal = dbgrid->getExtensionDiagonal();
  np = law_poisson(diagonal * _parparam.getIntensity() * GV_PI);
  if (np <= 0) return false;

  /* Generate the Poisson planes */

  planes = Plane::poissonPlanesGenerate(dbgrid, np);

  /* Assigning a value to the half-space that contains the center */

  for (int ip = 0; ip < np; ip++)
    planes[ip].setValue( (planes[ip].getRndval() > 0.5) ? -1 : 1);

  /* Simulating the directing function */

  VectorDouble cen(ndim);
  for (int iech = 0; iech < dbgrid->getSampleNumber(); iech++)
  {
    if (!dbgrid->isActive(iech)) continue;
    dbgrid->getCoordinatesInPlace(iech, cen);

    /* Loop on the planes */

    double valtot = 0.;
    for (int ip = 0; ip < np; ip++)
    {
      double prod = 0.;
      for (int i = 0; i < (int) cen.size(); i++)
        prod += planes[ip].getCoor(i) * cen[i];
      valtot += (prod + planes[ip].getIntercept() > 0) ?
          planes[ip].getRndval() : -planes[ip].getRndval();
    }
    dbgrid->setArray(iech, _iattOut, valtot);
  }

  /* Printout statistics on the information process */

  if (_verbose)
    message("Number of planes generated = %d\n", np);

  /******************/
  /* Coding process */
  /******************/

  for (int iech = 0; iech < dbgrid->getSampleNumber(); iech++)
  {
    if (!dbgrid->isActive(iech)) continue;
    double valref = dbgrid->getArray(iech, _iattOut);
    if (FFFF(valref)) continue;

    /* Check if the value has already been treated */

    double valsim = _stackSearch(stacks, valref);
    if (FFFF(valsim))
    {

      /* Not in the stack: read the value from the Gaussian field */

      valsim = dbgrid->getArray(iech, iattg);

      /* Add this new value to the stack */

      Stack stack;
      stack.valref = valref;
      stack.valsim = valsim;
      stacks.push_back(stack);
    }

    /* Substitute the values in the Gaussian field */

    dbgrid->setArray(iech, iattg, valsim);
  }

  // Delete the internal Simulation

  dbgrid->deleteColumnByUID(iattg);

  return true;
}

double CalcSimuPartition::_stackSearch(const std::vector<Stack>& stacks,
                                   double valref)
{
  for (int i = 0; i < (int) stacks.size(); i++)
  {
    if (stacks[i].valref == valref) return stacks[i].valsim;
  }
  return TEST;
}

bool CalcSimuPartition::_check()
{
  if (! ACalcSimulation::_check()) return false;

  if (! hasDbout())
  {
    messerr("The argument 'dbout' must be defined");
    return false;
  }
  if (! hasModel())
  {
    messerr("The argument 'model' must be defined");
    return false;
  }
  int ndim = _getNDim();
  if (ndim > 3)
  {
    messerr("The Turning Band Method is not a relevant simulation model");
    messerr("for this Space Dimension (%d)", ndim);
    return false;
  }
  if (! getDbout()->isGrid())
  {
    messerr("The argument 'dbout'  should be a grid");
    return false;
  }
  if (_mode != 1 && _mode != 2)
  {
    messerr("Argument 'mode'(%d) should be:");
    messerr(" 1 for Voronoi Tesselation");
    messerr(" 2 for Poisson Hyperplanes");
    return false;
  }
  return true;
}

bool CalcSimuPartition::_preprocess()
{
    _iattOut = _addVariableDb(2, 1, ELoc::SIMU, 1);
    if (_iattOut < 0) return false;
    return true;
}

bool CalcSimuPartition::_run()
{
  law_set_random_seed(getSeed());

  // Dispatch

  if (_mode == 1)
    return (_voronoi());
  else
    return (_poisson());
}

bool CalcSimuPartition::_postprocess()
{
  _renameVariable(ELoc::Z, 1, _iattOut, String(), getNbSimu());
  return true;
}

void CalcSimuPartition::_rollback()
{
  _cleanVariableDb(1);
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
int tessellation_poisson(DbGrid *dbgrid,
                         Model *model,
                         const SimuPartitionParam& parparam,
                         int seed,
                         int verbose,
                         const NamingConvention& namconv)
{
  CalcSimuPartition simpart(2, 1, seed, verbose);
  simpart.setDbout(dbgrid);
  simpart.setModel(model);
  simpart.setNamingConvention(namconv);
  simpart.setParparam(parparam);

  int error = (simpart.run()) ? 0 : 1;
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
int tessellation_voronoi(DbGrid *dbgrid,
                         Model *model,
                         const SimuPartitionParam& parparam,
                         int seed,
                         int verbose,
                         const NamingConvention& namconv)
{
  CalcSimuPartition simpart(1, 1, seed, verbose);
  simpart.setDbout(dbgrid);
  simpart.setModel(model);
  simpart.setNamingConvention(namconv);
  simpart.setParparam(parparam);

  int error = (simpart.run()) ? 0 : 1;
  return error;
}

