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
#include "Simulation/ASimulation.hpp"
#include "Simulation/SimuPartition.hpp"
#include "Simulation/SimuPartitionParam.hpp"
#include "Basic/Law.hpp"

#include <math.h>

#define COOR(iech,idim)    (coor[(iech) * ndim + (idim)])

SimuPartition::SimuPartition(int nbsimu, int seed)
    : ASimulation(nbsimu, seed)
{
}

SimuPartition::SimuPartition(const SimuPartition &r)
    : ASimulation(r)
{
}

SimuPartition& SimuPartition::operator=(const SimuPartition &r)
{
  if (this != &r)
  {
    ASimulation::operator=(r);
  }
  return *this;
}

SimuPartition::~SimuPartition()
{
}

/*****************************************************************************
 **
 ** Generate a simulation on a regular 3D grid using Voronoi Mosaic Model
 **
 ** \returns Error return code
 **
 ** \param[in]  dbgrid      Db structure (should be a grid)
 ** \param[in]  model       Model used for the valuation of tesselation
 ** \param[in]  dilate      Array of dilation radius (optional)
 ** \param[in]  seed        Seed
 ** \param[in]  intensity   Intensity of the Poisson Process
 ** \param[in]  nbtuba      Number of bands (for the gaussian field simulation)
 ** \param[in]  verbose     Verbose option
 **
 *****************************************************************************/
int SimuPartition::voronoi(DbGrid *dbgrid,
                           Model *model,
                           const SimuPartitionParam& parparam,
                           int iptr,
                           int verbose)
{
  law_set_random_seed(getSeed());
  int ndim = dbgrid->getNDim();
  VectorDouble simgrid(dbgrid->getSampleNumber());

  /************************************/
  /* Simulation of the Gaussian field */
  /************************************/

  double volume = 1.;
  VectorDouble field(ndim);
  VectorDouble origin(ndim);
  for (int i = 0; i < ndim; i++)
  {
    double dil = parparam.getDilate(i);
    field[i] = dbgrid->getDX(i) * dbgrid->getNX(i);
    origin[i] = dbgrid->getX0(i) - dbgrid->getDX(i) / 2. - dil * field[i] / 2.;
    field[i] = field[i] * (1. + dil);
    volume *= field[i];
  }

  /* Derive the number of points */

  int nbpoints = (int) (volume * parparam.getIntensity());
  if (verbose)
    message("Boolean simulation. Intensity = %lf - Nb. seeds = %d\n",
            parparam.getIntensity(), nbpoints);

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
  if (simtub(NULL, dbpoint, model, NULL, 1,
             getSeed(), parparam.getNbtuba())) return 1;

  /* Expand the data values over the grid nodes */

  int iattp = dbpoint->getColumnNumber() - 1;
  if (expand_point_to_grid(dbpoint, dbgrid,
                           iattp, -1, 0, -1, -1, -1, -1, 0,
                           VectorDouble(), simgrid)) return 1;

  /* Save the grid in dbgrid */

  dbgrid->setColumnByUID(simgrid, iptr);

  delete dbpoint;
  return 0;
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
int SimuPartition::poisson(DbGrid *dbgrid,
                           Model *model,
                           const SimuPartitionParam& parparam,
                           int iptr,
                           int verbose)
{
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

  if (simtub(NULL, dbgrid, model, NULL, 1, getSeed(),
             parparam.getNbtuba())) return 1;
  iattg = dbgrid->getColumnNumber() - 1;

  /***********************/
  /* Information process */
  /***********************/

  /* Calculate the number of planes */

  double diagonal = dbgrid->getExtensionDiagonal();
  np = law_poisson(diagonal * parparam.getIntensity() * GV_PI);
  if (np <= 0) return 1;

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
      for (int i = 0; i < 3; i++)
        prod += planes[ip].getCoor(i) * cen[i];
      valtot += (prod + planes[ip].getIntercept() > 0) ?
          planes[ip].getRndval() : -planes[ip].getRndval();
    }
    dbgrid->setArray(iech, iptr, valtot);
  }

  /* Printout statistics on the information process */

  if (verbose) message("Number of planes generated = %d\n", np);

  /******************/
  /* Coding process */
  /******************/

  for (int iech = 0; iech < dbgrid->getSampleNumber(); iech++)
  {
    if (!dbgrid->isActive(iech)) continue;
    double valref = dbgrid->getArray(iech, iptr);
    if (FFFF(valref)) continue;

    /* Check if the value has already been treated */

    double valsim = _stackSearch(stacks, valref);
    if (! FFFF(valsim))
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

  return 0;
}

double SimuPartition::_stackSearch(const std::vector<Stack>& stacks,
                                   double valref)
{
  for (int i = 0; i < (int) stacks.size(); i++)
  {
    if (stacks[i].valref == valref) return stacks[i].valsim;
  }
  return TEST;
}
