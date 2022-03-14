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
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Model/Model.hpp"

#include <math.h>

/*! \cond */

#define COOR(iech,idim)    (coor[(iech) * ndim + (idim)])

/*! \endcond */

typedef struct
{
  int nitem;
  VectorDouble valref;
  VectorDouble valsim;
} Stack;

/****************************************************************************/
/*!
 **  Generate the Poisson planes that cover the grid
 **
 ** \return Error return code
 **
 ** \param[in]  dbgrid   Db corresponding to the target grid
 **
 ** \param[out] splanes  SubPlanes structure
 **
 ** \remarks  The array 'planes' contains successively a,b,c,d such that
 ** \remarks  ax + by + cz + d = 0
 ** \remarks  The valuation of each line is assigned a uniform value [0,1]
 **
 *****************************************************************************/
int poisson_generate_planes(DbGrid *dbgrid, SubPlanes *splanes)
{
  double ap[3], diagonal;

  /* Determine the extension of the grid */

  int ndim = dbgrid->getNDim();
  VectorDouble mini(ndim);
  VectorDouble maxi(ndim);
  db_extension(dbgrid, mini, maxi);
  if (db_extension_diag(dbgrid, &diagonal)) return (1);

  /* Loop on the planes to be generated */

  for (int ip = 0; ip < splanes->nplan; ip++)
  {
    SubPlan &plan = splanes->plans[ip];
    double d0 = diagonal * law_uniform(-1., 1.) / 2.;
    double u = 0.;
    for (int idim = 0; idim < 3; idim++)
    {
      ap[idim] = law_gaussian();
      u += ap[idim] * ap[idim];
    }
    u = sqrt(u);
    for (int idim = 0; idim < 3; idim++)
    {
      ap[idim] /= u;
      d0 -= ap[idim] * (mini[idim] + maxi[idim]) / 2.;
    }
    if (d0 < 0)
    {
      for (int idim = 0; idim < 3; idim++)
        ap[idim] = -ap[idim];
      d0 = -d0;
    }

    /* Storing the plane */

    for (int idim = 0; idim < 3; idim++)
      plan.coor[idim] = ap[idim];
    plan.intercept = d0;
    plan.rndval = law_uniform(0., 1.);
  }
  return (0);
}

/*****************************************************************************
 **
 ** Management of the SubPlanes structure
 **
 ** \returns Pointer to the SubPlanes structure
 **
 ** \param[in]  mode        1 for allocation and -1 for deallocation
 ** \param[in]  np          Number of planes (for allocation only)
 ** \param[in]  splanes     SubPlanes structure to be deallocated
 **
 *****************************************************************************/
SubPlanes* poisson_manage_planes(int mode,
                                                 int np,
                                                 SubPlanes *splanes)
{

  /* Dispatch */

  if (mode > 0)
  {

    /* Allocation */

    splanes = new SubPlanes;
    splanes->plans.resize(np);
    splanes->nplan = np;
    for (int ip = 0; ip < np; ip++)
    {
      SubPlan &plan = splanes->plans[ip];
      for (int idim = 0; idim < 3; idim++)
        plan.coor[idim] = 0.;
      plan.intercept = 0.;
      plan.value = 0.;
      plan.rndval = 0.;
    }
    return splanes;
  }
  else
  {

    /* Deallocation */

    if (splanes == nullptr) return (splanes);
    delete splanes;
    return nullptr;
  }
}

/*****************************************************************************
 **
 ** Management of the Lookup Table
 **
 ** \returns Pointer to the Lookup Table
 **
 ** \param[in]  stack    Stack structure
 ** \param[in]  mode     Type of operation
 **                      1: initialize; 0: add an item; -1: free
 ** \param[in]  valref   Value of the tessellation
 ** \param[in]  valsim   Value of the gaussian field
 **
 *****************************************************************************/
static void st_manage_stack(Stack *stack,
                            int mode,
                            double valref,
                            double valsim)
{
  int nitem;

  /* Dispatch */

  switch (mode)
  {
    case 1: /* Creation */
      stack->nitem = 0;
      break;

    case 0: /* Adding an item */
      nitem = stack->nitem + 1;
      stack->valref.resize(nitem);
      stack->valsim.resize(nitem);
      stack->valref[stack->nitem] = valref;
      stack->valsim[stack->nitem] = valsim;
      stack->nitem = nitem;
      break;

    case -1: /* Free */
      stack->nitem = 0;
      break;
  }
}

/*****************************************************************************
 **
 ** Search for the value if already registered in the stack
 **
 ** \returns 1 if the value has already been stacked and 0 otherwise
 **
 ** \param[in]  stack       Stack structure
 ** \param[in]  valref      Value to be searched
 **
 ** \param[out] valsim      Simulated value (if matching is found)
 **
 *****************************************************************************/
static int st_stack_search(Stack *stack, double valref, double *valsim)
{
  int i;

  (*valsim) = TEST;
  for (i = 0; i < stack->nitem; i++)
    if (stack->valref[i] == valref)
    {
      (*valsim) = stack->valsim[i];
      return (1);
    }
  return (0);
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
int tessellation_poisson(DbGrid *dbgrid,
                         Model *model,
                         int seed,
                         double intensity,
                         int nbtuba,
                         int verbose)
{
  SubPlanes *splanes;
  Stack stack;
  double cen[3], prod, valtot, diagonal, valref, valsim;
  int *status, *indg;
  int i, ip, np, error, ndim, iatts, iattg, iech, idim;

  /* Initializations */

  law_set_random_seed(seed);
  error = 1;
  np = 0;
  iatts = iattg = -1;
  status = indg = nullptr;
  splanes = nullptr;
  st_manage_stack(&stack, 1, 0., 0.);

  /* Preliminary checks */

  ndim = dbgrid->getNDim();
  if (!is_grid(dbgrid))
  {
    messerr("The output Db file must be a grid");
    goto label_end;
  }
  if (!is_grid(dbgrid) || ndim > 3)
  {
    messerr(
        "The Poisson Tesselation is available for Grid File with dimension <= 3");
    goto label_end;
  }

  /* Add the attributes for storing the results */

  iatts = dbgrid->addColumnsByConstant(1, TEST);
  if (iatts < 0) goto label_end;
  indg = db_indg_alloc(dbgrid);

  /************************************/
  /* Simulation of the Gaussian field */
  /************************************/

  if (simtub(NULL, dbgrid, model, NULL, 1, seed, nbtuba)) goto label_end;
  iattg = dbgrid->getColumnNumber() - 1;

  /***********************/
  /* Information process */
  /***********************/

  /* Calculate the number of planes */

  if (db_extension_diag(dbgrid, &diagonal)) goto label_end;
  np = law_poisson(diagonal * intensity * GV_PI);
  if (np <= 0) goto label_end;

  /* Generate the Poisson planes */

  splanes = poisson_manage_planes(1, np, splanes);
  if (splanes == nullptr) goto label_end;
  if (poisson_generate_planes(dbgrid, splanes)) goto label_end;

  /* Assigning a value to the half-space that contains the center */

  for (ip = 0; ip < np; ip++)
  {
    SubPlan &plan = splanes->plans[ip];
    plan.value = (plan.rndval > 0.5) ? -1 :
                                       1;
    ;
  }

  /* Simulating the directing function */

  for (iech = 0; iech < dbgrid->getSampleNumber(); iech++)
  {
    if (!dbgrid->isActive(iech)) continue;
    db_index_sample_to_grid(dbgrid, iech, indg);
    for (idim = 0; idim < ndim; idim++)
      cen[idim] = dbgrid->getCoordinate(iech, idim);

    /* Loop on the planes */

    valtot = 0.;
    for (ip = 0; ip < np; ip++)
    {
      SubPlan &plan = splanes->plans[ip];
      prod = 0.;
      for (i = 0; i < 3; i++)
        prod += plan.coor[i] * cen[i];
      valtot += (prod + plan.intercept > 0) ? plan.rndval :
                                              -plan.rndval;
    }
    dbgrid->setArray(iech, iatts, valtot);
  }

  /* Core deallocation */

  splanes = poisson_manage_planes(-1, np, splanes);

  /* Printout statistics on the information process */

  if (verbose) message("Number of planes generated = %d\n", np);

  /******************/
  /* Coding process */
  /******************/

  for (iech = 0; iech < dbgrid->getSampleNumber(); iech++)
  {
    if (!dbgrid->isActive(iech)) continue;
    valref = dbgrid->getArray(iech, iatts);
    if (FFFF(valref)) continue;

    /* Check if the value has already been treated */

    if (!st_stack_search(&stack, valref, &valsim))
    {

      /* Not in the stack: read the value from the Gaussian field */

      valsim = dbgrid->getArray(iech, iattg);

      /* Add this new value to the stack */

      st_manage_stack(&stack, 0, valref, valsim);
    }

    /* Substitute the values in the gaussian field */

    dbgrid->setArray(iech, iattg, valsim);
  }

  /* Set the error return code */

  error = 0;

  /* Core deallocation */

  label_end: if (iatts >= 0) dbgrid->deleteColumnByUID(iatts);
  splanes = poisson_manage_planes(-1, np, splanes);
  st_manage_stack(&stack, -1, 0., 0.);
  indg = db_indg_free(indg);
  status = (int*) mem_free((char* ) status);
  return (error);
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
int tessellation_voronoi(DbGrid *dbgrid,
                         Model *model,
                         double *dilate,
                         int seed,
                         double intensity,
                         int nbtuba,
                         int verbose)
{
  int i, ip, nbpoints, error, iatts, iattp, ndim;
  double volume, origin[2], field[3], dil, oriloc, dimloc;
  Db *dbpoint;
  VectorDouble coor, simgrid, simpoint;

  /* Initializations */

  law_set_random_seed(seed);
  error = 1;
  iatts = iattp = -1;
  dbpoint = nullptr;

  /* Preliminary checks */

  ndim = dbgrid->getNDim();
  if (ndim > 3)
  {
    messerr(
        "The Poisson Tesselation is available for Grid File with dimension <= 3");
    goto label_end;
  }

  /* Add the attributes for storing the results */

  iatts = dbgrid->addColumnsByConstant(1, TEST);
  if (iatts < 0) goto label_end;
  simgrid.resize(dbgrid->getSampleNumber());

  /************************************/
  /* Simulation of the Gaussian field */
  /************************************/

  volume = 1.;
  for (i = 0; i < ndim; i++)
  {
    dil = (dilate != nullptr) ? dilate[i] :
                                0.;
    field[i] = dbgrid->getDX(i) * dbgrid->getNX(i);
    origin[i] = dbgrid->getX0(i) - dbgrid->getDX(i) / 2. - dil * field[i] / 2.;
    field[i] = field[i] * (1. + dil);
    volume *= field[i];
  }

  /* Derive the number of points */

  nbpoints = (int) (volume * intensity);
  if (verbose)
    message("Boolean simulation. Intensity = %lf - Nb. seeds = %d\n", intensity,
            nbpoints);

  /* Core allocation */

  coor.resize(nbpoints * ndim);

  /* Simulate the uniform points */

  for (i = 0; i < ndim; i++)
  {
    oriloc = origin[i];
    dimloc = field[i];
    for (ip = 0; ip < nbpoints; ip++)
      COOR(ip,i) = oriloc + dimloc * law_uniform(0., 1.);
  }

  /* Create the Point Data Base */

  dbpoint = db_create_point(nbpoints, ndim, ELoadBy::SAMPLE, 1, coor);
  dbpoint->setLocatorsByUID(ndim, 0, ELoc::X);
  simpoint.resize(dbpoint->getSampleNumber());

  /* Perform the simulation at the seed points */

  if (simtub(NULL, dbpoint, model, NULL, 1, seed, nbtuba)) goto label_end;

  /* Expand the data values over the grid nodes */

  iattp = dbpoint->getColumnNumber() - 1;
  if (expand_point_to_grid(dbpoint, dbgrid, iattp, -1, 0, -1, -1, -1, -1, 0,
                           VectorDouble(), simgrid)) goto label_end;

  /* Save the grid in dbgrid */

  dbgrid->setColumnByUID(simgrid, iatts);

  /* Set the error return code */

  error = 0;

  label_end: dbpoint = db_delete(dbpoint);
  return (error);
}

