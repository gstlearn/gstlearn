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
#include "Boolean/AShape.hpp"
#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"
#include "Simulation/ASimulation.hpp"
#include "Simulation/SimuSubstitution.hpp"
#include "Simulation/SimuSubstitutionParam.hpp"
#include "Basic/Law.hpp"

#include <math.h>

SimuSubstitution::SimuSubstitution(int nbsimu, int seed)
    : ASimulation(nbsimu, seed),
      AStringable()
{
}

SimuSubstitution::SimuSubstitution(const SimuSubstitution &r)
    : ASimulation(r),
      AStringable(r)
{
}

SimuSubstitution& SimuSubstitution::operator=(const SimuSubstitution &r)
{
  if (this != &r)
  {
    ASimulation::operator=(r);
    AStringable::operator =(r);
  }
  return *this;
}

SimuSubstitution::~SimuSubstitution()
{
}

String SimuSubstitution::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  return sstr.str();
}

/*****************************************************************************
 **
 ** Generate a simulation on a regular 3D grid using substitution method
 **
 ** \returns Error return code
 **
 ** \param[in]  dbgrid      Db structure (should be a grid)
 ** \param[in]  subparam    SimuSubstitutionParam structure
 ** \param[in]  iptr        Vector for storing the simulation
 ** \param[in]  verbose     Verbose option
 **
 *****************************************************************************/
int SimuSubstitution::simulate(DbGrid *dbgrid,
                               const SimuSubstitutionParam& subparam,
                               int iptr,
                               int verbose)
{
  SubPlanes *splanes;

  /* Initializations */

  law_set_random_seed(getSeed());
  int error = 1;
  int status = 0;
  VectorDouble props;
  bool flag_local = false;
  bool flag_angloc = false;
  int np = 0;

  /* Preliminary checks */

  if (subparam.isFlagCoding())
  {
    if (st_check_irreductibility(subparam.getNfacies(), verbose,
                                 subparam.getTrans())) goto label_end;
  }

  /* Check that the validity of the desorientation information */

  if (subparam.isFlagOrient())
  {
    for (int i = 0; i < 3; i++)
      if (subparam.getColang(i) >= 0) flag_angloc = 1;

    /* Check the (constant) angle */

    if (!flag_angloc &&
        st_check_orientation(subparam.getVector(), 1)) goto label_end;

    /* Check the (constant) desorientation factor */

    if (subparam.getColfac() < 0 &&
        st_check_factor(&subparam.getFactor(), 1)) goto label_end;

    flag_local = (flag_angloc || subparam.getColfac() >= 0);
  }

  /***********************/
  /* Information process */
  /***********************/

  if (subparam.isFlagDirect())
  {

    /* Calculate the number of planes */

    double diagonal dbgrid->getExtensionDiagonal();
    np = law_poisson(diagonal * subparam.getIntensity() * GV_PI);
    if (np <= 0) goto label_end;

    /* Generate the Poisson planes */

    splanes = poisson_manage_planes(1, np, splanes);
    if (splanes == nullptr) goto label_end;
    if (poisson_generate_planes(dbgrid, splanes)) goto label_end;

    /* Assigning a value to the half-space that contains the center */

    for (int ip = 0; ip < np; ip++)
    {
      SubPlan &plan = splanes->plans[ip];
      if (! subparam.isFlagOrient())
      {
        int ival = (plan.rndval > 0.5) ? -1 : 1;
        plan.value = ival;
      }
      else if (! flag_local)
      {
        st_calcul_value(plan, subparam.getFactor(), subparam.getVector());
      }
    }

    /* Simulating the directing function */

    for (int iech = 0; iech < dbgrid->getSampleNumber(); iech++)
    {
      db_index_sample_to_grid(dbgrid, iech, indg);
      for (idim = 0; idim < ndim; idim++)
        cen[idim] = dbgrid->getCoordinate(iech, idim);

      /* Loop on the planes */

      double valtot = 0.;
      for (int ip = 0; ip < np; ip++)
      {
        SubPlan &plan = splanes->plans[ip];
        double prod = 0.;
        for (int i = 0; i < 3; i++)
          prod += plan.coor[i] * cen[i];

        if (flag_local)
        {
          if (subparam.getColfac() >= 0)
          {
            factor = dbgrid->getArray(iech, colfac);
            (void) st_check_factor(&factor, 0);
          }
          if (flag_angloc)
          {
            for (i = 0; i < 3; i++)
              if (colang[i] >= 0) vector[i] = dbgrid->getArray(iech, colang[i]);
            (void) st_check_orientation(vector, 0);
          }
          st_calcul_value(plan, factor, vector);
        }

        valloc = plan.value / 2.;
        valtot += (prod + plan.intercept > 0) ? valloc :
                                                -valloc;
      }
      dbgrid->setArray(iech, iptr, valtot);
    }

    /* Core deallocation */

    splanes = poisson_manage_planes(-1, np, splanes);

    /* Printout statistics on the information process */

    if (verbose)
    {
      message("\n");
      message("Information Process:\n");
      message("Number of planes generated = %d\n", np);
    }
  }

  /***************************************/
  /* Determination of the extreme values */
  /***************************************/

  vmin = 1.e30;
  vmax = -1.e30;
  for (iech = 0; iech < dbgrid->getSampleNumber(); iech++)
  {
    if (!dbgrid->isActive(iech)) continue;
    value = (flag_direct) ? dbgrid->getArray(iech, iptr) :
                            dbgrid->getVariable(iech, 0);
    if (value < vmin) vmin = value;
    if (value > vmax) vmax = value;
  }
  if (vmin > vmax)
  {
    messerr("No Direction Function has been coded");
    messerr("before the Coding Process takes place");
    goto label_end;
  }
  np = (int) (vmax - vmin + 0.5);

  /******************/
  /* Coding process */
  /******************/

  if (flag_coding)
  {
    if (flag_auto) nstates = np;
    if (flag_direct && nstates != np)
    {
      message("You have used the internal information process\n");
      message("The number of states should be equal to %d\n", np);
      message("Nevertheless, your choice prevails\n");
    }
    props = trans_to_props(nfacies, verbose, trans);
    status = (int*) mem_alloc(sizeof(int) * nstates, 1);

    /* Simulation of the initial state */

    u = law_uniform(0., 1.);
    w0 = ie = 0;
    while (w0 < u)
      w0 += props[ie++];
    status[0] = ie - 1;

    /* Simulation of the current state */

    for (ip = 1; ip < nstates; ip++)
    {
      u = law_uniform(0., 1.);
      p0 = je = 0;
      ie = status[ip - 1];
      while (p0 < u)
      {
        p0 += TRANS(ie, je);
        je++;
      }
      status[ip] = je - 1;
    }

    /* Simulating the directing function */

    for (iech = 0; iech < dbgrid->getSampleNumber(); iech++)
    {
      if (!dbgrid->isActive(iech)) continue;
      value = (flag_direct) ? dbgrid->getArray(iech, iptr) :
                              dbgrid->getVariable(iech, 0);
      ival = (int) ((value - vmin) / (vmax - vmin) * nstates);
      if (ival < 0) ival = 0;
      if (ival >= nstates) ival = nstates - 1;
      dbgrid->setArray(iech, iptr, 1 + status[ival]);
    }

    /* Core deallocation */

    status = (int*) mem_free((char* ) status);
    props = (double*) mem_free((char* ) props);

    /* Printout statistics */

    if (verbose)
    {
      message("\nCoding process: \n");
      message("Number of coded states     = %d \n", nstates);
      message("Minimum information value  = %lf\n", vmin);
      message("Maximum information value  = %lf\n", vmax);
    }
  }

  /* Set the error return code */

  error = 0;

  /* Core deallocation */

  label_end: splanes = poisson_manage_planes(-1, np, splanes);
  indg = db_indg_free(indg);
  status = (int*) mem_free((char* ) status);
  props = (double*) mem_free((char* ) props);
  return (error);
}
