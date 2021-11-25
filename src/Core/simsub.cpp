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
#include "Basic/Law.hpp"
#include "geoslib_e.h"
#include "geoslib_old_f.h"

/*! \cond */

#define TRANS(i,j)     (trans[(j) + nfacies * (i)])
#define PROPS(i)       (props[i])
#define PROLD(i)       (propold[i])

/*! \endcond */

/****************************************************************************/
/*!
**  Check if the transition matrix is irreductible
**
** \return  1 if the transition matrix is not irreductible; 0 otherwise
**
** \param[in]  nfacies  Number of facies
** \param[in]  verbose  Verbose option
** \param[in]  trans    Transition matrix
**
*****************************************************************************/
static int st_check_irreductibility(int     nfacies,
                                    int     verbose,
                                    double *trans)
{
  int    *flag,i,j,ndeb,nend,error;
  double  total;

  /* Initializations */

  error = 1;
  flag  = nullptr;

  /* Check that the transition matrix is correct */
  
  for (i=0; i<nfacies; i++)
  {
    total = 0.;
    for (j=0; j<nfacies; j++)
    {
      if (TRANS(i,j) < 0. || TRANS(i,j) > 1.) goto label_end;
      total += TRANS(i,j);
    }
    if (total <= 0.) goto label_end;
    for (j=0; j<nfacies; j++) TRANS(i,j) /= total;
  }

  /* Check the irreductibility */

  flag = (int *) mem_alloc(sizeof(int) * nfacies,1);
  flag[0] = nend = ndeb = 0;
  for (i=1; i<nfacies; i++)
  {
    flag[i] = 0;
    if (TRANS(i,0) > 0)
    {
      flag[i] = 1;
      nend++;
    }
  }

  while (ndeb != nend)
  {
    for (i=0; i<nfacies; i++)
      if (flag[i])
        for (j=0; j<nfacies; j++)
          if (i != j && TRANS(j,i) > 0) flag[j] = 1;
    ndeb = nend;
    for (i=nend=0; i<nfacies; i++) nend += flag[i];
  }
  if (nend != nfacies) goto label_end;

  /* Printout (conditional) */

  if (verbose)
    print_matrix("Transitions",0,1,nfacies,nfacies,NULL,trans);

  /* Set the error return code */

  error = 0;

label_end:
  flag = (int *) mem_free((char *) flag);
  if (error) messerr("The transition matrix is not irreductible");
  return(error);
}

/****************************************************************************/
/*!
**  Derive proportions from the transition matrix
**
** \return The proprotion matrix
**
** \param[in]  nfacies  Number of facies
** \param[in]  verbose  Verbose option
** \param[in]  trans    Transition matrix
**
** \remarks The output proportion matrix must be freed afterwards 
**
*****************************************************************************/
static double *trans_to_props(int     nfacies,
                              int     verbose,
                              double *trans)
{
  double *props,*propold,diff,w0,val,total;
  int     i,j,flag_error;
  static  double eps = 1.e-5;

  /* Initializations */

  props = propold = nullptr;
  if (nfacies <= 0 || trans == nullptr) return(props);
  props   = (double *) mem_alloc(sizeof(double) * nfacies,1);
  propold = (double *) mem_alloc(sizeof(double) * nfacies,1);

  /* Checks the transition matrix */

  for (i=flag_error=0; i<nfacies && flag_error == 0; i++)
  {
    total = 0.;
    for (j=0; j<nfacies; j++)
      total += ABS(TRANS(i,j));

    if (total <= 0.) 
      flag_error = 1;
    else
      for (j=0; j<nfacies; j++) 
        TRANS(i,j) = ABS(TRANS(i,j)) / total;
  }

  /* Wrong transition matrix: initialize it */

  if (flag_error)
    for (i=0; i<nfacies; i++)
      for (j=0; j<nfacies; j++) 
        TRANS(i,j) = 1./ nfacies;
  
  /* Initialize the proportion matrix */
  
  for (i=0; i<nfacies; i++) PROPS(i) = 1. / nfacies;
  
  /* Loop to reach the stationarity of the propotions */

  diff = 2. * eps;
  while (diff > eps)
  {
    w0 = diff = 0.;
    for (i=0; i<nfacies; i++) PROLD(i) = PROPS(i);
    for (i=0; i<nfacies; i++)
    {
      val = 0.;
      for (j=0; j<nfacies; j++) val += TRANS(j,i) * PROLD(j);
      PROPS(i) = val;
      w0 += val;
    }
    if (w0 == 0.) w0 = 1.;
    for (i=0; i<nfacies; i++)
    {
      PROPS(i) /= w0;
      diff += ABS(PROLD(i) - PROPS(i));
    }
  }

  /* Core deallocation */

  propold = (double *) mem_free((char *) propold);

  /* Printout the proportions */

  if (verbose)
    print_matrix("Proportions",0,1,1,nfacies,NULL,props);

  return(props);
}

/*****************************************************************************
 **
 ** Checks the validity of an orientation factor
 **
 ** \returns 1 if the vector is not valid; 0 otherwise
 **
 ** \param[in]  factor      Disorientation factor
 ** \param[in]  verbose     Verbose option
 **
 ** \param[out]  factor     Possibly corrected Disorientation factor
 **
 *****************************************************************************/
static int st_check_factor(double *factor,
                           int     verbose)

{
  if ((*factor) < 0.) 
  {
    if (verbose) 
    {
      messerr("The desorientation factor cannot be negative");
      return(1);
    }
    (*factor) = 0.;
  }
  if ((*factor) > 1.)
  {
    if (verbose)
    {
      messerr("The desorientation factor cannot be larger than 1");
      return(1);
    }
    (*factor) = 1.;
  }
  return(0);
}

/*****************************************************************************
 **
 ** Checks the validity of an orientation vector
 **
 ** \returns 1 if the vector is not valid; 0 otherwise
 **
 ** \param[in]  vector      Disorientation vector
 ** \param[in]  verbose     Verbose option
 **
 ** \param[out] vector      Disorientation vector (normalized)
 **
 *****************************************************************************/
static int st_check_orientation(double *vector,
                                int     verbose)
{
  int i;
  double total;

  total = 0.;
  for (i=0; i<3; i++) total += vector[i] * vector[i];
  if (total <= 0.) 
  {
    if (verbose)
    {
      messerr("The desorientation vector should not be zero");
      return(1);
    }
    vector[0] = total = 1.;
  }
  for (i=0; i<3; i++) vector[i] /= sqrt(total);
  return(0);
}

/*****************************************************************************
 **
 ** Calculate the projected value
 **
 ** \param[in,out]  plan    SubPlan structure
 ** \param[in]  factor      Disorientation factor
 ** \param[in]  vector      Disorientation vector
 **
 *****************************************************************************/
static void st_calcul_value(SubPlan& plan,
                            double   factor,
                            double  *vector)
{
  int    ival,i;
  double cossin;

  ival = ((2. * plan.rndval) > (1. + factor)) ? -1 : 1;
  cossin = 0.;
  for (i=0; i<3; i++) cossin += plan.coor[i] * vector[i];
  if (cossin < 0) ival = -ival;
  plan.value = (double) ival;
}

/*****************************************************************************
 **
 ** Generate a simulation on a regular 3D grid using substitution method
 **
 ** \returns Error return code
 **
 ** \param[in]  dbgrid      Db structure (should be a grid)
 ** \param[in]  seed        Seed
 ** \param[in]  nfacies     Number of facies
 ** \param[in]  nstates     Number of states
 ** \param[in]  flag_direct 1 to perform the Directing step
 ** \param[in]  flag_coding 1 to perform the Coding step
 ** \param[in]  flag_orient 1 if disorientation must be used
 ** \param[in]  flag_auto   1 for an automatic number of states
 ** \param[in]  intensity   Intensity of the Poisson Process
 ** \param[in]  factor      Disorientation factor within [0,1]
 ** \param[in]  vector      Disorientation vector
 ** \param[in]  trans       Transition matrix (Dimension: nfacies * nfacies)
 ** \param[in]  verbose     Verbose option
 **
 *****************************************************************************/
GSTLEARN_EXPORT int substitution(Db      *dbgrid,
                             int      seed,
                             int      nfacies,
                             int      nstates,
                             int      flag_direct,
                             int      flag_coding,
                             int      flag_orient,
                             int      flag_auto,
                             double   intensity,
                             double   factor,
                             double   vector[3],
                             double  *trans,
                             int      colfac,
                             int      colang[3],
                             int      verbose)
{
  SubPlanes *splanes;
  double     cen[3],*props,w0,p0,u,prod,valloc,valtot,vmin,vmax,value,diagonal;
  int       *status,*indg,flag_local,flag_angloc;
  int        i,ip,np,ie,je,error,ndim,iptr,iech,idim,ival;

  /* Initializations */

  law_set_random_seed(seed);
  error    = 1;
  status   = indg = nullptr;
  props    = nullptr;
  flag_local = flag_angloc = np = 0;
  splanes  = nullptr;

  /* Preliminary checks */

  ndim = dbgrid->getNDim();
  if (! is_grid(dbgrid) || ndim > 3)
  {
    messerr("The substitution is available for Grid File with dimension <= 3");
    goto label_end;
  }
  if (flag_coding)
  {
    if (st_check_irreductibility(nfacies,verbose,trans)) goto label_end;
  }

  /* Check that the validity of the desorientation information */

  if (flag_orient)
  {
    for (i=0; i<3; i++) 
      if (colang[i] >= 0) flag_angloc = 1;

    /* Check the (constant) angle */

    if (! flag_angloc && st_check_orientation(vector,1)) goto label_end;

    /* Check the (constant) desorientation factor */

    if (colfac < 0 && st_check_factor(&factor,1)) goto label_end;

    flag_local = (flag_angloc || colfac >= 0);
  }

  /* Core allocation */

  indg = db_indg_alloc(dbgrid);

  /* Add the attributes for storing the results */

  iptr = dbgrid->addFields(1,0.);
  if (iptr < 0) goto label_end;

  /***********************/
  /* Information process */
  /***********************/

  if (flag_direct)
  {

    /* Calculate the number of planes */

    if (db_extension_diag(dbgrid,&diagonal)) goto label_end;
    np = law_poisson(diagonal * intensity * GV_PI);
    if (np <= 0) goto label_end;

    /* Generate the Poisson planes */

    splanes = poisson_manage_planes(1,np,splanes);
    if (splanes == nullptr) goto label_end;
    if (poisson_generate_planes(dbgrid,splanes)) goto label_end;
    
    /* Assigning a value to the half-space that contains the center */
    
    for (ip=0; ip<np; ip++)
    {
      SubPlan& plan = splanes->plans[ip];
      if (! flag_orient)
      {
        ival = (plan.rndval > 0.5) ? -1 : 1;
        plan.value = ival;
      }
      else if (! flag_local)
      {
        st_calcul_value(plan,factor,vector);
      }
    }

    /* Simulating the directing function */

    for (iech=0; iech<dbgrid->getSampleNumber(); iech++)
    {
      db_index_sample_to_grid(dbgrid,iech,indg);
      for (idim=0; idim<ndim; idim++) cen[idim] = dbgrid->getCoordinate(iech,idim);

      /* Loop on the planes */
      
      valtot = 0.;
      for (ip=0; ip<np; ip++)
      {
        SubPlan& plan = splanes->plans[ip];
        prod = 0.;
        for (i=0; i<3; i++) prod += plan.coor[i] * cen[i];

        if (flag_local)
        {
          if (colfac >= 0) 
          {
            factor = dbgrid->getArray(iech,colfac);
            (void) st_check_factor(&factor,0);
          }
          if (flag_angloc)
          {
            for (i=0; i<3; i++)
              if (colang[i] >= 0) vector[i] = dbgrid->getArray(iech,colang[i]);
            (void) st_check_orientation(vector,0);
          }	  
          st_calcul_value(plan,factor,vector);
        }

        valloc = plan.value / 2.;
        valtot += (prod + plan.intercept > 0) ? valloc : -valloc;
      }
      dbgrid->setArray(iech,iptr,valtot);
    }

    /* Core deallocation */

    splanes = poisson_manage_planes(-1,np,splanes);

    /* Printout statistics on the information process */

    if (verbose)
    {
      message("\n");
      message("Information Process:\n");
      message("Number of planes generated = %d\n",np);
    }
  }

  /***************************************/
  /* Determination of the extreme values */
  /***************************************/
  
  vmin =  1.e30;
  vmax = -1.e30;
  for (iech=0; iech<dbgrid->getSampleNumber(); iech++)
  {
    if (! dbgrid->isActive(iech)) continue;
    value = (flag_direct) ? 
      dbgrid->getArray(iech,iptr) : dbgrid->getVariable(iech,0);
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
      message("The number of states should be equal to %d\n",np);
      message("Nevertheless, your choice prevails\n");
    }
    props  = trans_to_props(nfacies,verbose,trans);
    status = (int *) mem_alloc(sizeof(int) * nstates,1);

    /* Simulation of the initial state */
    
    u = law_uniform(0.,1.);
    w0 = ie = 0;
    while (w0 < u) w0 += props[ie++];
    status[0] = ie - 1;
    
    /* Simulation of the current state */
    
    for (ip=1; ip<nstates; ip++)
    {
      u  = law_uniform(0.,1.);
      p0 = je = 0;
      ie = status[ip-1];
      while (p0 < u)
      {
        p0 += TRANS(ie,je);
        je++;
      }
      status[ip] = je - 1;
    }

    /* Simulating the directing function */
    
    for (iech=0; iech<dbgrid->getSampleNumber(); iech++)
    {
      if (! dbgrid->isActive(iech)) continue;
      value = (flag_direct) ? 
        dbgrid->getArray(iech,iptr) : dbgrid->getVariable(iech,0);
      ival = (int) ((value - vmin) / (vmax - vmin) * nstates);
      if (ival < 0)        ival = 0;
      if (ival >= nstates) ival = nstates - 1;
      dbgrid->setArray(iech,iptr,1 + status[ival]);
    }

    /* Core deallocation */

    status = (int    *) mem_free((char *) status);
    props  = (double *) mem_free((char *) props);

    /* Printout statistics */

    if (verbose)
    {
      message("\nCoding process: \n");
      message("Number of coded states     = %d \n",nstates);
      message("Minimum information value  = %lf\n", vmin);
      message("Maximum information value  = %lf\n", vmax);
    }
  }

  /* Set the error return code */

  error = 0;

  /* Core deallocation */

label_end:
  splanes = poisson_manage_planes(-1,np,splanes);
  indg = db_indg_free(indg);
  status = (int    *) mem_free((char *) status);
  props  = (double *) mem_free((char *) props);
  return(error);
}
