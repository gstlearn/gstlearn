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
#include "geoslib_d.h"
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "csparse_f.h"

#include "Basic/Law.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/OptCst.hpp"
#include "Basic/ECst.hpp"
#include "Space/Space.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/ELoadBy.hpp"
#include "Space/ASpaceObject.hpp"
#include "Covariances/CovLMC.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovContext.hpp"
#include "Model/Model.hpp"

/*****************************************************************************/
/*!
**  Save the resulting vectors in 'dbgrid' and save it in a neutral file
**
** \return Error returned code
**
** \param[in]  consmin  Array giving the minimum bounds per pixel (optional)
** \param[in]  consmax  Array giving the maximum bounds per pixel (optional)
** \param[in]  z        Array giving the simulated field
**
*****************************************************************************/
static int st_save(Db    *dbgrid,
                   double *consmin,
                   double *consmax,
                   double *z)
{
  int nech,iptr;

  /* Initializations */
  
  nech = dbgrid->getSampleNumber();

  /* Add the terms to 'dbgrid' */
  
  if (db_locator_attribute_add(dbgrid,ELoc::Z,1,0,0.,&iptr)) return(1);
  for (int i=0; i<nech; i++) dbgrid->setArray(i,iptr,z[i]);
  if (consmin != (double *) NULL)
  {
    if (db_locator_attribute_add(dbgrid,ELoc::L,1,0,0.,&iptr)) return(1);
    for (int i=0; i<nech; i++) dbgrid->setArray(i,iptr,consmin[i]);
  }
  if (consmax != (double *) NULL)
  {
    if (db_locator_attribute_add(dbgrid,ELoc::U,1,0,0.,&iptr)) return(1);
    for (int i=0; i<nech; i++) dbgrid->setArray(i,iptr,consmax[i]);
  }

  /* Save the resulting 'dbgrid' in a neutral file */

  if (dbgrid->dumpToNF("Colored_Gibbs",1)) return(1);
  return(0);
}

/*****************************************************************************/
/*!
**  Print all the imported information (used for check only)
**
** \param[in]  colors   Array of colors (same dimension as 'z')
** \param[in]  consmin  Array giving the minimum bounds per pixel (optional)
** \param[in]  consmax  Array giving the maximum bounds per pixel (optional)
** \param[in]  sigma    Array giving the st. dev. per pixel
** \param[in]  Q        Precision matrix (for all colors)
**
*****************************************************************************/
static void st_print_all(int    *colors,
                         double *consmin,
                         double *consmax,
                         double *sigma,
                         cs     *Q)
{
  if (consmin != (double *) NULL)
    print_matrix ("consmin",0,0,1,10,NULL,consmin);
  if (consmax != (double *) NULL)
    print_matrix ("consmax",0,0,1,10,NULL,consmax);
  print_matrix ("sigma"  ,0,0,1,10,NULL,sigma);
  print_imatrix("colors" ,0,0,1,10,NULL,colors);
  cs_print_nice("Q",Q,10,10);
}

/*****************************************************************************/
/*!
**  Compressing the input vector 'z' into 'zred' by selecting the only pixels
**  where the value of 'colors' is equal to 'colref'
**
** \return  Number of selected pixels
**
** \param[in]  nvertex  Dimension of the input array 'z'
** \param[in]  colref   Value of the reference color
** \param[in]  z        Array of input values
** \param[in]  colors   Array of colors (same dimension as 'z')
**
** \param[out]  ind     Array which contains address in 'z' of items of 'zred'
** \param[out]  zred    Array of output values retained
**
*****************************************************************************/
static int st_vector_compress(int     nvertex,
                              int     colref,
                              double *z,
                              int    *colors,
                              int    *ind,
                              double *zred)
{
  int up,nup;
  
  up = nup = 0;
  for (int i=0; i<nvertex; i++)
  {
    if (colors[i] != colref) 
      zred[nup++] = z[i];
    else
      ind[up++] = i;
  }
  return(up);
}
  
/*****************************************************************************/
/*!
**  Draw a random value from a Gaussian distribution under the constraint
**  to be larger to a minimum bound
**
** \return  Returned Gaussian value
**
** \param[in]  iter     Rank of the current iteration
** \param[in]  niter    Total count of iterations
** \param[in]  valmin   Value of the minimum bound
** \param[in]  valmax   Value of the maximum bound
** \param[in]  mean     Mean of the Gaussian distribution
** \param[in]  sigma    Standard deviation of the Gaussian distribution
**
*****************************************************************************/
static double st_simcond(int    iter,
                         int    niter,
                         double valmin,
                         double valmax,
                         double mean,
                         double sigma)
{
  double locmin,locmax,x,delta,ratio;

  ratio  = (double) (niter-iter-1) / (double) (iter+1);
  delta  = valmax - valmin;
  if (valmin == -10.)
    locmin = -10.;
  else 
  {    
    locmin = valmin - delta * ratio;
    locmin = (locmin - mean) / sigma;
  }
  if (valmax == +10.)
    locmax = valmax;
  else
  {
    locmax = valmax + delta * ratio;
    locmax = (locmax - mean) / sigma;
  }

  x = law_gaussian_between_bounds(locmin,locmax);
  return(mean + x * sigma);
}

/*****************************************************************************/
/*!
**  Perform the Gibbs sampler algorithm on the contents of the vector 'z'
**  This algorithm treats simultaneously all pixels sharing the same 'color'
**  (due to Markovian property, their influence do not overlap)
**
** \return  Error returned code
**
** \param[in]  niter    Maximum number of Gibbs iterations
** \param[in]  ncolor   Total number of different colors (in 'colors')
** \param[in]  nvertex  Dimension of the arrays 'z', 'colors', 'cons', 'sigma'
** \param[in]  colors   Array giving the colors per pixel
** \param[in]  Qcols    Array of 'cs' structures (sparse matrices per color)
** \param[in]  consmin  Array giving the minimum bounds per pixel (optional)
** \param[in]  consmax  Array giving the maximum bounds per pixel (optional)
** \param[in]  sigma    Array giving the st. dev. per pixel
**
** \param[in/out] z     Array of gaussian values simulated
** \param[out] ind      Working array
** \param[out] krig     Working array
** \param[out] zred     Working array
**
** \remarks Arrays 'colors' and 'ind' are integer 
** \remarks        (Dimension 'nvertex')
** \remarks Arrays 'cons', 'sigma', 'z', 'krig' and 'zred' are double precision
** \remarks        (Dimension 'nvertex')
**
*****************************************************************************/
static int st_gibbs(int  niter,
                    int  ncolor,
                    int  nvertex,
                    int *colors,
                    cs **Qcols,
                    double *consmin,
                    double *consmax,
                    double *sigma,
                    double *z,
                    int    *ind,
                    double *krig,
                    double *zred)
{
  double valmin,valmax;
  int nc,i,ic;
  
  mestitle(1,"Entering in Gibbs algorithm with niter=%d and ncolor=%d",
           niter,ncolor);

  for (int iter=0; iter<niter; iter++)
  {
    if (iter % 1000 == 0) message("Iteration %d\n",iter);
    for (int icol=0; icol<ncolor; icol++)
    {
      nc = st_vector_compress(nvertex,icol,z,colors,ind,zred);
      cs_mulvec(Qcols[icol],nc,zred,krig);

      for (ic=0; ic<nc; ic++)
      {
        i      = ind[ic];
        valmin = (consmin != (double *) NULL) ? consmin[i] : -10.;
        valmax = (consmax != (double *) NULL) ? consmax[i] : +10.;
        z[i]   = st_simcond(iter,niter,valmin,valmax,krig[ic],sigma[i]);
      }
    }
  }
  return(0);
}

/****************************************************************************/
/*!
** Main Program
**
*****************************************************************************/
int main(int /*argc*/, char */*argv*/[])

{
  DbGrid   *dbgrid;
  Model    *model1,*model2;
  SPDE_Option    s_option;
  cs            *Q,**Qcols;

  double range_spde =   30.;
  double param_spde =    1.;
  double sill_spde  =    1.;
  double range_cons =   50.;
  double param_cons =    2.;
  double sill_cons  =    1.;

  VectorInt nx = { 100, 100 };
  VectorDouble x0 = { 0., 0. };
  VectorDouble dx = { 1., 1. };

  int    niter      = 10000;
  int    flag_print =     0;
  int    flag_save  =     1;
  const char triswitch[] = "nqQ";
  int     verbose, seed, ndim, nvar, iptr, nvertex, ncolor;
  int    *colors, *ind, rank;
  double *z, *krig, *zred, *consmin, *consmax, *sigma, diag;
  
  /***********************/
  /* 1 - Initializations */
  /***********************/

  dbgrid   = (DbGrid      *) NULL;
  model1   = (Model       *) NULL;
  model2   = (Model       *) NULL;
  colors   = (int         *) NULL;
  ind      = (int         *) NULL;
  z        = (double      *) NULL;
  krig     = (double      *) NULL;
  zred     = (double      *) NULL;
  consmin  = (double      *) NULL;
  consmax  = (double      *) NULL;
  sigma    = (double      *) NULL;
  Q        = (cs          *) NULL;
  Qcols    = (cs         **) NULL;
  verbose  = ncolor = 0;
  seed     = 31415;
  ndim     = 2;
  nvar     = 1;

  ASpaceObject::defineDefaultSpace(SPACE_RN, ndim);
  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Gibbs-");

  // Setup constants

  OptDbg::reset();
  law_set_random_seed(seed);
  OptCst::define(ECst::NTCAR,10.);
  OptCst::define(ECst::NTDEC,6.);
  
  // 2-D grid output file

  dbgrid = db_create_grid(0,ndim,0,ELoadBy::COLUMN,1,nx,x0,dx);
  db_locator_attribute_add(dbgrid,ELoc::X,ndim,0,0.,&iptr);
  db_grid_define_coordinates(dbgrid);
  db_extension_diag(dbgrid,&diag);
  CovContext ctxt(nvar,ndim,1000,diag);
    
  // Model for SPDE

  model1 = new Model(ctxt);
  CovLMC covs1(ctxt.getSpace());
  CovAniso cova1(ECov::BESSEL_K,range_spde,param_spde,sill_spde,ctxt);
  covs1.addCov(&cova1);
  model1->setCovList(&covs1);

  // Model for constraints

  model2 = new Model(ctxt);
  CovLMC covs2(ctxt.getSpace());
  CovAniso cova2(ECov::BESSEL_K,range_cons,param_cons,sill_cons,ctxt);
  covs2.addCov(&cova2);
  model2->setCovList(&covs2);

  // Creating the meshing for extracting Q

  s_option = spde_option_alloc();
  spde_option_update(s_option,triswitch);
  if (spde_check(NULL,dbgrid,model1,NULL,verbose,VectorDouble(),
                 1,1,0,0,0,0,0)) goto label_end;
  if (spde_prepar(NULL,dbgrid,VectorDouble(),s_option)) goto label_end;
  {
    SPDE_Matelem& Matelem = spde_get_current_matelem(0);
    Q = Matelem.QC->Q;
    if (Q == (cs *) NULL) goto label_end;
    nvertex = Matelem.s_mesh->nvertex;
  }

  // Create the color coding 

  colors = cs_color_coding(Q,0,&ncolor);
  if (colors == (int *) NULL) goto label_end;

  // Core allocation
  
  ind     = (int    *) mem_alloc(sizeof(int)    * nvertex,0);
  if (ind     == (int    *) NULL) goto label_end;
  z       = (double *) mem_alloc(sizeof(double) * nvertex,0);
  if (z       == (double *) NULL) goto label_end;
  krig    = (double *) mem_alloc(sizeof(double) * nvertex,0);
  if (krig    == (double *) NULL) goto label_end;
  zred    = (double *) mem_alloc(sizeof(double) * nvertex,0);
  if (zred    == (double *) NULL) goto label_end;
  consmin = (double *) mem_alloc(sizeof(double) * nvertex,0);
  if (consmin == (double *) NULL) goto label_end;
  consmax = (double *) mem_alloc(sizeof(double) * nvertex,0);
  if (consmax == (double *) NULL) goto label_end;

  // Creating the constraints

  if (spde_f(NULL, dbgrid, model2, VectorDouble(), s_option, 1, 1, seed, 2, 0,
             0, 0, 0, 0, 0, 0, 0)) goto label_end;
  rank  = dbgrid->getFieldNumber();
  for (int i=0; i<nvertex; i++)
  {
    consmin[i] = MIN(dbgrid->getArray(i,rank-1), dbgrid->getArray(i,rank-2));
    consmax[i] = MAX(dbgrid->getArray(i,rank-1), dbgrid->getArray(i,rank-2));
    z[i] = (consmin[i] + consmax[i]) / 2.;
  }
  
  // Creating the variance
  
  sigma = csd_extract_diag(Q,-2);
  if (sigma == (double *) NULL) goto label_end;

  // Scaling the Q matrix
  
  if (cs_scale(Q)) goto label_end;

  // Check the imported information

  if (flag_print) st_print_all(colors,consmin,consmax,sigma,Q);

  //----------------//
  // Main Algorithm //
  //----------------//

  Qcols = (cs **) mem_alloc(sizeof(cs *) * ncolor,1);
  for (int icol=0; icol<ncolor; icol++) Qcols[icol] = (cs *) NULL;
  for (int icol=0; icol<ncolor; icol++)
  {
    Qcols[icol] = cs_extract_submatrix_by_color(Q,colors,icol,1,0);
    if (Qcols[icol] == (cs *) NULL) goto label_end;
  }

  // Perform the Gibbs sampler
  if (st_gibbs(niter,ncolor,nvertex,colors,Qcols,consmin,consmax,sigma,
               z,ind,krig,zred)) goto label_end;

  // Add the newly created field to the grid for printout
  if (flag_save)
  {
    if (st_save(dbgrid,consmin,consmax,z)) goto label_end;
  }
  
label_end:
  for (int icol=0; icol<ncolor; icol++) Qcols[icol] = cs_spfree(Qcols[icol]);
  dbgrid   = db_delete(dbgrid);
  model1   = model_free(model1);
  model2   = model_free(model2);
  colors   = (int    *) mem_free((char *) colors);
  ind      = (int    *) mem_free((char *) ind);
  z        = (double *) mem_free((char *) z);
  krig     = (double *) mem_free((char *) krig);
  zred     = (double *) mem_free((char *) zred);
  consmin  = (double *) mem_free((char *) consmin);
  consmax  = (double *) mem_free((char *) consmax);
  sigma    = (double *) mem_free((char *) sigma);
  return(0);
}
