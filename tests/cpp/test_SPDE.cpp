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

#include "Basic/OptDbg.hpp"
#include "Db/DbGrid.hpp"
#include "Db/ELoadBy.hpp"
#include "Space/ASpaceObject.hpp"
#include "Model/Model.hpp"
#include "Covariances/CovContext.hpp"

#define VERBOSE 0

/****************************************************************************/
/*!
** Main Program
**
*****************************************************************************/
int main(int /*argc*/, char */*argv*/[])

{
  DbGrid      *dbgrid;
  Model       *model = nullptr;
  SPDE_Option  s_option;
  CovContext   ctxt;
  const char triswitch[] = "nqQ";
  int verbose, seed, ndim, iptr, nsimu;
  double diag,range,param;
  VectorInt nx = { 400, 300 };
  VectorDouble dx = { 1., 1. };
  VectorDouble x0 = { 0., 0. };
  
  /***********************/
  /* 1 - Initializations */
  /***********************/

  dbgrid   = nullptr;
  model    = nullptr;
  verbose  = 0;
  seed     = 31415;
  range    = 79.8;
  param    = 1.;
  ndim     = 2;
  nsimu    = 10;
  VectorDouble sill{1.};
  VectorDouble gext{ 2*79.8, 2*79.8 };

  /* 1.b - Setup the license */

  if (setup_license("Demonstration")) goto label_end;
  ASpaceObject::defineDefaultSpace(SPACE_RN, ndim);

  /* 1.c - Setup constants */

  OptDbg::reset();
  
  // Create the 2-D grid output file

  dbgrid = db_create_grid(0,ndim,0,ELoadBy::COLUMN,1,nx,x0,dx);
  if (dbgrid == nullptr) goto label_end;
  if (db_locator_attribute_add(dbgrid,ELoc::X,ndim,0,0.,
                               &iptr)) goto label_end;
  if (db_grid_define_coordinates(dbgrid)) goto label_end;
  if (db_extension_diag(dbgrid,&diag)) goto label_end;
    
  // Model 

  ctxt = CovContext(1, ndim, 1000, diag);
  model = new Model(ctxt);
  if (model == nullptr) goto label_end;
  if (model_add_cova(model,ECov::BESSEL_K,0,0,range,param,
                     VectorDouble(),VectorDouble(),sill,0.)) goto label_end;

  // Perform the non-conditional simulation

  s_option = spde_option_alloc();
  spde_option_update(s_option,triswitch);
  if (spde_f(NULL,dbgrid,model,gext,s_option,1,1,
           seed,nsimu,0,0,0,0,0,0,0,verbose))
    goto label_end;
  
  // Print statistics on the results

  db_print(dbgrid,1,1,1,1);
  
label_end:
  dbgrid   = db_delete(dbgrid);
  model    = model_free(model);
  return(0);
}
