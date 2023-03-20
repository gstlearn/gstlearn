/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "geoslib_d.h"
#include "geoslib_f.h"
#include "geoslib_old_f.h"

#include "Enum/ECov.hpp"
#include "Enum/ELoadBy.hpp"

#include "Basic/OptDbg.hpp"
#include "Basic/File.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
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
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());
  DbStringFormat dbfmt;

  DbGrid      *dbgrid;
  Model       *model = nullptr;
  SPDE_Option  s_option;
  CovContext   ctxt;
  const char triswitch[] = "nqQ";
  int verbose, seed, ndim, nsimu;
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

  /* 1.b - Setup the default space */

  defineDefaultSpace(ESpaceType::RN, ndim);

  /* 1.c - Setup constants */

  OptDbg::reset();
  
  // Create the 2-D grid output file

  dbgrid = DbGrid::create(nx, dx, x0, VectorDouble(), ELoadBy::COLUMN,
                          VectorDouble(), VectorString(), VectorString(), 1);
  if (dbgrid == nullptr) goto label_end;
  if (db_extension_diag(dbgrid,&diag)) goto label_end;
    
  // Model 

  ctxt = CovContext(1, ndim, diag);
  model = new Model(ctxt);
  if (model == nullptr) goto label_end;
  if (model_add_cova(model,ECov::BESSEL_K,0,0,range,param,
                     VectorDouble(),VectorDouble(),sill,0.)) goto label_end;

  // Perform the non-conditional simulation

  s_option = spde_option_alloc();
  spde_option_update(s_option,triswitch);
  if (spde_check(NULL,dbgrid,model,NULL,verbose,gext,
                 true, true, true, false, false, false, false)) goto label_end;
  if (spde_f(NULL,dbgrid,model,gext,s_option,1,1,
             seed,nsimu,0,0,0,0,0,0,0,verbose))
    goto label_end;
  
  // Print statistics on the results

  dbfmt.setFlags(true, true, true, true, true);
  dbgrid->display(&dbfmt);
  
label_end:
  delete dbgrid;
  delete model;
  return(0);
}
