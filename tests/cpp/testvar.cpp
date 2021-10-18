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
#include "geoslib_d.h"
#include "geoslib_f.h"
#include "geoslib_old_f.h"

#include <stdlib.h>

/*********************/
/* Program principal */
/*********************/

int main(int argc, char *argv[])
{
  char      *filename = new char[BUFFER_LENGTH];
  Db        *dbout;
  Vario     *vario;
  Model     *model;
  Option_AutoFit mauto;
  Option_VarioFit options;
  Constraints constraints;
  int       nbsimu,nbtuba,seed,flag_norm_sill,flag_goulard_used;
  static int verbose = 0;

  /* Initializations */

  dbout = (Db    *) NULL;
  vario = (Vario *) NULL;
  model = (Model *) NULL;
  flag_norm_sill = 0;
  flag_goulard_used = 1;

  /* Connect the Geoslib Library */

  if (setup_license("Demonstration")) goto label_end;

  /* Setup constants */

  debug_reset();
  constant_reset();

  /* Getting the Study name */

  if (argc != 2) messageAbort("Wrong number of arguments");
  ascii_study_define(argv[1]);

  /* Define the environment */

  ascii_filename("Environ",0,0,filename);
  ascii_environ_read(filename,verbose);

  /* Define the options */

  ascii_filename("Option",0,0,filename);
  ascii_option_defined(filename,0,"Norm_sill",0,&flag_norm_sill);
  ascii_option_defined(filename,0,"Goulard_used",0,&flag_goulard_used);

  /* Define the output grid file */

  ascii_filename("Grid",0,0,filename);
  dbout = ascii_db_read(filename,1,verbose);

  /* Look for simulations */

  ascii_filename("Simu",0,0,filename);
  ascii_simu_read(filename,verbose,&nbsimu,&nbtuba,&seed);

  /* Define the model */

  ascii_filename("Model",0,0,filename);
  model = ascii_model_read(filename,verbose);
  if (model == (Model *) NULL) goto label_end;
  
  /* Perform the non-conditional Simulation */
  
  if (dbout != (Db *) NULL)
  {
    if (simtub((Db *) NULL,dbout,model,(Neigh *) NULL,nbsimu,seed,nbtuba,0))
      messageAbort("Simulations");
    /* Set the current variable to the conditional expectation */
    dbout->setLocatorByAttribute(dbout->getFieldNumber()-1,ELoc::Z);
  }
  seed = law_get_random_seed();
  
  /* Define the variogram */
  
  ascii_filename("Vario",0,0,filename);
  vario = ascii_vario_read(filename,verbose);
  if (vario == (Vario *) NULL) goto label_end;
  if (dbout != (Db *) NULL)
  {
    vario->attachDb(dbout);
    vario->compute("vg");
    ascii_filename("Vario",0,1,filename);
    if (vario->serialize(filename,verbose))
      messageAbort("ascii_vario_write");
  }
  
  /* Fit the model */

  if (flag_norm_sill) opt_mauto_add_unit_constraints(mauto);
  options.setFlagGoulardUsed(flag_goulard_used);
  (void) model_auto_fit(vario,model,verbose,mauto,constraints,options);
  model->display();
  ascii_filename("Model",0,1,filename);
  if (model->serialize(filename,verbose))
    messageAbort("ascii_model_write");
  
/* Core deallocation */

label_end:
  model = model_free(model);
  dbout = db_delete(dbout);
  vario = variogram_delete(vario);
  delete [] filename;
  return(0);
}
