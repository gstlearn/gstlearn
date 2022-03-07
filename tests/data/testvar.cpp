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
#include "geoslib_f.h"

#include "Variogram/Vario.hpp"
#include "Basic/Law.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/File.hpp"
#include "Model/Model.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"

#include <iostream>
#include <fstream>

/*********************/
/* Program principal */
/*********************/

int main(int argc, char *argv[])
{
  char      *filename = new char[BUFFER_LENGTH];
  DbGrid    *dbout;
  Vario     *vario;
  Model     *model;
  Option_AutoFit mauto;
  Option_VarioFit options;
  Constraints constraints;
  int       nbsimu,nbtuba,seed,flag_norm_sill,flag_goulard_used;
  double    gof;
  static int verbose = 0;

  /* Initializations */

  dbout = (DbGrid *) NULL;
  vario = (Vario *) NULL;
  model = (Model *) NULL;
  flag_norm_sill = 0;
  flag_goulard_used = 1;
  double gofThresh = 2.;
  char string[1000];

  /* Standard output redirection to file */

  StdoutRedirect sr("Result.out");

  /* Setup the license */

  if (setup_license("Demonstration")) goto label_end;

  /* Setup constants */

  OptDbg::reset();

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
  dbout = DbGrid::createFromNF(filename,verbose);

  /* Look for simulations */

  ascii_filename("Simu",0,0,filename);
  ascii_simu_read(filename,verbose,&nbsimu,&nbtuba,&seed);

  /* Define the model */

  ascii_filename("Model",0,0,filename);
  model = Model::createFromNF(filename,verbose);
  if (model == (Model *) NULL) goto label_end;

  // Define and store the Space

  ASpaceObject::defineDefaultSpace(SPACE_RN,model->getDimensionNumber());

  /* Perform the non-conditional Simulation */
  
  if (dbout != (Db *) NULL)
  {
    if (simtub(nullptr,dbout,model,nullptr,nbsimu,seed,nbtuba,0))
      messageAbort("Simulations");
    /* Set the current variable to the conditional expectation */
    dbout->setLocatorByUID(dbout->getColumnNumber()-1,ELoc::Z);
  }
  seed = law_get_random_seed();
  
  /* Define the variogram */
  
  ascii_filename("Vario",0,0,filename);
  vario = Vario::createFromNF(filename,verbose);
  if (vario == (Vario *) NULL) goto label_end;
  if (dbout != (Db *) NULL)
  {
    vario->attachDb(dbout);
    vario->computeByKey("vg");
    ascii_filename("Vario",0,1,filename);
    if (vario->dumpToNF(filename,verbose))
      messageAbort("ascii_vario_write");
  }
  
  /* Fit the model */

  if (flag_norm_sill) opt_mauto_add_unit_constraints(mauto);
  options.setFlagGoulardUsed(flag_goulard_used);
  (void) model_auto_fit(vario,model,verbose,mauto,constraints,options);
  // Model is not printed any more to avoid differences among platforms
  //  model->display();
  ascii_filename("Model",0,1,filename);
  if (model->dumpToNF(filename,verbose))
    messageAbort("ascii_model_write");
  
  // produce the Goodness-of-fit score

  gof = model->gofToVario(vario);

  if (gof > gofThresh)
    message("Goodness-of-Fit (as a percentage of the Sill) may demonstrate an issue: %5.2lf\n",gof);
  else
    message("Goodness-of-Fit (<%5.2lf percent of the Sill): AutoFit is a success\n",gofThresh);

/* Core deallocation */

label_end:
  model = model_free(model);
  dbout = db_delete(dbout);
  vario = variogram_delete(vario);
  delete [] filename;
  return(0);
}
