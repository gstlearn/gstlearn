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

#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Variogram/Vario.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Neigh/NeighBench.hpp"
#include "Neigh/NeighImage.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/File.hpp"
#include "Space/ASpaceObject.hpp"

/****************************************************************************
**
** FUNCTION: st_modify
**
** PURPOSE:  Duplicate and transform the model
**
** RETURNS:  The duplicated model
**
** IN_ARGS:  model: Input model
** IN_ARGS:  db   : Input Db structure
**
*****************************************************************************/
static Model *st_modify(Model *model,
                        Db    *db)
{
  Model *new_model;
  double ball_radius = 0.01;

  /* Initializations */

  new_model = nullptr;

  /* Modify the model */

  if (db->getGradientNumber() > 0)
  {
    /* Modify the gradients into standard variables */

    if (db_gradient_update(db)) return(new_model);

    /* Create the new Model */
    new_model = model_duplicate(model,ball_radius,1);
  }
  else
  {
    new_model = model_duplicate(model,0.,0);
  }
  return(new_model);
}

/*********************/
/* Program principal */
/*********************/

int main(int argc, char *argv[])

{
  char       filename[BUFFER_LENGTH];
  Db        *dbin;
  DbGrid    *dbout;
  Vario     *vario;
  Model     *model,*new_model;
  ANeighParam *neighparam;
  Option_AutoFit mauto;
  int        nbsimu,seed,nbtuba;
  static int    nboot   = 10;
  static int    niter   = 10;
  static int    verbose = 0;

  /* Initializations */

  dbin  = nullptr;
  dbout = nullptr;
  vario = nullptr;
  model = new_model = nullptr;
  neighparam = nullptr;

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

  /* Define the data */

  ascii_filename("Data",0,0,filename);
  dbin = Db::createFromNF(filename,verbose);
  if (dbin == nullptr) goto label_end;
  db_print(dbin,1,0,1,1,1);

  /* Define the Default Space according to the Dimension of the Input Db */

  ASpaceObject::defineDefaultSpace(SPACE_RN,dbin->getNDim());

  /* Define the output grid file */

  ascii_filename("Grid",0,0,filename);
  dbout = DbGrid::createFromNF(filename,verbose);

  /* Define the variogram */

  ascii_filename("Vario",0,0,filename);
  vario = Vario::createFromNF(filename,verbose);
  if (vario != nullptr)
  {
    vario->attachDb(dbin);
    vario->computeByKey("vg");
    vario->display();
    ascii_filename("Vario",0,1,filename);
    if (vario->dumpToNF(filename,verbose))
      messageAbort("ascii_vario_write");
  }

  /* Define the model */

  ascii_filename("Model",0,0,filename);
  model = Model::createFromNF(filename,verbose);
  if (model == nullptr) goto label_end;
  if (vario != nullptr) 
  {
    if (model_fitting_sills(vario,model,mauto)) goto label_end;
    ascii_filename("Model",0,1,filename);
    if (model->dumpToNF(filename,verbose))
      messageAbort("ascii_model_write");
  }
  new_model = st_modify(model,dbin);
  if (new_model == nullptr) goto label_end;

  /* Define the neighborhood */

  ascii_filename("Neigh",0,0,filename);
  neighparam = NeighUnique::createFromNF(filename,verbose);
  if (neighparam == nullptr)
    neighparam = NeighImage::createFromNF(filename, verbose);
  if (neighparam == nullptr)
    neighparam = NeighBench::createFromNF(filename, verbose);
  if (neighparam == nullptr)
    neighparam = NeighMoving::createFromNF(filename, verbose);

  /* Look for simulations */

  ascii_filename("Simu",0,0,filename);
  ascii_simu_read(filename,verbose,&nbsimu,&nbtuba,&seed);

  /* Conditional expectation */

  if (dbin->getIntervalNumber() > 0)
  {
    dbin->clearLocators(ELoc::Z);
    if (gibbs_sampler(dbin,new_model,nullptr,
                      1,seed,nboot,niter,false,true,false,false,0,
                      5.,true,true,true))
      messageAbort("gibbs_sampler");
    /* Set the current variable to the conditional expectation */
    dbin->setLocatorByUID(dbin->getColumnNumber()-1,ELoc::Z);
  }

  /* Perform the estimation */

  if (neighparam != nullptr)
  {
    if (nbsimu > 0)
    {
    
      /* Simulation case */

      if (simtub(dbin,dbout,new_model,neighparam,nbsimu,seed,nbtuba,0))
        messageAbort("Simulations");
      db_print(dbout,1,0,1,1,1);
    }
    else
    {
      if (neighparam->getFlagXvalid())
      {
        
        /* Cross-validation */

        if (xvalid(dbin,new_model,neighparam,
                   0,1,1,0)) messageAbort("xvalid");
        db_print(dbin,1,0,1,1,1);
      }
      else
      {

        /* Estimation case */

        if (dbout == nullptr) goto label_end;
        if (kriging(dbin,dbout,new_model,neighparam,EKrigOpt::PONCTUAL,
                    1,1,0)) messageAbort("kriging");
        db_print(dbout,1,0,1,1,1);
        dbout = db_delete(dbout);
      }
    }
  }

  /* Core deallocation */

label_end:
  dbin  = db_delete(dbin);
  dbout = db_delete(dbout);
  vario = variogram_delete(vario);
  model = model_free(model);
  new_model = model_free(new_model);
  delete neighparam;
  return(0);
}
