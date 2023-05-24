/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_d.h"
#include "geoslib_f.h"
#include "geoslib_old_f.h"

#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Variogram/Vario.hpp"
#include "Neigh/ANeigh.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Neigh/NeighBench.hpp"
#include "Neigh/NeighImage.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/File.hpp"
#include "Space/ASpaceObject.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Estimation/CalcKriging.hpp"

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

  if (db->getLocNumber(ELoc::G) > 0)
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
  ANeigh    *neigh;
  Option_AutoFit mauto;
  Constraints constraints;
  DbStringFormat dbfmt;
  int        nbsimu,seed,nbtuba;
  static int    nboot   = 10;
  static int    niter   = 10;
  static int    verbose = 0;

  /* Initializations */

  dbin  = nullptr;
  dbout = nullptr;
  vario = nullptr;
  model = new_model = nullptr;
  neigh = nullptr;

  /* Standard output redirection to file */

  StdoutRedirect sr("Result.out");

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
  dbfmt.setFlags(true, false, true, true, true);
  dbin->display(&dbfmt);

  /* Define the Default Space according to the Dimension of the Input Db */

  defineDefaultSpace(ESpaceType::RN,dbin->getNDim());

  /* Define the output grid file */

  ascii_filename("Grid",0,0,filename);
  dbout = DbGrid::createFromNF(filename,verbose);
//  if (dbout != nullptr) dbout->display(&dbfmt);

  /* Define the variogram */

  ascii_filename("Vario",0,0,filename);
  vario = Vario::createFromNF(filename,verbose);
  if (vario != nullptr)
  {
    vario->attachDb(dbin);
    vario->computeByKey("vg");
    vario->display();
    ascii_filename("Vario",0,1,filename);
    if (! vario->dumpToNF(filename,verbose))
      messageAbort("ascii_vario_write");
  }

  /* Define the model */

  ascii_filename("Model",0,0,filename);
  model = Model::createFromNF(filename,verbose);
  if (model == nullptr) goto label_end;
  if (vario != nullptr) 
  {
    if (model_fitting_sills(vario,model,constraints,mauto)) goto label_end;
    ascii_filename("Model",0,1,filename);
    if (! model->dumpToNF(filename,verbose))
      messageAbort("ascii_model_write");
  }
  new_model = st_modify(model,dbin);
  if (new_model == nullptr) goto label_end;

  /* Define the neighborhood */

  ascii_filename("Neigh",0,0,filename);
  neigh = NeighUnique::createFromNF(filename,verbose);
  if (neigh == nullptr)
    neigh = NeighImage::createFromNF(filename, verbose);
  if (neigh == nullptr)
    neigh = NeighBench::createFromNF(filename, verbose);
  if (neigh == nullptr)
    neigh = NeighMoving::createFromNF(filename, verbose);

  /* Look for simulations */

  ascii_filename("Simu",0,0,filename);
  ascii_simu_read(filename,verbose,&nbsimu,&nbtuba,&seed);

  /* Conditional expectation */

  if (dbin->getIntervalNumber() > 0)
  {
    dbin->clearLocators(ELoc::Z);
    if (gibbs_sampler(dbin,new_model,
                      1,seed,nboot,niter,false,false,true,false,false,0,
                      5.,true,true,true))
      messageAbort("gibbs_sampler");
    /* Set the current variable to the conditional expectation */
    dbin->setLocatorByUID(dbin->getColumnNumber()-1,ELoc::Z);
  }

  /* Perform the estimation */

  if (neigh != nullptr)
  {
    if (nbsimu > 0)
    {
    
      /* Simulation case */

      if (simtub(dbin,dbout,new_model,neigh,nbsimu,seed,nbtuba,0))
        messageAbort("Simulations");
      dbfmt.setFlags(true, false, true, true, true);
      dbout->display(&dbfmt);
    }
    else
    {
      if (dbout == nullptr)
      {
        
        /* Cross-validation */

        if (xvalid(dbin,new_model,neigh,0,1,0,0)) messageAbort("xvalid");
        dbfmt.setFlags(true, false, true, true, true);
        dbin->display(&dbfmt);
      }
      else
      {

        /* Estimation case */

        if (kriging(dbin,dbout,new_model,neigh,EKrigOpt::POINT,
                    1,1,0)) messageAbort("kriging");
        dbfmt.setFlags(true, false, true, true, true);
        dbout->display(&dbfmt);
        delete dbout;
        dbout = nullptr;
      }
    }
  }

  /* Core deallocation */

label_end:
  delete dbin;
  delete dbout;
  delete model;
  delete new_model;
  delete neigh;
  delete vario;
  return(0);
}
