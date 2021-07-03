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

  /* Initializations */

  new_model = (Model *) NULL;

  /* Modify the model */

  if (db->getGradientNumber() > 0)
  {
    /* Modify the gradients into standard variables */

    if (db_gradient_update(db)) return(new_model);

    /* Create the new Model */
    new_model = model_duplicate(model,model->getContext().getBallRadius(),1);
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
  Db        *dbin,*dbout;
  Vario     *vario;
  Model     *model,*new_model;
  Neigh     *neigh;
  Option_AutoFit mauto;
  int        nbsimu,seed,nbtuba;
  static int    nboot   = 10;
  static int    niter   = 10;
  static double toleps  = 1.;
  static int    verbose = 0;

  /* Initializations */

  dbin  = (Db    *) NULL;
  dbout = (Db    *) NULL;
  vario = (Vario *) NULL;
  model = new_model = (Model *) NULL;
  neigh = (Neigh *) NULL;

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

  /* Define the data */

  ascii_filename("Data",0,0,filename);
  dbin = ascii_db_read(filename,0,verbose);
  if (dbin == (Db *) NULL) goto label_end;
  db_print(dbin,1,0,1,1,1);

  /* Define the Default Space according to the Dimension of the Input Db */

  ASpaceObject::createGlobalSpace(SPACE_RN,dbin->getNDim());

  /* Define the output grid file */

  ascii_filename("Grid",0,0,filename);
  dbout = ascii_db_read(filename,1,verbose);

  /* Define the variogram */

  ascii_filename("Vario",0,0,filename);
  vario = ascii_vario_read(filename,verbose);
  if (vario != (Vario *) NULL)
  {
    if (variogram_compute(dbin,vario)) messageAbort("variogram_calcul");
    variogram_print(vario,1); 
    ascii_filename("Vario",0,1,filename);
    if (ascii_vario_write(filename,vario,verbose,1)) 
      messageAbort("ascii_vario_write");
  }

  /* Define the model */

  ascii_filename("Model",0,0,filename);
  model = ascii_model_read(filename,verbose);
  if (model == (Model *) NULL) goto label_end;
  if (vario != (Vario *) NULL) 
  {
    if (model_fitting_sills(vario,model,mauto)) goto label_end;
    ascii_filename("Model",0,1,filename);
    if (ascii_model_write(filename,model,verbose))
      messageAbort("ascii_model_write");
  }
  new_model = st_modify(model,dbin);
  if (new_model == (Model *) NULL) goto label_end;

  /* Define the neighborhood */

  ascii_filename("Neigh",0,0,filename);
  neigh = ascii_neigh_read(filename,verbose);

  /* Look for simulations */

  ascii_filename("Simu",0,0,filename);
  ascii_simu_read(filename,verbose,&nbsimu,&nbtuba,&seed);

  /* Conditional expectation */

  if (dbin->getIntervalNumber() > 0)
  {
    dbin->clearLocators(LOC_Z);
    if (gibbs_sampler(dbin,new_model,1,seed,nboot,niter,0,0,5.,toleps,1,1,1))
      messageAbort("krineq");
    /* Set the current variable to the conditional expectation */
    dbin->setLocatorByAttribute(dbin->getFieldNumber()-1,LOC_Z);
  }

  /* Perform the estimation */

  if (neigh != (Neigh *) NULL)
  {
    if (nbsimu > 0)
    {
    
      /* Simulation case */

      if (simtub(dbin,dbout,new_model,neigh,nbsimu,seed,nbtuba,0))
        messageAbort("Simulations");
      db_print(dbout,1,0,1,1,1);
    }
    else
    {
      if (neigh->getFlagXvalid())
      {
        
        /* Cross-validation */

        if (kriging(dbin,dbin,new_model,neigh,KOPTION_PONCTUAL,
                    1,1,0)) messageAbort("kriging");
        db_print(dbin,1,0,1,1,1);
      }
      else
      {

        /* Estimation case */

        if (dbout == (Db *) NULL) goto label_end;
        if (kriging(dbin,dbout,new_model,neigh,KOPTION_PONCTUAL,
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
  neigh = neigh_free(neigh);

  return(0);
}
