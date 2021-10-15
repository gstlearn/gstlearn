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
#include "Basic/Limits.hpp"
#include "geoslib_d.h"
#include "geoslib_f.h"
#include <stdlib.h>

/*********************/
/* Program principal */
/*********************/

int main(int argc, char *argv[])

{
  char   filename[BUFFER_LENGTH];
  Db        *dbin,*dbout;
  Vario     *vario;
  Model     *model[2][2];
  Neigh     *neigh;
  Rule      *rule[2];
  Option_VarioFit options;
  double  delta,total;
  int     i,j,lec,nbsimu,seed,nbtuba,npgs,ntot,nfac[2];
  int     flag_vario,flag_grid,iatt_z,iatt_ind,ifac,nclass;
  VectorDouble props;
  RuleProp ruleprop;
  static int    niter   = 100;
  static int    nboot   = 10;
  static double toleps  = 1.;
  static int    verbose = 1;

  /* Initializations */

  npgs    = ntot = 0;
  dbin    = (Db     *) NULL;
  dbout   = (Db     *) NULL;
  vario   = (Vario *) NULL;
  neigh   = (Neigh  *) NULL;
  for (i=0; i<2; i++)
  {
    rule[i] = (Rule *) NULL;
    for (j=0; j<2; j++)
      model[i][j] = (Model *) NULL;
  }

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
  iatt_z = db_attribute_identify(dbin,ELoc::Z,0);
  if (verbose) db_print(dbin,1,0,1,1,1);

  /* Define the Default Space according to the Dimension of the Input Db */

  ASpaceObject::defineDefaultSpace(SPACE_RN,dbin->getNDim());

  /* Define the variogram (optional) */
  
  ascii_filename("Vario",0,0,filename);
  vario = ascii_vario_read(filename,verbose);
  flag_vario = (vario != (Vario *) NULL);

  /* Define the output grid file */

  ascii_filename("Grid",0,0,filename);
  dbout = ascii_db_read(filename,1,verbose);
  flag_grid = (dbout != (Db *) NULL);

  /* Define the rules */

  for (i=lec=0; i<2; i++)
  {

    /* Read the rule */

    ascii_filename("Rule",i,0,filename);
    rule[i] = ascii_rule_read(filename,verbose);
    if (rule[i] == (Rule *) NULL) continue;

    npgs++;
    if (verbose) rule[i]->display(false, false);
    nfac[i] = rule[i]->getFaciesNumber();

    /* Define the models */
    
    for (j=0; j<2; j++,lec++)
    {
      if (! rule[i]->isYUsed(j)) continue;
      ascii_filename("Model",lec,0,filename);
      model[i][j] = ascii_model_read(filename,verbose);
      if (model[i][j] == (Model *) NULL) goto label_end;
    }

    /* Calculate the experimental variogram of indicators */
    
    if (flag_vario)
    {

      /* Define the indicators */

      nclass   = nfac[i];
      iatt_ind = dbin->getFieldNumber();
      Limits limits = Limits(nclass);
      limits.toIndicator(dbin);
      dbin->setLocatorsByAttribute(nclass,iatt_ind,ELoc::Z);
      
      /* Calculate the experimental variograms */
      
      vario->attachDb(dbin);
      vario->compute("vg");
      if (verbose) variogram_print(vario,1);
      ascii_filename("Vario",0,1,filename);
      if (vario->serialize(filename,verbose))
        messageAbort("ascii_vario_write");
      
      /* Delete the indicator variables */
      
      dbin->clearLocators(ELoc::Z);
      for (ifac=0; ifac<nclass; ifac++)
        dbin->deleteFieldByAttribute(iatt_ind+ifac);
      dbin->setLocatorByAttribute(iatt_z,ELoc::Z);
    }
  }

  /* Look for simulations */

  ascii_filename("Simu",0,0,filename);
  ascii_simu_read(filename,verbose,&nbsimu,&nbtuba,&seed);

  /* Create the proportions */

  ntot  = (npgs == 1) ? nfac[0] : nfac[0] * nfac[1];
  props.resize(ntot);
  total = 0.;
  for (i=0; i<ntot; i++)
  {
    props[i] = law_uniform(0.,1.);
    total   += props[i];
  }
  for (i=0; i<ntot; i++) props[i] /= total;

  /* Define the neighborhood */

  neigh = neigh_init_unique(dbout->getNDim());

  /* Perform the Pluri-Gaussian Simulations */

  if (flag_grid)
  {
    delta = dbout->getDX(0);
    if (npgs == 1)
    {
      ruleprop = RuleProp(rule[0],props);
      if (simpgs(dbin,dbout,&ruleprop,model[0][0],model[0][1],
                 neigh,nbsimu,seed,0,0,0,0,nbtuba,nboot,niter,1,toleps,
                 delta)) goto label_end;
    }
    else
    {
      ruleprop = RuleProp(rule[0],rule[1],props);
      if (simbipgs(dbin,dbout,&ruleprop,
                   model[0][0],model[0][1],model[1][0],model[1][1],
                   neigh,nbsimu,seed,0,0,0,0,nbtuba,nboot,niter,1,toleps))
        goto label_end;
    }
    if (verbose) db_print(dbout,1,0,1,1,1);
  }

  /* Serialization of results (for visual check) */

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("testPGS-");
  dbout->serialize("Result");

  /* Core deallocation */

label_end:
  dbin   = db_delete(dbin);
  dbout  = db_delete(dbout);
  for (i=0; i<2; i++)
  {
    rule[i] = rule_free(rule[i]);
    for (j=0; j<2; j++)
      model[i][j] = model_free(model[i][j]);
  }
  neigh = neigh_free(neigh);
  return(0);
}
