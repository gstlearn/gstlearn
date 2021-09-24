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
#include <stdlib.h>

/*********************/
/* Program principal */
/*********************/

int main(int argc, char *argv[])

{
  char   filename[BUFFER_LENGTH];
  Db     *db;
  Vario  *vario;
  Model  *model;
  Neigh  *neigh;
  Rule   *rule;
  int     rank;
  static int verbose = 0;

  /* Initializations */

  db    = (Db    *) NULL;
  vario = (Vario *) NULL;
  model = (Model *) NULL;
  neigh = (Neigh *) NULL;
  rule  = (Rule  *) NULL;

  /* Connect the Geoslib Library */

  if (setup_license("Demonstration")) return(0);
  ASerializable::setSerializedContainerName("");

  /* Setup constants */

  debug_reset();
  constant_reset();
  if (argc != 2) messageAbort("Wrong number of arguments.\n"
                              "Please provide the data directory");
  ascii_study_define(argv[1]);

  /* Loop on the Rank of different files */

  for (rank=0; rank<4; rank++)
  {

    /* Update the files */
    
    ascii_filename("Data",rank,0,filename);
    db = ascii_db_read(filename,0,verbose);
    if (db != (Db *) NULL) 
    {
      db->serialize(filename,verbose);
      db  = db_delete(db);
    }

    ascii_filename("Grid",rank,0,filename);
    db = ascii_db_read(filename,1,verbose);
    if (db != (Db *) NULL) 
    {
      db->serialize(filename,verbose);
      db  = db_delete(db);
    }

    ascii_filename("Vario",rank,0,filename);
    vario = ascii_vario_read(filename,verbose);
    if (vario != (Vario *) NULL) 
    {
      vario->serialize(filename,verbose);
      vario  = variogram_delete(vario);
    }

    ascii_filename("Model",rank,0,filename);
    model = ascii_model_read(filename,verbose);
    if (model != (Model *) NULL) 
    {
      model->serialize(filename,verbose);
      model = model_free(model);
    }

    ascii_filename("Neigh",rank,0,filename);
    neigh = ascii_neigh_read(filename,verbose);
    if (neigh != (Neigh *) NULL) 
    {
      neigh->serialize(filename,verbose);
      neigh = neigh_free(neigh);
    }

    ascii_filename("Rule",rank,0,filename);
    rule = ascii_rule_read(filename,verbose);
    if (rule != (Rule *) NULL) 
    {
      rule->serialize(filename,verbose);
      rule = rule_free(rule);
    }
  }

  return(0);
}
