/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Basic/Law.hpp"
#include "Enum/ECov.hpp"
#include "Enum/ELoadBy.hpp"

#include "API/SPDE.hpp"

#include "Basic/OptDbg.hpp"
#include "Basic/File.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Space/ASpaceObject.hpp"
#include "Model/Model.hpp"

#define VERBOSE 0

/****************************************************************************/
/*!
** Main Program
**
*****************************************************************************/
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  /***********************/
  /* 1 - Initializations */
  /***********************/

  /* 1.b - Setup the default space */

  int ndim     = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);
  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("SPDE-");

  /* 1.c - Setup constants */

  OptDbg::reset();
  
  // Create the 2-D grid output file

  VectorInt    nx = { 400, 300 };
  VectorDouble dx = { 1., 1. };
  VectorDouble x0 = { 0., 0. };
  DbGrid *dbgrid = DbGrid::create(nx, dx, x0, VectorDouble(), ELoadBy::COLUMN,
                                  VectorDouble(), VectorString(),
                                  VectorString(), 1);
    
  // Model 

  double range    = 79.8;
  double sill     = 1.;
  double param    = 1.;
  Model* model = Model::createFromParam(ECov::MATERN,range,sill,param);

  // Perform the non-conditional simulation

  bool verbose    = false;
  int seed        = 31415;
  int nsimu       = 10;
  int useCholesky = 1;
  law_set_random_seed(seed);
  (void) simulateSPDE(NULL, dbgrid, model, nullptr, nsimu, NULL, useCholesky,
                      SPDEParam(), verbose);

  // Print statistics on the results

  DbStringFormat dbfmt;
  dbfmt.setFlags(true, true, true, true, true);
  dbgrid->display(&dbfmt);
  (void) dbgrid->dumpToNF("pgs.ascii");
  
  delete dbgrid;
  delete model;
  return(0);
}
