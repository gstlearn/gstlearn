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
#include "geoslib_d.h"
#include "geoslib_f.h"

#include "Enum/ESpaceType.hpp"
#include "Enum/ECov.hpp"
#include "Enum/EKrigOpt.hpp"

#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Basic/File.hpp"
#include "Basic/Timer.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int argc, char *argv[])
{
  bool verbose = false;
  Timer timer;
  DbStringFormat* dbfmt;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Migrate-");

  // Global parameters
  int ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Generate the data set
  int nech = 10000;
  int nvar = 0;
  Db* data = Db::createFillRandom(nech, ndim, nvar);
  if (verbose) data->display();

  // Generate the output grid
  int nx = 200;
  VectorInt nodes = {nx,nx};
  DbGrid* grid = DbGrid::createCoveringDb(data, nodes);
  if (verbose) grid->display();

  // Generate a structured variable on Grid
  double range = 0.3;
  Model* model = Model::createFromParam(ECov::CUBIC,range);
  NeighUnique* neigh = NeighUnique::create();
  (void) simtub(nullptr, grid, model, neigh);
  (void) grid->dumpToNF("Init.ascii");

  // Testing several migrations
  mestitle(1, "Migrations");
  message("- with a Point data base with %d points\n",nech);
  message("- with a Grid data base with %d nodes\n",nx*nx);
  message("Note: The 'dmax' argument has been added to limitate extrapolation\n");
  message("      This option creates some empty gaps on the results on Grid\n");

  // Migrate the information Grid -> Point
  timer.reset();
  (void) migrate(grid, data, "Simu", 1, VectorDouble(), false, false,
                 NamingConvention("G2P",false, false));
  timer.displayIntervalMilliseconds("Migrate Grid to Point", 5);
  (void) data->dumpToNF("G2P.ascii");

  // Produce some statistics for comparison
  dbfmt = DbStringFormat::create(FLAG_STATS, {"G2P"});
  data->display(dbfmt);
  delete dbfmt;

  // Migrate the information Point -> Grid
  timer.reset();
  VectorDouble dmax = {0.01, 0.01};
  (void) migrateByLocator(data, grid, ELoc::Z, 1, dmax, true, false,
                          NamingConvention("P2G",false,false));
  timer.displayIntervalMilliseconds("Migrate Point to Grid", 2000);
  (void) grid->dumpToNF("P2G.ascii");

  // Produce some statistics for comparison
  dbfmt = DbStringFormat::create(FLAG_STATS, {"P2G"});
  grid->display(dbfmt);
  delete dbfmt;

  if (data  != nullptr) delete data;
  if (grid  != nullptr) delete grid;
  if (model != nullptr) delete model;
  if (neigh != nullptr) delete neigh;

  return (0);
}
