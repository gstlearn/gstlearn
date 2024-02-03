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
#include "Tree/Ball.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int argc, char *argv[])
{
  Timer timer;
  VectorDouble vec;
  VectorDouble vecb;
  VectorDouble diff;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Migrate-");

  // Global parameters
  int ndim = 2;
//  VectorDouble dmax = {0.01, 0.01};
  VectorDouble dmax;
  bool flag_stat = true;
  bool flag_simu = false;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Generate the Model and the Neighborhood (for simulation)
  double range = 0.3;
  Model* model = Model::createFromParam(ECov::CUBIC,range);
  NeighUnique* neigh = NeighUnique::create();

  // Generate the data set
  int nech = 10000;
  Db* data = Db::createFillRandom(nech, ndim, 0);
  if (! flag_simu)
  {
    VectorDouble resd = VH::add(data->getColumn("x-1"), data->getColumn("x-2"));
    data->addColumns(resd, "Simu", ELoc::Z);
  }
  else
  {
    (void) simtub(nullptr, data, model, neigh);
  }
  if (flag_stat) data->display();

  // Generate the output grid
  int nx = 200;
  DbGrid* grid = DbGrid::createCoveringDb(data, {nx,nx});
  if (! flag_simu)
  {
    VectorDouble resg = VH::add(grid->getColumn("x1"), grid->getColumn("x2"));
    grid->addColumns(resg, "Simu", ELoc::Z);
  }
  else
  {
    (void) simtub(nullptr, grid, model, neigh);
  }
  if (flag_stat) grid->display();

  // Testing several migrations
  mestitle(1, "Migrations");
  message("- with a Point data base with %d points\n",nech);
  message("- with a Grid  data base with %d nodes\n",nx*nx);
  if (! dmax.empty()) VH::display("- using Dmax criterion", dmax);
  message("\n");

  // Migrate Grid -> Point (no interpolation)
  timer.reset();
  (void) migrate(grid, data, "Simu", 1, dmax, false, false, false,
                 NamingConvention("G2Pn",false, false));
  timer.displayIntervalMilliseconds("Migrate Grid to Point (no interpolation)", 5);

  // Migrate Grid -> Point (with interpolation)
  timer.reset();
  (void) migrate(grid, data, "Simu", 1, dmax, false, true, false,
                 NamingConvention("G2Py",false, false));
  timer.displayIntervalMilliseconds("Migrate Grid to Point (with interpolation)", 5);

  // Migrate Grid -> Grid
  timer.reset();
  (void) migrate(grid, grid, "Simu", 1, dmax, false, false, false,
                 NamingConvention("G2Gn",false, false));
  timer.displayIntervalMilliseconds("Migrate Grid to Grid", 30);

  // Expand Grid -> Grid
  timer.reset();
  (void) migrate(grid, grid, "Simu", 1, dmax, true, false, false,
                 NamingConvention("G2Gy",false, false));
  timer.displayIntervalMilliseconds("Expand Grid to Grid", 30);

  // Migrate Point -> Grid
  timer.reset();
  (void) migrate(data, grid, "Simu", 1, dmax, false, false, false,
                 NamingConvention("P2Gn",false,false));
  timer.displayIntervalMilliseconds("Migrate Point to Grid", 1500);

  // Expand Point -> Grid
  timer.reset();
  (void) migrate(data, grid, "Simu", 1, dmax, true, false, false,
                 NamingConvention("P2Gyr",false,false));
  timer.displayIntervalMilliseconds("Expand Point to Grid", 1500);
  vec = grid->getColumn("P2Gyr*");

  // Expand Point -> Grid (with balltree)
  timer.reset();
  (void) migrate(data, grid, "Simu", 1, dmax, true, false, true,
  NamingConvention("P2Gyb",false,false));
  timer.displayIntervalMilliseconds("Expand Point to Grid (with balltree)", 1500);
  vecb = grid->getColumn("P2Gyb*");

  // Compare impact of balltree option on P2Gy
  diff = VH::subtract(vec, vecb);
  VH::getMostSignificant(diff, EPSILON6, 10);
  grid->addColumns(diff, "Diff_P2Gy");

  // Expand Point -> Point
  timer.reset();
  (void) migrate(data, data, "Simu", 1, dmax, true, false, false,
                 NamingConvention("P2Pyr",false,false));
  timer.displayIntervalMilliseconds("Expand Point to Point", 19000);
  vec = data->getColumn("P2Pyr*");

  // Expand Point -> Point (with balltree)
  timer.reset();
  (void) migrate(data, data, "Simu", 1, dmax, true, false, true,
                 NamingConvention("P2Pyb",false,false));
  timer.displayIntervalMilliseconds("Expand Point to Point (with balltree)", 19000);
  vecb = data->getColumn("P2Pyb*");

  // Compare impact of balltree option on P2Py
  diff = VH::subtract(vec, vecb);
  VH::getMostSignificant(diff, EPSILON6, 10);
  data->addColumns(diff, "Diff_P2Py");

  // Produce some neutral file
  (void) data->dumpToNF("data.ascii");
  (void) grid->dumpToNF("grid.ascii");

  // Printout some stats
  if (flag_stat)
  {
    data->display();
    grid->display();
  }

  if (data     != nullptr) delete data;
  if (grid     != nullptr) delete grid;
  if (model    != nullptr) delete model;
  if (neigh    != nullptr) delete neigh;

  return (0);
}
