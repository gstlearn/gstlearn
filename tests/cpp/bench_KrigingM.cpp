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
#include "Neigh/NeighMoving.hpp"
#include "Estimation/CalcKriging.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int argc, char *argv[])
{
  bool graphic = true;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("BenchKrigingM-");

  // Global parameters
  int ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Generate the data base
  int nech = 1000;
  int nvar = 1;
  Db* data = Db::createFillRandom(nech, ndim, nvar);

  // Generate the output grid
  int ncell       = 100;
  VectorInt nx    = {ncell, ncell};
  VectorDouble dx = {1. / ncell, 1. / ncell};
  DbGrid* grid    = DbGrid::create(nx, dx);

  // Create the Model
  double range = 1. / 5.;
  double sill  = 2.;
  Model* model = Model::createFromParam(ECov::SPHERICAL, range, sill);

  // Moving Neighborhood
  int nmaxi = 20;
  int nmini = 2;
  int nsect = 8;
  int nsmax = 3;
  double radius = 1.;

  NeighMoving* neighM =
    NeighMoving::create(false, nmaxi, radius, nmini, nsect, nsmax);
  int leaf_size = 30;
  neighM->setBallSearch(true, leaf_size);

  // Print the test environment
  message("This test is mean to test Kriging using Moving Neighborhood\n");
  message("- the Data Set contains %d samples\n", data->getSampleNumber(true));
  message("- the Output Grid contains %d nodes\n", grid->getSampleNumber(true));
  message("- the Moving Neighborhood is required:\n");
  message("  . Radius dimension = %lf\n", radius);
  message("  . Maximum number of neighbors = %d\n", nmaxi);
  message("  . Minimum number of neighbors = %d\n", nmini);
  message("  . Number of angular sectors   = %d\n", nsect);
  message("  . Maxmimum number of neighbors per sector = %d\n", nsmax);
  message("  . Leaf Size for Ball Tree algorithm = %d\n", leaf_size);

  Timer timer;
  kriging(data, grid, model, neighM, EKrigOpt::POINT, true, false);
  timer.displayIntervalMilliseconds("Kriging in Moving Neighborhood", 1500);

  // Produce some stats for comparison
  DbStringFormat* dbfmt = DbStringFormat::create(FLAG_STATS, {"*estim"});
  grid->display(dbfmt);
  delete dbfmt;

  if (graphic)
    (void) grid->dumpToNF("Grid.ascii");

  delete neighM;
  delete data;
  delete grid;
  delete model;

  return (0);
}
