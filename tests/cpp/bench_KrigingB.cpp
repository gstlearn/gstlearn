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
#include "Basic/OptCustom.hpp"
#include "Neigh/NeighBench.hpp"
#include "Estimation/CalcKriging.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int argc, char *argv[])
{
  bool graphic = false;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("BenchKrigingB-");
  OptCustom::define("oldStyle", 0.);

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

  // Bench Neighborhood
  double width = 0.5;
  NeighBench* neighB = NeighBench::create(false, width);

  // Print the test environment
  message("This test is meant to test Kriging using Moving Neighborhood\n");
  message("- the Data Set contains %d samples\n", data->getSampleNumber(true));
  message("- the Output Grid contains %d nodes\n", grid->getSampleNumber(true));
  message("- the Bench Neighborhood is required:\n");
  message("  . Bench width = %lf\n", width);

  Timer timer;
  kriging(data, grid, model, neighB, EKrigOpt::POINT, true, false);
  timer.displayIntervalMilliseconds("Kriging in Bench Neighborhood", 1800);
  // Produce some stats for comparison
  DbStringFormat* dbfmt = DbStringFormat::create(FLAG_STATS, {"*estim"});
  grid->display(dbfmt);
  delete dbfmt;

  if (graphic)
    (void) grid->dumpToNF("Grid.ascii");
  delete neighB;
  delete data;
  delete grid;
  delete model;

  return (0);
}
