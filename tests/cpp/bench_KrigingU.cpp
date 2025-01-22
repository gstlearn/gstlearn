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
#include "Basic/OptDbg.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Estimation/CalcKriging.hpp"

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

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("BenchKrigingU-");

  // Global parameters
  int ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);
  OptCustom::define("oldStyle", 1.);

  // Generate the data base
  int nech = 100;
  int nvar = 1;
  bool verbose = false;
  Db* data = Db::createFillRandom(nech, ndim, nvar);

  // Generate the output grid
  int ncell = 100;
  VectorInt nx = {ncell, ncell};
  VectorDouble dx = {1. / ncell, 1. / ncell};
  DbGrid* grid = DbGrid::create(nx, dx);

  // Create the Model
  double range = 1. / 5.;
  double sill = 2.;
  Model* model = Model::createFromParam(ECov::SPHERICAL, range, sill);

  // Unique Neighborhood
  NeighUnique* neighU = NeighUnique::create();

  // Set the verbose option
  if (verbose) OptDbg::setReference(1);

  // Print the test environment
  message("This test is meant to test Kriging using Unique Neighborhood\n");
  message("- the Data Set contains %d samples\n", data->getSampleNumber(true));
  message("- the output Grid contains %d nodes\n", grid->getSampleNumber(true));
  message("- the Unique Neighborhood is required\n");

  Timer timer;
  kriging(data, grid, model, neighU, EKrigOpt::POINT, true, false, false);
  timer.displayIntervalMilliseconds("Kriging in Unique Neighborhood", 700);

  // Produce some statistics for comparison
  DbStringFormat* gridfmt = DbStringFormat::create(FLAG_STATS, {"*estim"});
  grid->display(gridfmt);
  delete gridfmt;

  delete neighU;
  delete data;
  delete grid;
  delete model;

  return (0);
}
