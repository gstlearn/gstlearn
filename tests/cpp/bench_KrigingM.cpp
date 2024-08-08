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
  bool verbose = false;
  bool graphic = true;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("BenchKrigingM-");

  // Global parameters
  defineDefaultSpace(ESpaceType::RN, 2);

  // Generate the data base
  String filename = ASerializable::getTestData("benchmark","sic_obs.dat");
  CSVformat csv(false, 6);
  Db* data = Db::createFromCSV(filename, csv);
  data->setName("New.1","ID");
  data->setName("New.2","X");
  data->setName("New.3","Y");
  data->setName("New.4","rainfall");
  data->setLocators({"X","Y"},ELoc::X);
  data->setLocator("rainfall",ELoc::Z);
  if (verbose) data->display();

  // Generate the output grid
  bool flagSmall = false;
  VectorInt nx;
  if (flagSmall)
    nx = {50,60};
  else
    nx = {360,240};
  VectorDouble dx = {1000, 1000};
  VectorDouble x0 = {-180000, -120000};
  DbGrid* grid = DbGrid::create(nx, dx, x0);
  if (verbose) grid->display();
  if (graphic) (void) data->dumpToNF("Data.ascii");

  // Create the Model
  Model* model = Model::createFromParam(ECov::SPHERICAL, 80000, 14000);
  if (verbose) model->display();

  // Moving Neighborhood
  int nmaxi = 20;
  int nmini = 2;
  int nsect = 8;
  int nsmax = 3;
  double radius = 50000;

  NeighMoving* neighM = NeighMoving::create(false, nmaxi, radius, nmini, nsect, nsmax);
  if (verbose) neighM->display();

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
