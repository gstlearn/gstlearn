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

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int argc, char *argv[])
{
  bool verbose  = false;
  bool graphic  = false;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("BenchKrigingU3D-");

  // Global parameters
  defineDefaultSpace(ESpaceType::RN, 3);

  // Generate the data base

  Db* data = Db::createFillRandom(100, 3, 1, 0,0,0,0,VectorDouble(),{0,0,0}, {100,100,100});

  if (graphic)
    (void) data->dumpToNF("Data.ascii");

  // Generate the output grid
  VectorInt nx = {100,100,100};

  DbGrid* grid = DbGrid::create(nx);
  if (verbose) grid->display();

  // Create the Model
  Model* model = Model::createFromParam(ECov::EXPONENTIAL, 20);
  if (verbose) model->display();

  // Unique Neighborhood
  NeighUnique* neighU = NeighUnique::create();
  if (verbose) neighU->display();

  Timer timer;
  kriging(data, grid, model, neighU, EKrigOpt::POINT, true, false, false);
  timer.displayIntervalMilliseconds("Kriging in Unique Neighborhood", 2400);

  // Produce some statistics for comparison
  DbStringFormat* gridfmt = DbStringFormat::create(FLAG_STATS, {"*estim"});
  grid->display(gridfmt);
  delete gridfmt;

  if (graphic)
    (void) grid->dumpToNF("Grid.ascii");

  if (neighU    != nullptr) delete neighU;
  if (data      != nullptr) delete data;
  if (grid      != nullptr) delete grid;
  if (model     != nullptr) delete model;

  return (0);
}
