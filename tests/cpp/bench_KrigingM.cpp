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
#include "Basic/VectorNumT.hpp"
#include "Enum/ESpaceType.hpp"
#include "Enum/ECov.hpp"

#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Basic/File.hpp"
#include "Basic/Timer.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Estimation/CalcKriging.hpp"

void st_test(Db* grid, Model* model, int nech, int leaf_size, bool verbose)
{
  // Generate the data base

  int ndim = 2;
  int nvar = 1;
  Db* data = Db::createFillRandom(nech, ndim, nvar);

  // Moving Neighborhood
  int nmaxi     = 20;
  int nmini     = 2;
  int nsect     = 8;
  int nsmax     = 3;
  double radius = 1.;
  bool useBall = (leaf_size > 0);

  NeighMoving* neighM =
    NeighMoving::create(false, nmaxi, radius, nmini, nsect, nsmax,
    VectorDouble(), VectorDouble(), useBall, leaf_size);

  if (verbose) neighM->display();

  Timer timer;
  kriging(data, grid, model, neighM, true, false);

  if (verbose)
    timer.displayIntervalMilliseconds("Kriging in Moving Neighborhood", 1500);
  else
  {
    double msec = timer.getIntervalMilliseconds(true);
    message("Nsample = %7d - Leaf = %3d - Time = %10d ms\n", nech, leaf_size,
            (int) msec);
  }

  delete neighM;
  delete data;
}

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
  ASerializable::setPrefixName("BenchKrigingM-");

  // Global parameters
  int ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);
  int leaf_size = 30;

  // Generate the output grid
  int ncell       = 100;
  VectorInt nx    = {ncell, ncell};
  VectorDouble dx = {1. / ncell, 1. / ncell};
  DbGrid* grid    = DbGrid::create(nx, dx);

  // Create the Model
  double range = 1. / 5.;
  double sill  = 2.;
  Model* model = Model::createFromParam(ECov::SPHERICAL, range, sill);

  message("This test is meant to test Kriging Efficiency using Moving "
          "Neighborhood\n");
  message("- the Output Grid contains %d nodes\n",
          grid->getNSample(true));
  message("- Various number of samples in the Data Set\n");
  message("- Leaf size for the Ball Tree search (0: no Ball Tree search)\n");

  VectorInt nechs = {300, 1000, 10000, 100000};

  for (int icas = 0; icas < (int)nechs.size(); icas++)
  {
    int nloop = (icas < 2) ? 2 : 1;
    for (int iloop = 0; iloop < nloop; iloop++)
    {
      int ileaf = (iloop == 0) ? leaf_size : 0;
      st_test(grid, model, nechs[icas], ileaf, false);
    }
  }

  grid->deleteColumn("Kriging.*");
  st_test(grid, model, nechs[3], leaf_size, true);

  // Produce some stats for comparison
  DbStringFormat* dbfmt = DbStringFormat::create(FLAG_STATS, {"Kriging*estim"});
  grid->display(dbfmt);

  delete grid;
  delete model;
  delete dbfmt;

  return (0);
}
