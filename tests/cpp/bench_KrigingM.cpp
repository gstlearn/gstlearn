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
#include "Basic/OptCustom.hpp"
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
  OptCustom::define("ompthreads",5);
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

  // Create the neighborhood
  int nmaxi     = 20;
  int nmini     = 2;
  int nsect     = 8;
  int nsmax     = 3;
  double radius = 1.;

  // Print the environment
  message("This test is meant to test Kriging Efficiency using Moving Neighborhood\n");
  message("- the Output Grid contains %d nodes\n",
          grid->getNSample(true));
  message("- the Neighborhood is fixed at nmaxi=%d nsect=%d radius=%lf\n",
          nmaxi, nsect, radius);
  message("- Various number of samples in the Data Set\n");
  message("- Leaf size for the Ball Tree search (when used) = %d\n", leaf_size);

  VectorInt nechs = {300, 1000, 10000, 100000};
  VectorInt timeb = {460, 490, 730, 2100};
  VectorInt times = {900, 2540, 0, 0};
  NeighMoving* neighM;
  Db* data;
  Timer timer;
  for (int icas = 0; icas < (int)nechs.size(); icas++)
  {
    int nloop = (icas < 2) ? 2 : 1;
    int nech = nechs[icas];
    for (int iloop = 0; iloop < nloop; iloop++)
    {
      int ileaf = (iloop == 0) ? leaf_size : 0;
      data      = Db::createFillRandom(nech, 2, 1);
      neighM    = NeighMoving::create(false, nmaxi, radius, nmini, nsect, nsmax,
                                      VectorDouble(), VectorDouble(), ileaf > 0, ileaf);

      timer.reset();
      message("Nsample = %7d - Leaf = %3d\n", nech, ileaf);
      kriging(data, grid, model, neighM, true, false);
      int time = (iloop == 0) ? timeb[icas] : times[icas];
      timer.displayIntervalMilliseconds("Kriging in Moving Neighborhood", time);

      delete data;
      delete neighM;
    }
  }

  // Produce some stats for comparison
  data = Db::createFillRandom(100000, 2, 1);
  grid->deleteColumn("Kriging.*");
  neighM = NeighMoving::create(false, nmaxi, radius, nmini, nsect, nsmax,
                        VectorDouble(), VectorDouble(), true, leaf_size);
  kriging(data, grid, model, neighM, true, false);
  DbStringFormat* dbfmt = DbStringFormat::create(FLAG_STATS, {"Kriging*estim"});
  grid->display(dbfmt);

  delete data;
  delete neighM;
  delete grid;
  delete model;
  delete dbfmt;

  return (0);
}
