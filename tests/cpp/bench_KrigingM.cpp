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

  NeighMoving* neighM =
    NeighMoving::create(false, nmaxi, radius, nmini, nsect, nsmax);
  if (leaf_size > 0) neighM->setBallSearch(true, leaf_size);

  if (verbose)
  {
    // Print the test environment
    message("This test is meant to test Kriging using Moving Neighborhood\n");
    message("- the Data Set contains %d samples\n",
            data->getSampleNumber(true));
    message("- the Output Grid contains %d nodes\n",
            grid->getSampleNumber(true));
    message("- the Moving Neighborhood is required:\n");
    message("  . Radius dimension = %lf\n", radius);
    message("  . Maximum number of neighbors = %d\n", nmaxi);
    message("  . Minimum number of neiNumber of Rowsghbors = %d\n", nmini);
    message("  . Number of angular sectors   = %d\n", nsect);
    message("  . Maxmimum number of neighbors per sector = %d\n", nsmax);
    if (leaf_size > 0)
      message("  . Leaf Size for Ball Tree algorithm = %d\n", leaf_size);
  }

  Timer timer;
  kriging(data, grid, model, neighM, EKrigOpt::POINT, true, false);

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
  bool graphic = true;
  bool onlyOne = true;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("BenchKrigingM-");

  // Global parameters
  int ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);
  OptCustom::define("oldStyle", 0.);

  // Generate the output grid
  int ncell       = 100;
  VectorInt nx    = {ncell, ncell};
  VectorDouble dx = {1. / ncell, 1. / ncell};
  DbGrid* grid    = DbGrid::create(nx, dx);

  // Create the Model
  double range = 1. / 5.;
  double sill  = 2.;
  bool verbose = false;
  Model* model = Model::createFromParam(ECov::SPHERICAL, range, sill);

  if (onlyOne)
  {
    int nech      = 1000;
    int leaf_size = 30;
    if (verbose) OptDbg::setReference(1);
    st_test(grid, model, nech, leaf_size, true);

    // Produce some stats for comparison
    DbStringFormat* dbfmt = DbStringFormat::create(FLAG_STATS, {"*estim"});
    grid->display(dbfmt);
    delete dbfmt;
    if (graphic) (void)grid->dumpToNF("Grid.ascii");
  }
  else
  {
    message("This test is meant to test Kriging Efficiency using Moving "
            "Neighborhood\n");
    message("- the Output Grid contains %d nodes\n",
            grid->getSampleNumber(true));
    message("- Various number of samples in the Data Set\n");
    message("- Various dimensions of Leaf sizes (for Ball Tree search)\n");
    message("  (0: no Ball Tree search)\n");

    VectorInt nechs = {300, 1000, 10000, 100000};
    VectorInt leafs = {0, 10, 20, 30, 50, 100, 200};

    for (int iech = 0; iech < (int)nechs.size(); iech++)
    {
      for (int ileaf = 0; ileaf < (int)leafs.size(); ileaf++)
        st_test(grid, model, nechs[iech], leafs[ileaf], false);
      message("\n");
    }
  }

  delete grid;
  delete model;

  return (0);
}
