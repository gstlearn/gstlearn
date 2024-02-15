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
/**
 * This function is meant to evaluate the bench marks on the SPDE functionalities
 *
 */
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
#include "Basic/OptCst.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Neigh/NeighBench.hpp"
#include "Stats/Classical.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "API/SPDE.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int argc, char *argv[])
{
  Timer timer;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(false);
  ASerializable::setPrefixName("BenchSPDE-");

  // Global parameters
  defineDefaultSpace(ESpaceType::RN, 2);
  int seed = 123;
  int nsim = 10;
  int ndat = 50;
  int nxref = 101;
  double matern_param = 1.0;

  setGlobalFlagEigen(true);
  message("Use of Eigen Library = %d\n",isGlobalFlagEigen());

  OptCst::define(ECst::NTDEC, 2);
  OptCst::define(ECst::NTROW, -1);
  bool flagExhaustiveTest = false;

  // Feature to be tested:
  // -1: all of them
  //  0: Establishing Shift Operator only
  //  1: Kriging
  //  2: non-conditional simulations
  //  3: conditional simulations
  int mode = -1;

  int nfois = 2;
  // Feature to be tested:
  // -1: all cases
  //  0: not using the Cholesky option
  //  1: using the Cholesky option
  int ifois_ref = -1;

  int ncov_tot = 2;
  // Feature to be tested:
  // -1: all the covariances
  //  0: only the case with one covariance
  //  1: only the case with two covariances
  int ncov_ref = -1;

  bool verbose = false;
  bool showStats = false;

  // Generate the data base
  Db* dat = Db::createFillRandom(ndat);

  // Generate the output grid
  DbGrid* grid = DbGrid::create({nxref, nxref}, {1./(nxref-1), 1./(nxref-1)});

  // Printout of general environment
  if (showStats)
  {
    mestitle(1, "Statistics");
    message("- Number of active data  = %d\n", ndat);
    message("- Number of target sites = %d\n", nxref * nxref);
    message("- Number of simulations  = %d\n", nsim);
  }

  // Loop for usage of Cholesky

  for (int ncov = 0; ncov < ncov_tot; ncov++)
  {
    if (ncov_ref >= 0 && ncov_ref != ncov) continue;

    // Generate the Model
    Model *model;
    model = Model::createFromParam(ECov::BESSEL_K, TEST, 1, matern_param,
                                   { 0.1, 0.3 }, VectorDouble(), { 30., 0. });
    if (ncov >= 1)
      model->addCovFromParam(ECov::BESSEL_K, TEST, 1, matern_param,
                             { 0.3, 0.2 }, VectorDouble(), { -10., 0.});
    String sncov = (ncov == 0) ? "1" : "2";

    // Printout of general environment
    if (showStats)
      message("- Number of covariances  = %d\n", model->getCovaNumber());

    // Building Shift Operator
    if (mode < 0 || mode == 0)
    {
      MeshETurbo mesh(grid);
      for (int icov = 0; icov <= ncov; icov++)
      {
        timer.reset();
        ShiftOpCs shiftop(&mesh, model, nullptr, 0, icov);
        timer.displayIntervalMilliseconds("Establishing S", 150);
      }
    }

    for (int ifois = 0; ifois < nfois; ifois++)
    {
      if (ifois_ref >= 0 && ifois != ifois_ref) continue;

      int useCholesky = ifois;
      String option = (ifois == 0) ? ".NoChol" : ".Chol";
      if (showStats)
        message("- Cholesky Option        = %d\n", useCholesky);

       // Kriging
      if (mode < 0 || mode == 1)
      {
        timer.reset();
        String namconv = "Kriging" + option + sncov;
        (void) krigingSPDE(dat, grid, model, true, false, false, nullptr,
                           useCholesky, SPDEParam(), verbose, showStats,
                           NamingConvention(namconv));
        timer.displayIntervalMilliseconds(namconv, 400);
      }

      // Non-conditional simulations
      if (mode < 0 || mode == 2)
      {
        timer.reset();
        String namconv = "Simu.NC" + option + sncov;
        (void) simulateSPDE(NULL, grid, model, nsim, NULL, useCholesky,
                            SPDEParam(), seed, verbose, showStats,
                            NamingConvention(namconv));
        timer.displayIntervalMilliseconds(namconv, 1350);
      }

      // Conditional simulations
      if (mode < 0 || mode == 3)
      {
        timer.reset();
        String namconv = "Simu.CD" + option + sncov;
        (void) simulateSPDE(dat, grid, model, nsim, NULL, useCholesky,
                            SPDEParam(), seed, verbose, showStats,
                            NamingConvention(namconv));
        timer.displayIntervalMilliseconds(namconv, 3130);
      }
    }
    delete model;
  }

  // Produce some statistics for comparison
  dbStatisticsPrint(grid, { "Kriging*", "Simu*" },
                    EStatOption::fromKeys( { "MINI", "MAXI", "MEAN", "STDV" }));
  if (flagExhaustiveTest)
  {
    DbStringFormat *dbfmt = DbStringFormat::createFromFlags(false, false, false, false, true);
    grid->display(dbfmt);
  }
  (void) grid->dumpToNF("Grid.ascii");

  if (dat       != nullptr) delete dat ;
  if (grid      != nullptr) delete grid;

  return (0);
}
