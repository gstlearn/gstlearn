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
#include "geoslib_d.h"
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
#include "Neigh/NeighBench.hpp"
#include "Stats/Classical.hpp"
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

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("BenchSPDE-");

  // Global parameters
  defineDefaultSpace(ESpaceType::RN, 2);
  int seed = 123;
  int nsim = 10;
  int ndat = 50;
  int nxref = 101;
  double matern_param = 1.0;

  setGlobalFlagEigen(true);
  OptCst::define(ECst::NTDEC, 2);
  OptCst::define(ECst::NTROW, -1);
  bool flagExhaustiveTest = false;

  // Feature to be tested:
  // 0: all of them
  // 1: Kriging
  // 2: non-conditional simulations
  // 3: conditional simulations
  int mode = 0;
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

  for (int ncov = 0; ncov < 2; ncov++)
  {
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

    for (int ifois = 0; ifois < 2; ifois++)
    {
      int useCholesky;
      String option;
      if (ifois == 0)
      {
        useCholesky = 0;
        option = ".NoChol";
      }
      else
      {
        useCholesky = 1;
        option = ".Chol";
      }
      if (showStats)
        message("- Cholesky Option        = %d\n", useCholesky);

      // Kriging
      if (mode == 0 || mode == 1)
      {
        timer.reset();
        String namconv = "Kriging" + option + sncov;
        (void) krigingSPDE(dat, grid, model, true, false, false, nullptr,
                           useCholesky, SPDEParam(), verbose, showStats,
                           NamingConvention(namconv));
        timer.displayIntervalMilliseconds(namconv, 400);
      }

      // Non-conditional simulations
      if (mode == 0 || mode == 2)
      {
        timer.reset();
        String namconv = "Simu.NC" + option + sncov;
        (void) simulateSPDE(NULL, grid, model, nsim, NULL, useCholesky,
                            SPDEParam(), seed, verbose, showStats,
                            NamingConvention(namconv));
        timer.displayIntervalMilliseconds(namconv, 1350);
      }

      // Conditional simulations
      if (mode == 0 || mode == 3)
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
