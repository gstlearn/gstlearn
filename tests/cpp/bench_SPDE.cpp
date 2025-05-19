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
#include "Basic/Law.hpp"
#include "Enum/ESpaceType.hpp"
#include "Enum/ECov.hpp"

#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Basic/File.hpp"
#include "Basic/Timer.hpp"
#include "Basic/OptCst.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Neigh/NeighBench.hpp"
#include "Stats/Classical.hpp"
#include "LinearOp/ShiftOpMatrix.hpp"
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
  int seed  = 123;
  int nsim  = 10;
  int ndat  = 50;
  int nxref = 101;
  double matern_param = 1.0;
  setGlobalFlagEigen(true);
  message("Use of Eigen Library = %d\n",isGlobalFlagEigen());

  OptCst::define(ECst::NTDEC, 2);
  OptCst::define(ECst::NTROW, -1);
  bool flagExhaustiveTest = false;
  bool flagStatistics     = true;

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
  //  0: one covariance
  //  1: two covariances
  int ncov_ref = -1;

  bool showStats = false;
  bool verbose   = false;
  bool flagOld   = false;
  bool flagStd   = false;
  int nbMC       = 10;

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
    model = Model::createFromParam(ECov::MATERN, TEST, 1, matern_param,
                                   { 0.1, 0.3 }, MatrixSymmetric(), { 30., 0. });
    if (ncov >= 1)
      model->addCovFromParam(ECov::MATERN, TEST, 1, matern_param,
                             { 0.3, 0.2 }, MatrixSymmetric(), { -10., 0.});
    String sncov = (ncov == 0) ? "1" : "2";

    // Printout of general environment
    if (showStats)
      message("- Number of covariances  = %d\n", model->getNCov());

    // Building Shift Operator
    if (mode < 0 || mode == 0)
    {
      MeshETurbo mesh(grid);
      for (int icov = 0; icov <= ncov; icov++)
      {
        timer.reset();
        ShiftOpMatrix shiftop(&mesh, model->getCovAniso(icov), nullptr);
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
        String namconv;
        namconv.append("Kriging");
        namconv.append(option);
        namconv.append(sncov);
        law_set_random_seed(13243);
        if (flagOld)
          (void)krigingSPDEOld(dat, grid, model, nullptr, true, flagStd, nullptr,
                            useCholesky, SPDEParam(), nbMC, verbose, showStats,
                            NamingConvention(namconv));
        else
          (void)krigingSPDE(dat, grid, model, true, true, useCholesky,
                            VectorMeshes(), nullptr, SPDEParam(),
                            NamingConvention(namconv));
        timer.displayIntervalMilliseconds(namconv, 400);
      }

      // Non-conditional simulations
      if (mode < 0 || mode == 2)
      {
        timer.reset();
        String namconv;
        namconv.append("Simu.NC");
        namconv.append(option);
        namconv.append(sncov);
        law_set_random_seed(seed);
        if (flagOld)
          (void)simulateSPDEOld(NULL, grid, model, nullptr, nsim, NULL, useCholesky,
                                SPDEParam(), verbose, showStats,
                                NamingConvention(namconv));
        else
          (void)simulateSPDE(nullptr, grid, model, nsim, useCholesky,
                             VectorMeshes(), nullptr, SPDEParam(),
                             NamingConvention(namconv));
        timer.displayIntervalMilliseconds(namconv, 1350);
      }

      // Conditional simulations
      if (mode < 0 || mode == 3)
      {
        timer.reset();
        String namconv;
        namconv.append("Simu.CD");
        namconv.append(option);
        namconv.append(sncov);
        law_set_random_seed(seed);
        if (flagOld)
          (void)simulateSPDEOld(dat, grid, model, nullptr, nsim, NULL, useCholesky,
                                SPDEParam(), verbose, showStats,
                                NamingConvention(namconv));
        else
          (void)simulateSPDE(dat, grid, model, nsim, useCholesky,
                             VectorMeshes(), nullptr, SPDEParam(),
                             NamingConvention(namconv));
        timer.displayIntervalMilliseconds(namconv, 3130);
      }
    }
    delete model;
  }

  // Produce some statistics for comparison
  if (flagStatistics)
    dbStatisticsPrint(grid, { "Kriging*", "Simu*" },
                      EStatOption::fromKeys( { "MINI", "MAXI", "MEAN", "STDV" }));
  if (flagExhaustiveTest)
  {
    DbStringFormat *dbfmt = DbStringFormat::createFromFlags(false, false, false, false, true);
    grid->display(dbfmt);
  }
  (void) grid->dumpToNF("Grid.ascii");

  delete dat ;
  delete grid;

  return (0);
}
