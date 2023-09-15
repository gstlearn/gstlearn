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

  // Global parameters
  defineDefaultSpace(ESpaceType::RN, 2);
  int seed = 123;
  int nsim = 10;
  int ndat = 50;
  int nxref = 101;
  double matern_param = 1.0;

  // Feature to be tested:
  // 0: all of them
  // 1: Kriging
  // 2: non-conditional simulations
  // 3: conditional simulations
  int mode = 1;

  // Generate the data base
  Db* dat = Db::createFillRandom(ndat);

  // Generate the Model
  Model* model = Model::createFromParam(ECov::BESSEL_K, TEST, 1, matern_param,
                                        { 0.1, 0.3 }, VectorDouble(), { 30., 0. });

  // Generate the output grid
  DbGrid* grid = DbGrid::create({nxref, nxref}, {1./(nxref-1), 1./(nxref-1)});

  // Evaluate Kriging
  if (mode == 0 || mode == 1)
  {
    timer.reset();
    (void) krigingSPDE(dat, grid, model, true, false, false);
    timer.displayIntervalMilliseconds("Kriging", 400);
  }

  // Evaluate non-conditional simulations
  if (mode == 0 || mode == 2)
  {
    timer.reset();
    (void) simulateSPDE(NULL, grid, model, nsim, NULL, 0, 11, 18, 8, seed, 1.e-2,
                        false, NamingConvention("Simu.NC"));
    timer.displayIntervalMilliseconds("Non-conditional simulations", 1350);
  }

  // Evaluate conditional simulations
  if (mode == 0 || mode == 3)
  {
    timer.reset();
    (void) simulateSPDE(dat, grid, model, nsim, NULL, 0, 11, 18, 8, seed, 1.e-2,
                        false, NamingConvention("Simu.CD"));
    timer.displayIntervalMilliseconds("Conditional simulations", 3130);
  }

  // Produce some stats for comparison
  dbStatisticsPrint(grid, {"Simu*"}, EStatOption::fromKeys({"MEAN", "STDV"}));

  if (dat       != nullptr) delete dat ;
  if (grid      != nullptr) delete grid;
  if (model     != nullptr) delete model;

  return (0);
}
