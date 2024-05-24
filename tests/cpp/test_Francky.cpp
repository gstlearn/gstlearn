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
#include "Basic/VectorHelper.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "Basic/FunctionalSpirale.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "LinearOp/PrecisionOpMultiConditional.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Model/NoStatArray.hpp"
#include "Model/NoStatFunctional.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Estimation/CalcKriging.hpp"
#include "API/SPDE.hpp"

#include <math.h>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#define __USE_MATH_DEFINES
#include <cmath>

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int argc, char *argv[])

{
  int seed = 10355;
  law_set_random_seed(seed);

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  DbStringFormat dbfmt(FLAG_STATS,{"Kriging*"});

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("test_Francky-");

  // Creating the 2-D Grid
  auto nx = { 101, 101 };
  DbGrid* grid = DbGrid::create(nx);

  // Creating the 2-D Data Db with a Normal Variable
  auto ndata = 100;
  Db* dat = Db::createFromBox(ndata, {0.,0.}, {100.,100.}, 3243);

  // Creating the Neighborhood (Unique)
  auto neighU = std::shared_ptr<ANeigh>(NeighUnique::create());

  // Creating the Non-stationary Model
  Model* model = Model::createFromParam(ECov::BESSEL_K, 1., 1., 1., {10., 40.}, VectorDouble(), {30., 0.});
  FunctionalSpirale spirale(0., -1.4, 1., 1., 50., 50.);
  NoStatFunctional NoStat(&spirale);
  model->addNoStat(&NoStat);

  // Simulating variable at data location (using SPDE)
  int useCholesky = 0;
  (void) simulateSPDE(nullptr, dat, model, 1, nullptr, useCholesky, SPDEParam(), 13256, false, false,
                      NamingConvention("Data", true, false));
  (void) dat->dumpToNF("Data.ascii");

  // Testing Kriging (with SPDE)
  (void) krigingSPDE(dat, grid, model, true, false, nullptr, useCholesky, SPDEParam());

  // Testing Kriging (traditional method)
  (void) kriging(dat, grid, model, neighU);

  // Printout (optional)
  (void) grid->dumpToNF("Grid.ascii");
  grid->display(&dbfmt);

  message("Test performed successfully\n");

  delete dat;
  delete grid;
  delete model;

  return 0;
}
