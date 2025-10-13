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
#include "API/SPDE.hpp"
#include "Basic/File.hpp"
#include "Basic/FunctionalSpirale.hpp"
#include "Basic/Law.hpp"
#include "Basic/OptCustom.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighUnique.hpp"

#define __USE_MATH_DEFINES
#include <cmath>

using namespace gstlrn;

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int argc, char* argv[])

{
  OptCustom::define("ompthreads", 1);
  Id seed = 10355;
  law_set_random_seed(seed);

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  DbStringFormat dbfmt(FLAG_STATS, {"Kriging*"});

  ASerializable::setPrefixName("test_Francky-");

  // Creating the 2-D Grid
  VectorInt nx = {101, 101};
  DbGrid* grid = DbGrid::create(nx);

  // Creating the 2-D Data Db with a Normal Variable
  auto ndata = 100;
  Db* dat    = Db::createFromBox(ndata, {0., 0.}, {100., 100.}, 3243);

  // Creating the Neighborhood (Unique)
  NeighUnique* neighU = NeighUnique::create();

  // Creating the Non-stationary Model
  Model* model = Model::createFromParam(ECov::MATERN, 1., 1., 1., {10., 40.},
                                        MatrixSymmetric(), {30., 0.});

  FunctionalSpirale spirale(0., -1.4, 1., 1., 50., 50.);

  model->getCovAniso(0)->makeAngleNoStatFunctional(&spirale);

  // Simulating variable at data location (using SPDE)
  Id useCholesky = 0;
  law_set_random_seed(13256);
  (void)simulateSPDE(nullptr, dat, model, 1, useCholesky,
                     nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, SPDEParam(), false,
                     NamingConvention("Data", true, false));
  (void)dat->dumpToNF("Data.NF");

  // Testing Kriging (traditional method)
  (void)kriging(dat, grid, model, neighU, true, false);

  // Testing Kriging (with SPDE)
  (void)krigingSPDE(dat, grid, model, true, false, useCholesky);

  // Printout (optional)
  (void)grid->dumpToNF("Grid.NF");
  grid->display(&dbfmt);

  message("Test performed successfully\n");

  delete dat;
  delete grid;
  delete neighU;
  delete model;

  return 0;
}
