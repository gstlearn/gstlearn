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
#include "Basic/FunctionalSpirale.hpp"
#include "Basic/File.hpp"
#include "Basic/OptCst.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Db/DbGrid.hpp"
#include "API/SPDE.hpp"
#include "Model/Model.hpp"

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
  ASerializable::setPrefixName("test_SPDEAPI-");
  int seed = 10355;
  int nbsimu = 3;
  law_set_random_seed(seed);

  // Creating the resulting Grid
  auto nx = { 101, 101 };
  DbGrid *grid = DbGrid::create(nx);

  // Creating the Model
  Model* model = Model::createFromParam(ECov::MATERN, 1., 1., 1., {10.,30.});

  // Creating non-stationarity field (spiral) and attaching it to the Model
  FunctionalSpirale spirale(0., -1.4, 1., 1., 50., 50.);
  model->getCovAniso(0)->makeAngleNoStatFunctional(&spirale);
  model->display();

  // Creating Data
  int ndata = 100;
  Db* dat = Db::createFromBox(ndata, {0.,0.}, {100.,100.}, 43246);
  VectorDouble z = VH::simulateGaussian(ndata);
  int useCholesky = 0;
  law_set_random_seed(132341);
  (void)simulateSPDE(nullptr, dat, model, nullptr, 1, nullptr, useCholesky,
                     SPDEParam(), false, false,
                     NamingConvention("variable", false, false));
  dat->display();

  // Estimation and simulations
  (void)krigingSPDE(dat, grid, model, nullptr, true, false, nullptr,
                    useCholesky, SPDEParam(), 0, false, false,
                    NamingConvention("K-spirale"));
  law_set_random_seed(132341);

  (void)simulateSPDE(nullptr, grid, model, nullptr, nbsimu, nullptr,
                     useCholesky, SPDEParam(), false, false,
                     NamingConvention("NCS-spirale"));
  law_set_random_seed(132341);
  (void)simulateSPDE(dat, grid, model, nullptr, nbsimu, nullptr, useCholesky,
                     SPDEParam(), false, false,
                     NamingConvention("CDS-spirale"));

  (void) grid->dumpToNF("grid.ascii");
  DbStringFormat dbfmt(FLAG_STATS,{"spde*"});
  // To prevent diff between some platforms (round to 10^-2)
  OptCst::define(ECst::NTDEC, 2);
  grid->display(&dbfmt);

  delete dat;
  delete grid;
  delete model;
  return 0;
}

