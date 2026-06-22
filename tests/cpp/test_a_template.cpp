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
#include "Basic/ASerializable.hpp"
#include "Basic/OptCst.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Estimation/Estimations.hpp"
#include "Estimation/Vecchia.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighUnique.hpp"
#include "geoslib_define.h"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Do not remove
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-"); // Here set the test name

  OptCst::defineByKey("ASP", 0);
  auto ndim = 2;
  (void)defineDefaultSpace(ESpaceType::RN, ndim);

  auto ndat = 3;
  Db* data = Db::createFillRandom(ndat);

  Model* model = Model::createFromParam(ECov::MATERN, 0.3, 2., 1.);
  (void)model->setMean(5.);

  auto nx = 3;
  DbGrid* grid = DbGrid::create({nx, nx}, {1. / nx, 1. / nx});

  // auto* neighU = NeighUnique::create();
  // (void)kriging(data, grid, model, neighU, true, true);

  // auto* neighM = NeighMoving::create(false, 20, 1., 20);
  // (void)kriging(data, grid, model, neighM, true, true);

  auto nb_vecchia = 50;
  (void)krigingVecchia(data, grid, model, nb_vecchia);
  grid->getStatsAsTable({"Vecchia.*"}).display();
}
