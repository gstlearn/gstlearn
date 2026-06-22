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
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Estimation/Vecchia.hpp"
#include "Model/Model.hpp"
#include "geoslib_define.h"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Do not remove
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-"); // Here set the test name

  auto ndim = 2;
  (void)defineDefaultSpace(ESpaceType::RN, ndim);

  auto ndat = 3;
  Db* data = Db::createFillRandom(ndat);
  data->getStatsAsTable().display();

  auto flagKS = false;
  Model* model = Model::createFromParam(ECov::MATERN, 0.3, 2., 1.);
  if (flagKS)
    (void)model->setMean(5.);
  else
    (void)model->setDriftIRF(1);

  auto nx = 3;
  DbGrid* grid = DbGrid::create({nx, nx}, {1. / nx, 1. / nx});

  auto nb_vecchia = 10;
  bool verbose = false;
  (void)krigingVecchia(data, grid, model, nb_vecchia, verbose);
  grid->getStatsAsTable({"Vecchia.*"}).display();
}
