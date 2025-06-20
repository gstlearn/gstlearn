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

#include "Basic/AStringable.hpp"
#include "Basic/File.hpp"
#include "Db/Db.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Model/Model.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Variogram/Vario.hpp"
#include "Variogram/VarioParam.hpp"
#include "Variogram/VMap.hpp"

/**
 * This file is meant to parametrized the ModelGeneric in terms of ParamInfo
 * and to fit the values of these parameters starting from a variogram Map
 * and using the NlOpt machinery
 */
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("OptimVMap-");

  Model* model    = Model::createFromParam(ECov::EXPONENTIAL, TEST, 2., 1., {0.1, 0.3}, MatrixSymmetric(), {30., 0});
  Model* modelfit = Model::createFromParam(ECov::EXPONENTIAL, TEST, 1, 1, {1., 1}, MatrixSymmetric(), {0, 0});

  mestitle(0, "Test fit from Variogram");
  mestitle(1, "True Model");
  model->display();

  // Simulating the Input file
  int nx         = 100;
  double dx      = 1. / nx;
  DbGrid* dbgrid = DbGrid::create({nx, nx}, {dx, dx});
  (void)simtub(nullptr, dbgrid, model);
  (void)dbgrid->dumpToNF("dbgrid.ascii");

  // Calculating the experimental variogram Map
  DbGrid* dbmap = db_vmap(dbgrid, ECalcVario::VARIOGRAM, {50, 50});
  (void)dbmap->dumpToNF("VMap.ascii");

  mestitle(1, "Initial Model");
  modelfit->display();

  // Fit the Model
  ModelOptimParam mop = ModelOptimParam();
  mop.setFlagGoulard(true);
  bool verbose = true;
  bool trace = true;
  modelfit->fitNew(nullptr, nullptr, dbmap, nullptr, mop,
                   ITEST, verbose, trace);

  // Fitting procedure
  mestitle(1, "Fitted Model");
  modelfit->display();
  modelfit->dumpToNF("model.ascii");

  delete dbgrid;
  delete dbmap;
  delete model;
  delete modelfit;
  return (0);
}
