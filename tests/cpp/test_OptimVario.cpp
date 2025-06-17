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

/**
 * This file is meant to parametrized the ModelGeneric in terms of ParamInfo
 * and to fit the values of these parameters starting from an experimental variogram
 * and using the NlOpt machinery
 */
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("OptimVario-");

  Db* db          = Db::createFillRandom(1000, 2, 0);
  Model* model    = Model::createFromParam(ECov::EXPONENTIAL, TEST, 2., 1., {0.1, 0.3}, MatrixSymmetric(), {30., 0});
  Model* modelfit = Model::createFromParam(ECov::EXPONENTIAL, TEST, 1, 1, {1., 1}, MatrixSymmetric(), {0, 0});

  mestitle(0, "Test fit from Variogram");

  mestitle(1, "True Model");
  model->display();

  simtub(nullptr, db, model, nullptr, 1, 234555, 3000);

  // Calculating the experimental variogram
  double diagonal = db->getExtensionDiagonal();
  int nlag = 10;
  double dlag = diagonal / 2. / nlag;
  VarioParam* varioparam = VarioParam::createMultiple(4, nlag, dlag);
  Vario* vario = Vario::computeFromDb(*varioparam, db);
  vario->dumpToNF("vario.ascii");

  mestitle(1, "Initial Model");
  modelfit->display();

  // Fitting procedure
  ModelOptimParam mop = ModelOptimParam();
  mop.setFlagGoulard(true);
  modelfit->fitNew(nullptr, vario, nullptr, nullptr, mop,
                   ITEST, false, false);

  mestitle(1, "Fitted Model");
  modelfit->display();
  modelfit->dumpToNF("model.ascii");

  delete db;
  delete vario;
  delete varioparam;
  delete model;
  delete modelfit;
  return (0);
}
