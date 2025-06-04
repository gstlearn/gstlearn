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

/**
 * This file is meant to parametrized the ModelGeneric in terms of ParamInfo
 * and to fit the values of these parameters according to the Maximum LogLikelihood
 * method and using the Vecchia approximation.
 */
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  bool useVecchia = false;
  Db* db          = Db::createFillRandom(100, 2, 0);
  Model* model    = Model::createFromParam(ECov::EXPONENTIAL, TEST, 2., 1., {0.1, 0.3}, MatrixSymmetric(), {30., 0});
  Model* modelfit = Model::createFromParam(ECov::EXPONENTIAL, TEST, 1, 1, {1., 1.}, MatrixSymmetric(), {0., 0});

  mestitle(0, "Test fit likelihood");

  mestitle(1, "True Model");
  model->display();

  simtub(nullptr, db, model, nullptr, 1, 234555, 3000);

  modelfit->fitLikelihood(db, useVecchia, false);

  mestitle(1,"Fitted Model");
  modelfit->display();

  delete db;
  delete model;
  delete modelfit;
  return (0);
}
