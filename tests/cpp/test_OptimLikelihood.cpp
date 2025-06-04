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
#include <iostream>
#include "LinearOp/CholeskySparseInv.hpp"

/**
 * This file is meant to perform any test that needs to be coded for a quick trial
 * It will be compiled but not run nor diff'ed.
 */
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  bool useVecchia = false;
  Db* db = Db::createFillRandom(100, 2, 0);
  Model* model = Model::createFromParam(ECov::EXPONENTIAL, TEST, 2., 1., {0.1, 0.3}, MatrixSymmetric(), {30., 0});
  Model* modelfit = Model::createFromParam(ECov::EXPONENTIAL, TEST, 1, 1, {1., 1.}, MatrixSymmetric(), {0., 0});


  message("Test fit likelihood\n");
  message("True Model\n");
  model->display();
  simtub(nullptr, db, model, nullptr, 1, 234555, 3000);
  db->display();
  modelfit->fitLikelihood(db, useVecchia, true);
  message("Fitted Model\n");
  modelfit->display();

  delete db;
  delete model;
  delete modelfit;
  return (0);
}
