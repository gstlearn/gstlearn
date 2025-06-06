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

  Db* db          = Db::createFillRandom(100, 2, 0);
  Model* model    = Model::createFromParam(ECov::EXPONENTIAL, TEST, 2., 1., {0.1, 0.3}, MatrixSymmetric(), {30., 0});
  Model* modelfit1 = Model::createFromParam(ECov::EXPONENTIAL, TEST, 1, 1, {1., 1.}, MatrixSymmetric(), {0., 0});
  Model* modelfit2 = modelfit1->clone();
  mestitle(0, "Test fit likelihood");

  mestitle(1, "True Model");
  model->display();

  simtub(nullptr, db, model, nullptr, 1, 234555, 3000);
  message("Start Fitting Model with Vecchia Approximation\n");
  // Do not use 'verbose' for cross-platforms comparison
  modelfit1->fitNew(db, nullptr, nullptr, Option_AutoFit(), Option_VarioFit(),
                   true, 30, false);

  mestitle(1,"Fitted Model");
  modelfit1->display();

  message("Start Fitting Model with Likelihood\n");
  // Do not use 'verbose' for cross-platforms comparison
  modelfit2->fitNew(db, nullptr, nullptr, Option_AutoFit(), Option_VarioFit(),
                    false, 30, false);

  mestitle(1,"Fitted Model");
  modelfit2->display();
  delete db;
  delete model;
  delete modelfit1;
  delete modelfit2;
  return (0);
}
