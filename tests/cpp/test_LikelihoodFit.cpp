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
/**
 * This file is meant to parametrized the ModelGeneric in terms of ParamInfo
 * and to fit the values of these parameters according to the Maximum LogLikelihood
 * method and using the Vecchia approximation.
 */
#include "Basic/VectorHelper.hpp"
#include "Enum/ELoc.hpp"
#include "Estimation/ALikelihood.hpp"
#include "Estimation/AModelOptimFactory.hpp"
#include "Model/Model.hpp"
#include "geoslib_define.h"
#include "utils.hpp"
using namespace gstlrn;

int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  bool verbose = true;
  bool trace   = false;
  Model* model = Model::createFromParam(ECov::EXPONENTIAL, 20, 1, 1);
  String filename;
  Db* db = nullptr;

  model->setDriftIRF(0);
  Id caset = 0;
  if (caset == 0 || caset == 1)
  {
    message("-----------------------------------------------------\n");
    message("Test on Pollution dataset\n");
    filename = getTestData("Pollution", "Pollution.dat");
    db       = Db::createFromCSV(filename, CSVformat(), false);
    db->setLocators({"X", "Y"}, ELoc::X);
    db->setLocator("Zn", ELoc::Z);
    db->setLocator("Pb");
    model->setDriftIRF(0);
    model->fitNew(db, nullptr, nullptr, nullptr, ModelOptimParam(),
                  ITEST, verbose, trace);
    model->display();
    delete db;
  }
  if (caset == 0 || caset == 2)
  {
    message("-----------------------------------------------------\n");
    message("Test on Scotland Temperatures dataset\n");
    filename = getTestData("Scotland", "Scotland_Temperatures.csv");
    db       = Db::createFromCSV(filename, CSVformat(), false);

    db->setLocators({"Longitude", "Latitude"}, ELoc::X);
    db->setLocator("January_temp", ELoc::Z);

    model->setDriftIRF(0);
    model->fitNew(db, nullptr, nullptr, nullptr, ModelOptimParam(),
                 ITEST, verbose, trace);
    model->display();
    delete db;
  }
  if (caset == 0 || caset == 3)
  {
    message("-----------------------------------------------------\n");
    message("Test on Scotland Temperatures dataset without API\n");
    filename = getTestData("Scotland", "Scotland_Temperatures.csv");
    db       = Db::createFromCSV(filename, CSVformat(), false);

    db->setLocators({"Longitude", "Latitude"}, ELoc::X);
    db->setLocator("January_temp", ELoc::Z);

    model->setDriftIRF(0);

    bool reml            = false;
    Id nb_neighVecchia   = ITEST;
    ModelOptimParam* mop = ModelOptimParam::create(false);

    auto* ll = AModelOptimFactory::create(model, db, nullptr, nullptr, nullptr,
                                          *mop, nb_neighVecchia, reml);
    ll->setVerbose(verbose, trace);

    ll->run();

    model->display();
    delete mop;
    delete db;
  }
  delete model;
  return 0;
}
