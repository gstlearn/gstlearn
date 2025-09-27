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
#include "Covariances/CovAniso.hpp"
#include "Enum/ELoc.hpp"
#include "Enum/ESpaceType.hpp"
#include "Estimation/ALikelihood.hpp"
#include "Estimation/AModelOptimFactory.hpp"
#include "Model/Model.hpp"
#include "Space/ASpaceObject.hpp"
#include "geoslib_define.h"
#include "utils.hpp"
using namespace gstlrn;

int main(int argc, char* argv[])
{
  DECLARE_UNUSED(argc);
  DECLARE_UNUSED(argv);

  Model* model = Model::createFromParam(ECov::EXPONENTIAL, 20, 1, 1);
  model->setDriftIRF(0);
  Id caset = 3;
  if (caset == 0 || caset == 1)
  {
    String filename = getTestData("Pollution", "Pollution.dat");
    Db* db          = Db::createFromCSV(filename, CSVformat(), false);
    db->setLocators({"X","Y"}, ELoc::X);
    db->setLocator("Zn", ELoc::Z);
    db->setLocator("Pb");
    Model* model = Model::createFromParam(ECov::EXPONENTIAL, 20, 1, 1);
    model->setDriftIRF(0);
    model->fitNew(db, nullptr, nullptr, nullptr, ModelOptimParam(),
                  ITEST, true, true);

  }
  else if (caset == 0 || caset == 2)
  {
    String filename = getTestData("Scotland", "Scotland_Temperatures.csv");
    Db* db          = Db::createFromCSV(filename, CSVformat(), false);

    db->setLocators({"Longitude", "Latitude"}, ELoc::X);
    db->setLocator("January_temp", ELoc::Z);

    Model* model = Model::createFromParam(ECov::EXPONENTIAL, 20, 1, 1);
    model->setDriftIRF(0);
    model->fitNew(db, nullptr, nullptr, nullptr, ModelOptimParam(),
                  ITEST, false, false);
    model->display();
  }
  else if (caset == 0 || caset == 3)
  {
    String filename = getTestData("Scotland", "Scotland_Temperatures.csv");
    Db* db          = Db::createFromCSV(filename, CSVformat(), false);

    db->setLocators({"Longitude", "Latitude"}, ELoc::X);
    db->setLocator("January_temp", ELoc::Z);

    Model* model = Model::createFromParam(ECov::EXPONENTIAL, 20, 1, 1);
    model->setDriftIRF(0);

    bool verbose =  false;
    bool trace = false;
    bool reml = false;
    Id nb_neighVecchia = ITEST;
    ModelOptimParam* mop = ModelOptimParam::create(false); 

    auto* ll = AModelOptimFactory::create(model,db,nullptr,nullptr, nullptr,
                                          *mop, nb_neighVecchia, reml);
    ll->setVerbose(verbose, trace);
  
    double cost = ll->run(); 
    message("Cost %lf\n", cost);


    VH::dump("beta", (dynamic_cast<ALikelihood*>(ll))->getBeta());
    model->display();
  }
}
