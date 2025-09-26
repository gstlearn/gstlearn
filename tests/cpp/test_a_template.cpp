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
#include "Model/Model.hpp"
#include "Space/ASpaceObject.hpp"
#include "geoslib_define.h"
#include "utils.hpp"
using namespace gstlrn;

int main(int argc, char* argv[])
{
  DECLARE_UNUSED(argc);
  DECLARE_UNUSED(argv);

  String filename = getTestData("Scotland", "Scotland_Temperatures.csv");
  Db* db          = Db::createFromCSV(filename, CSVformat(), false);

  db->setLocators({"Longitude", "Latitude"}, ELoc::X);
  db->setLocator("January_temp", ELoc::Z);

  Model* model = Model::createFromParam(ECov::EXPONENTIAL, 20, 1, 1);
  model->setDriftIRF(0);
  model->fitNew(db, nullptr, nullptr, nullptr, ModelOptimParam(),
                ITEST, true, false);
}
