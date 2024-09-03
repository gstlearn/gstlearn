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
#include "API/SPDE.hpp"
#include "geoslib_old_f.h"

#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbHelper.hpp"
#include "Model/Model.hpp"
#include "Space/ASpaceObject.hpp"

/**
 * This file is meant to perform any test that needs to be coded for a quick trial
 * It will be compiled but not run nor diff'ed.
 */
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("AAA_");

  defineDefaultSpace(ESpaceType::SN);

  Db* db = Db::createFillRandom(2, 2, 1);

  DbGrid* grid = DbGrid::create({2, 2});

  Model* model = Model::createFromParam(ECov::MATERN, 1., 1., 1.);

  message("value = %lf\n", logLikelihoodSPDE(db, grid, model));

  delete db;
  delete grid;
  delete model;

  return(0);
}
