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

#include "Basic/File.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Estimation/CalcGlobal.hpp"

/**
 * This file is meant to perform any test that needs to be coded for a quick trial
 * It will be compiled but not run nor diff'ed.
 */
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  Db myDb = Db();
  myDb.addColumns({0, 1}, "longitude", ELoc::X, 0);
  myDb.addColumns({0, 1}, "latitude", ELoc::X, 1);
  myDb.addColumns({3, 7}, "density", ELoc::Z, 0);
  myDb.display();

  DbGrid* myGrid = DbGrid::create({10, 10});
  myGrid->display();

  Model* myModel = Model::createFromParam(ECov::LINEAR, 10, 1);

  Global_Result res = global_kriging(&myDb, myGrid, myModel, 0, true);

  delete myGrid;
  delete myModel;
  return (0);
}
