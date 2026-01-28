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
#include "geoslib_define.h"

#include "Basic/ASerializable.hpp"
#include "Basic/OptDbg.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighMoving.hpp"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Do not remove
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-"); // Here set the test name

  Db* dbin = Db::createFromNF("/home/drenard/Rtest/DR/Meryem/dbin.NF");
  dbin->display();
  DbGrid* dbout = DbGrid::createFromNF("/home/drenard/Rtest/DR/Meryem/dbout.NF");
  dbout->display();
  Model* model = Model::createFromNF("/home/drenard/Rtest/DR/Meryem/model.NF");
  model->display();
  NeighMoving* neigh = NeighMoving::createFromNF("/home/drenard/Rtest/DR/Meryem/neigh.NF");
  neigh->display();

  OptDbg::setReference(202);
  (void)kriging(dbin, dbout, model, neigh);

  delete dbout;
  delete model;
  delete neigh;
  return 0;
}
