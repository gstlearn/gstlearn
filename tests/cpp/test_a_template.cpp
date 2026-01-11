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
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Model/Model.hpp"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Do not remove
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-"); // Here set the test name

  auto* dbnostat1 = DbGrid::create({10, 10});
  auto* dbnostat2 = DbGrid::create({10, 10});
  auto tab1       = dbnostat1->getColumns({"x1"});
  auto tab2       = dbnostat2->getColumns({"x1"});
  tab1.addCst(1.);
  dbnostat1->addColumns(tab1, "sills");
  dbnostat2->addColumns(tab2, "angles");
  tab2.addCst(2.);
  dbnostat2->addColumns(tab2, "ranges");

  auto* model1 = Model::createFromParam(ECov::NUGGET);
  model1->addCovFromParam(ECov::EXPONENTIAL, EPSILON6, 1., 1.2, {10, 20}, {}, {90, 0});

  model1->attachNoStatDb(dbnostat1);
  model1->getCovAniso(0)->makeSillNoStatDb("sills", 0, 0, dbnostat1);
  model1->getCovAniso(1)->makeAngleNoStatDb("angles", 0, dbnostat2);
  model1->getCovAniso(1)->makeRangeNoStatDb("ranges", 1, dbnostat2);

  Id ndat  = 5;
  auto* db = Db::createFillRandom(ndat);

  mestitle(0, "Writing the contents of the Model in a Neutral File");
  model1->dumpToNF("dummy");
  model1->evalCovMatSym(db).display();
  model1->display();

  mestitle(0, "Reading the contents of the Model from a Neutral File");
  auto* model2 = Model::createFromNF("dummy");
  model2->evalCovMat(db).display();
  model2->display();

  delete dbnostat1;
  delete dbnostat2;
  delete db;
  delete model1;
  delete model2;
  return 0;
}
