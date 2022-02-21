/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_d.h"
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "Db/Db.hpp"
#include "Variogram/VarioParam.hpp"
#include "Variogram/Vario.hpp"
#include "Model/Model.hpp"
#include "Basic/Table.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Polygon/Polygons.hpp"
#include "LithoRule/Rule.hpp"

/****************************************************************************/
/*!
** Main Program for testing the sparse matrix algebra
**
*****************************************************************************/
int main(int /*argc*/, char */*argv*/[])

{
  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("TS-");

  // =======================
  // Checking Db
  // =======================

  // ===== Create the Db db1
  int nech = 20;
  int ndim = 2;
  bool verbose = false;

  Db* db1 = Db::createFromBox(nech,VectorDouble(),VectorDouble(),ndim);
  VectorDouble vec1 = ut_vector_simulate_gaussian(nech);
  db1->addColumns(vec1,"myvar1",ELoc::Z, 0);
  VectorDouble vec2 = ut_vector_simulate_gaussian(nech);
  db1->addColumns(vec2,"myvar2");
  db1->display();
  
  // Serialize db1
  (void) db1->dumpToNF("Neutral.Db.ascii",verbose);
  (void) db1->dumpToNF("Neutral2.Db.ascii",verbose);

  // Deserialize db2
  Db* db2 = Db::createFromNF("Neutral.Db.ascii",verbose);
  db2->display();
  delete db2;
  db2 = Db::createFromNF("Neutral2.Db.ascii",verbose);
  db2->display();
  (void) db2->dumpToNF("Neutral22.Db.ascii",verbose);

  // =======================
  // Checking Db (grid)
  // =======================

  // ===== Create the Grid Db
  DbGrid* dbg1 = DbGrid::create({12,10},{0.1,0.3},{0.2,0.4});
  vec1 = ut_vector_simulate_gaussian(dbg1->getSampleNumber());
  dbg1->addColumns(vec1,"myvar1",ELoc::Z, 0);
  vec2 = ut_vector_simulate_gaussian(dbg1->getSampleNumber());
  vec2[2] = TEST;
  vec2[5] = TEST;
  dbg1->addColumns(vec2,"myvar2",ELoc::Z, 1);
  dbg1->display();

  // Serialize dbg1
  (void) dbg1->dumpToNF("Neutral.Dbg.ascii",verbose);

  // Deserialize dbg2
  Db* dbg2 = DbGrid::createFromNF("Neutral.Dbg.ascii",verbose);
  dbg2->display();

  // =======================
  // Checking Polygons
  // =======================

  // ===== Create the Polygon poly1
  Polygons poly1;
  poly1.resetFromDb(db1);
  Polygons polyb;
  polyb.resetFromDb(dbg1);
  poly1.addPolySet(polyb.getPolySet(0));
  poly1.display();

  // Serialize poly1
  (void) poly1.dumpToNF("Neutral.Polygon.ascii",verbose);

  // Deserialize poly2
  Polygons* poly2 = Polygons::createFromNF("Neutral.Polygon.ascii",verbose);
  poly2->display();
  delete poly2;

  // =======================
  // Checking Vario
  // =======================

  // ===== Compute an experimental variogram
  VarioParam varioparam1;
  DirParam dirparam(2, 10, 0.02);
  varioparam1.addDirs(dirparam);
  Vario vario1 = Vario(&varioparam1,db1);
  vario1.computeByKey("vg");
  vario1.display();

  // Serialize vario1
  (void) vario1.dumpToNF("Neutral.Vario.ascii",verbose);

  // Deserialize vario2
  Vario* vario2 = Vario::createFromNF("Neutral.Vario.ascii",verbose);
  vario2->display();

  // =======================
  // Checking Model
  // =======================

  // ===== Create a Model
  db1->display();
  Model model1(db1);
  CovContext ctxt = model1.getContext();
  CovLMC covs(ctxt.getSpace());
  CovAniso cova(ECov::EXPONENTIAL, 0.3, 1., 0.2, ctxt);
  covs.addCov(&cova);
  model1.setCovList(&covs);
  model1.display();

  // Serialize model1
  (void) model1.dumpToNF("Neutral.Model.ascii",verbose);

  // Deserialize model2
  Model* model2 = Model::createFromNF("Neutral.Model.ascii",verbose);
  model2->display();

  // =======================
  // Checking Table
  // =======================

  // ===== Create a Table
  VectorVectorDouble table;
  int ncols = 3;
  int nrows = 10;
  table.resize(ncols);
  for (int icol = 0; icol < ncols; icol++)
    table[icol] = ut_vector_simulate_uniform(nrows);
  Table* table1 = Table::createFromArray(table);
  table1->display();

  // Serialize table
  (void) table1->dumpToNF("Neutral.Table.ascii",verbose);

  // Deserialize table1
  Table* table2 = Table::createFromNF("Neutral.Table.ascii",verbose);
  table2->display();

  // =======================
  // Checking Rule
  // =======================

  Rule* rule = Rule::createFromNames({"S","F1","T","F2","S","F3","F4"});
  rule->display();

  // Serialize
  (void) rule->dumpToNF("Neutral.Rule.ascii",verbose);

  // Deserialize
  Rule* rule2 = Rule::createFromNF("Neutral.Rule.ascii",verbose);
  rule2->display();

  delete db1;
  delete db2;
  delete dbg1;
  delete dbg2;
  delete vario2;
  delete model2;
  delete table1;
  delete table2;
  delete rule;
  delete rule2;
  return(0);
}
