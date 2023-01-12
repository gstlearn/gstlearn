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
#include "geoslib_old_f.h"

#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Variogram/VarioParam.hpp"
#include "Variogram/Vario.hpp"
#include "Model/Model.hpp"
#include "Basic/Table.hpp"
#include "Basic/File.hpp"
#include "Basic/PolyLine2D.hpp"
#include "Basic/VectorHelper.hpp"
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
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("TS-");

  // =======================
  // Checking Db
  // =======================

  // ===== Create the Db db1
  int nech = 20;
  bool verbose = false;

  Db* db1 = Db::createFromBox(nech,{0.,0.},{1.,1.}, 32432);
  VectorDouble vec1 = VH::simulateGaussian(nech);
  db1->addColumns(vec1,"myvar1",ELoc::Z, 0);
  VectorDouble vec2 = VH::simulateGaussian(nech);
  db1->addColumns(vec2,"myvar2");
  db1->display();
  
  // Serialize db1
  (void) db1->dumpToNF("Neutral.Db.ascii");
  (void) db1->dumpToNF("Neutral2.Db.ascii");

  // Deserialize db2
  Db* db2 = Db::createFromNF("Neutral.Db.ascii",verbose);
  db2->display();
  delete db2;
  db2 = Db::createFromNF("Neutral2.Db.ascii",verbose);
  db2->display();
  (void) db2->dumpToNF("Neutral22.Db.ascii");

  // =======================
  // Checking Db (grid)
  // =======================

  // ===== Create the Grid Db
  DbGrid* dbg1 = DbGrid::create({12,10},{0.1,0.3},{0.2,0.4});
  vec1 = VH::simulateGaussian(dbg1->getSampleNumber());
  dbg1->addColumns(vec1,"myvar1",ELoc::Z, 0);
  vec2 = VH::simulateGaussian(dbg1->getSampleNumber());
  vec2[2] = TEST;
  vec2[5] = TEST;
  dbg1->addColumns(vec2,"myvar2",ELoc::Z, 1);
  dbg1->display();

  // Serialize dbg1
  (void) dbg1->dumpToNF("Neutral.Dbg.ascii");

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
  (void) poly1.dumpToNF("Neutral.Polygon.ascii");

  // Deserialize poly2
  Polygons* poly2 = Polygons::createFromNF("Neutral.Polygon.ascii",verbose);
  poly2->display();
  delete poly2;

  // =======================
  // Checking Vario
  // =======================

  // ===== Compute an experimental variogram
  VarioParam varioparam1;
  DirParam dirparam(10, 0.02);
  varioparam1.addDir(dirparam);
  Vario vario1 = Vario(&varioparam1,db1);
  vario1.computeByKey("vg");
  vario1.display();

  // Serialize vario1
  (void) vario1.dumpToNF("Neutral.Vario.ascii");

  // Deserialize vario2
  Vario* vario2 = Vario::createFromNF("Neutral.Vario.ascii",verbose);
  vario2->display();

  // =======================
  // Checking Model
  // =======================

  // ===== Create a Model
  Model* model1 = Model::createFromParam(ECov::EXPONENTIAL, 0.3, 0.2, 1.);
  model1->display();

  // Serialize model1
  (void) model1->dumpToNF("Neutral.Model.ascii");

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
    table[icol] = VH::simulateUniform(nrows);
  Table* table1 = Table::createFromVVD(table, false);
  table1->display();

  // Serialize table
  (void) table1->dumpToNF("Neutral.Table.ascii");

  // Deserialize table1
  Table* table2 = Table::createFromNF("Neutral.Table.ascii",verbose);
  table2->display();

  // =======================
  // Checking Rule
  // =======================

  Rule* rule = Rule::createFromNames({"S","F1","T","F2","S","F3","F4"});
  rule->display();

  // Serialize
  (void) rule->dumpToNF("Neutral.Rule.ascii");

  // Deserialize
  Rule* rule2 = Rule::createFromNF("Neutral.Rule.ascii",verbose);
  rule2->display();

  // ======================
  // Checking PolyLine2D
  // ======================

  int npolyline = 100;
  VectorDouble xpolyline = VH::simulateGaussian(npolyline);
  VectorDouble ypolyline = VH::simulateGaussian(npolyline);
  PolyLine2D* polyline = new PolyLine2D(xpolyline, ypolyline);
  polyline->display();

  // Serialize
  (void) polyline->dumpToNF("Neutral.Polyline.ascii");

  // Deserialize
  PolyLine2D* polyline2 = PolyLine2D::createFromNF("Neutral.Polyline.ascii", verbose);
  polyline2->display();

  delete db1;
  delete db2;
  delete dbg1;
  delete dbg2;
  delete vario2;
  delete model1;
  delete model2;
  delete table1;
  delete table2;
  delete rule;
  delete rule2;
  delete polyline;
  delete polyline2;

  return(0);
}
