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
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Variogram/VarioParam.hpp"
#include "Variogram/Vario.hpp"
#include "Matrix/Table.hpp"
#include "Model/Model.hpp"
#include "Basic/File.hpp"
#include "Basic/Law.hpp"
#include "Basic/PolyLine2D.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/AStringFormat.hpp"
#include "Polygon/Polygons.hpp"
#include "LithoRule/Rule.hpp"

/****************************************************************************/
/*!
** Main Program for testing the sparse matrix algebra
**
*****************************************************************************/
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

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
  vec1 = VH::simulateGaussian(dbg1->getNSample());
  dbg1->addColumns(vec1,"myvar1",ELoc::Z, 0);
  vec2 = VH::simulateGaussian(dbg1->getNSample());
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
  poly1.addPolyElem(polyb.getPolyElem(0));
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
  Vario vario1 = Vario(varioparam1);
  vario1.compute(db1, ECalcVario::VARIOGRAM);
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
  Table* table1 = Table::create(nrows, ncols);
  for (int irow = 0; irow < nrows; irow++)
    for (int icol = 0; icol < ncols; icol++)
      table1->setValue(irow, icol, law_uniform());
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
  PolyLine2D* polyline   = new PolyLine2D(xpolyline, ypolyline);
  AStringFormat afmt(3);
  polyline->display(&afmt);

  // Serialize
  (void) polyline->dumpToNF("Neutral.Polyline.ascii");

  // Deserialize
  PolyLine2D* polyline2 =
    PolyLine2D::createFromNF("Neutral.Polyline.ascii", verbose);
  polyline2->display(&afmt);

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
