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
#include "Basic/AStringFormat.hpp"
#include "Basic/File.hpp"
#include "Basic/Law.hpp"
#include "Basic/PolyLine2D.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "LithoRule/Rule.hpp"
#include "Matrix/Table.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Polygon/Polygons.hpp"
#include "Variogram/Vario.hpp"
#include "Variogram/VarioParam.hpp"

using namespace gstlrn;
/****************************************************************************/
/*!
** Main Program for testing Serialize/Deserialize
** in Neutral or HDF5 format
**
*****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setPrefixName("test_Serialize-");

  // Next flag indicates if the format is NF (true) or H5 (false)
  bool flagNeutral = false;
  bool verbose     = false;
  Id mode          = 0;

  // =================
  // Preliminary tasks
  // =================

  Id nech           = 20;
  Db* db            = Db::createFromBox(nech, {0., 0.}, {1., 1.}, 32432);
  VectorDouble vec1 = VH::simulateGaussian(nech);
  db->addColumns(vec1, "myvar1", ELoc::Z, 0);
  VectorDouble vec2 = VH::simulateGaussian(nech);
  db->addColumns(vec2, "myvar2");

  DbGrid* dbg = DbGrid::create({12, 10}, {0.1, 0.3}, {0.2, 0.4});
  vec1        = VH::simulateGaussian(dbg->getNSample());
  dbg->addColumns(vec1, "myvar1", ELoc::Z, 0);
  vec2    = VH::simulateGaussian(dbg->getNSample());
  vec2[2] = TEST;
  vec2[5] = TEST;
  dbg->addColumns(vec2, "myvar2", ELoc::Z, 1);

  // ==
  // Db
  // ==

  if (mode == 0 || mode == 1)
  {
    Db* db1 = db->clone();
    db1->display();

    // Serialize
    if (flagNeutral)
      (void)db1->dumpToNF("Db.NF.ascii", EFormatNF::ASCII);
    else
      (void)db1->dumpToNF("Db.NF.h5", EFormatNF::H5);

    // Deserialize
    Db* db2 = nullptr;
    if (flagNeutral)
      db2 = Db::createFromNF("Db.NF.ascii", verbose);
    else
      db2 = Db::createFromNF("Db.NF.h5", verbose);
    db2->display();

    delete db1;
    delete db2;
  }

  // ======
  // DbGrid
  // ======

  if (mode == 0 || mode == 2)
  {
    DbGrid* dbg1 = dbg->clone();
    dbg1->display();

    // Serialize
    if (flagNeutral)
      (void)dbg1->dumpToNF("Dbg.NF.ascii", EFormatNF::ASCII);
    else
      (void)dbg1->dumpToNF("Dbg.NF.h5", EFormatNF::H5);

    // Deserialize
    Db* dbg2 = nullptr;
    if (flagNeutral)
      dbg2 = DbGrid::createFromNF("Dbg.NF.ascii", verbose);
    else
      dbg2 = DbGrid::createFromNF("Dbg.NF.h5", verbose);
    dbg2->display();

    delete dbg1;
    delete dbg2;
  }

  // ========
  // Polygons
  // ========

  if (mode == 0 || mode == 3)
  {
    auto* poly1 = new Polygons();
    poly1->resetFromDb(db);
    auto* polyb = new Polygons();
    polyb->resetFromDb(dbg);
    poly1->addPolyElem(polyb->getPolyElem(0));
    poly1->display();

    // Serialize
    if (flagNeutral)
      (void)poly1->dumpToNF("Polygon.NF.ascii", EFormatNF::ASCII);
    else
      (void)poly1->dumpToNF("Polygon.NF.h5", EFormatNF::H5);

    // Deserialize
    Polygons* poly2 = nullptr;
    if (flagNeutral)
      poly2 = Polygons::createFromNF("Polygon.NF.ascii", verbose);
    else
      poly2 = Polygons::createFromNF("Polygon.NF.h5", verbose);
    poly2->display();

    delete polyb;
    delete poly1;
    delete poly2;
  }

  // =====
  // Vario
  // =====

  if (mode == 0 || mode == 4)
  {
    VarioParam varioparam;
    DirParam dirparam(10, 0.02);
    varioparam.addDir(dirparam);
    Vario vario1(varioparam);
    vario1.compute(db, ECalcVario::VARIOGRAM);
    vario1.display();

    // Serialize
    if (flagNeutral)
      (void)vario1.dumpToNF("Vario.NF.ascii", EFormatNF::ASCII);
    else
      (void)vario1.dumpToNF("Vario.NF.h5", EFormatNF::H5);

    // Deserialize
    Vario* vario2 = nullptr;
    if (flagNeutral)
      vario2 = Vario::createFromNF("Vario.NF.ascii", verbose);
    else
      vario2 = Vario::createFromNF("Vario.NF.h5", verbose);
    vario2->display();

    delete vario2;
  }

  // ==========
  // Covariance
  // ==========

  if (mode == 0 || mode == 5)
  {
    VarioParam varioparam;
    DirParam dirparam(10, 0.02);
    varioparam.addDir(dirparam);
    Vario covariance1(varioparam);
    covariance1.compute(db, ECalcVario::COVARIANCE);
    covariance1.display();

    // Serialize
    if (flagNeutral)
      // (void)covariance1.dumpToNF("Covariance.NF.ascii", EFormatNF::ASCII);
      messerr("Useless trying to save a Covariance in old format");
    else
      (void)covariance1.dumpToNF("Covariance.NF.h5", EFormatNF::H5);

    // Deserialize
    Vario* covariance2 = nullptr;
    if (flagNeutral)
      // covariance2 = Vario::createFromNF("Covariance.NF.ascii", verbose);
      messerr("Useless trying to recover a Covariance in old format");
    else
      covariance2 = Vario::createFromNF("Covariance.NF.h5", verbose);
    covariance2->display();

    delete covariance2;
  }

  // =====
  // Model
  // =====

  if (mode == 0 || mode == 6)
  {
    Model* model1 = Model::createFromParam(ECov::EXPONENTIAL, 0.3, 0.2, 1.);
    model1->display();

    // Serialize model1
    if (flagNeutral)
      (void)model1->dumpToNF("Model.NF.ascii", EFormatNF::ASCII);
    else
      (void)model1->dumpToNF("Model.NF.h5", EFormatNF::H5);

    // Deserialize model2
    Model* model2 = nullptr;
    if (flagNeutral)
      model2 = Model::createFromNF("Model.NF.ascii", verbose);
    else
      model2 = Model::createFromNF("Model.NF.h5", verbose);
    model2->display();

    delete model1;
    delete model2;
  }

  // =====
  // Table
  // =====

  if (mode == 0 || mode == 7)
  {
    VectorVectorDouble table;
    Id ncols      = 3;
    Id nrows      = 10;
    Table* table1 = Table::create(nrows, ncols);
    for (Id irow = 0; irow < nrows; irow++)
      for (Id icol = 0; icol < ncols; icol++)
        table1->setValue(irow, icol, law_uniform());
    table1->display();

    // Serialize table
    if (flagNeutral)
      (void)table1->dumpToNF("Table.NF.ascii", EFormatNF::ASCII);
    else
      (void)table1->dumpToNF("Table.NF.h5", EFormatNF::H5);

    // Deserialize table1
    Table* table2 = nullptr;
    if (flagNeutral)
      table2 = Table::createFromNF("Table.NF.ascii", verbose);
    else
      table2 = Table::createFromNF("Table.NF.h5", verbose);
    table2->display();

    delete table1;
    delete table2;
  }

  // ====
  // Rule
  // ====

  if (mode == 0 || mode == 8)
  {
    Rule* rule1 = Rule::createFromNames({"S", "F1", "T", "F2", "S", "F3", "F4"});
    rule1->display();

    // Serialize
    if (flagNeutral)
      (void)rule1->dumpToNF("Rule.NF.ascii", EFormatNF::ASCII);
    else
      (void)rule1->dumpToNF("Rule.NF.h5", EFormatNF::H5);

    // Deserialize
    Rule* rule2 = nullptr;
    if (flagNeutral)
      rule2 = Rule::createFromNF("Rule.NF.ascii", verbose);
    else
      rule2 = Rule::createFromNF("Rule.NF.h5", verbose);
    rule2->display();

    delete rule1;
    delete rule2;
  }

  // ==========
  // PolyLine2D
  // ==========

  if (mode == 0 || mode == 9)
  {
    Id npolyline           = 100;
    VectorDouble xpolyline = VH::simulateGaussian(npolyline);
    VectorDouble ypolyline = VH::simulateGaussian(npolyline);
    auto* polyline1        = new PolyLine2D(xpolyline, ypolyline);
    AStringFormat afmt(3);
    polyline1->display(&afmt);

    // Serialize
    if (flagNeutral)
      (void)polyline1->dumpToNF("Polyline.NF.ascii", EFormatNF::ASCII);
    else
      (void)polyline1->dumpToNF("Polyline.NF.h5", EFormatNF::H5);

    // Deserialize
    PolyLine2D* polyline2 = nullptr;
    if (flagNeutral)
      polyline2 = PolyLine2D::createFromNF("Polyline.NF.ascii", verbose);
    else
      polyline2 = PolyLine2D::createFromNF("Polyline.NF.h5", verbose);
    polyline2->display(&afmt);

    delete polyline1;
    delete polyline2;
  }

  // ===================
  // Moving Neighborhood
  // ===================

  if (mode == 0 || mode == 10)
  {
    Id nmaxi            = 20;
    double radius       = 4.;
    Id nmini            = 2;
    Id nsect            = 5;
    Id nsmax            = 3;
    VectorDouble coeffs = {2., 3.};
    VectorDouble angles = {25., 0.};
    bool useBallTree    = true;

    NeighMoving* neigh1 = NeighMoving::create(false, nmaxi, radius,
                                              nmini, nsect, nsmax,
                                              coeffs, angles, useBallTree);
    neigh1->display();

    // Serialize
    if (flagNeutral)
      (void)neigh1->dumpToNF("NeighMoving.NF.ascii", EFormatNF::ASCII);
    else
      (void)neigh1->dumpToNF("NeighMoving.NF.h5", EFormatNF::H5);

    // Deserialize
    NeighMoving* neigh2 = nullptr;
    if (flagNeutral)
      neigh2 = NeighMoving::createFromNF("NeighMoving.NF.ascii", verbose);
    else
      neigh2 = NeighMoving::createFromNF("NeighMoving.NF.h5", verbose);
    neigh2->display();

    delete neigh1;
    delete neigh2;
  }

  // ===================
  // Unique Neighborhood
  // ===================

  if (mode == 0 || mode == 11)
  {
    NeighUnique* neigh1 = NeighUnique::create();
    neigh1->display();

    // Serialize
    if (flagNeutral)
      (void)neigh1->dumpToNF("NeighUnique.NF.ascii", EFormatNF::ASCII);
    else
      (void)neigh1->dumpToNF("NeighUnique.NF.h5", EFormatNF::H5);

    // Deserialize
    NeighUnique* neigh2 = nullptr;
    if (flagNeutral)
      neigh2 = NeighUnique::createFromNF("NeighUnique.NF.ascii", verbose);
    else
      neigh2 = NeighUnique::createFromNF("NeighUnique.NF.h5", verbose);
    neigh2->display();

    delete neigh1;
    delete neigh2;
  }

  // ===============
  // Meshing (Turbo)
  // ===============

  if (mode == 0 || mode == 12)
  {
    MeshETurbo* mesh1 = MeshETurbo::create({10, 10});
    mesh1->display();

    // Serialize
    if (flagNeutral)
      (void)mesh1->dumpToNF("Mesh.NF.ascii", EFormatNF::ASCII);
    else
      (void)mesh1->dumpToNF("Mesh.NF.h5", EFormatNF::H5);

    // Deserialize
    MeshETurbo* mesh2 = nullptr;
    if (flagNeutral)
      mesh2 = MeshETurbo::createFromNF("Mesh.NF.ascii", verbose);
    else
      mesh2 = MeshETurbo::createFromNF("Mesh.NF.h5", verbose);
    mesh2->display();

    delete mesh1;
    delete mesh2;
  }

  // Cleaning procedure
  delete db;
  delete dbg;

  return (0);
}
