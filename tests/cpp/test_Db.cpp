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
#include "Enum/ECst.hpp"

#include "API/SPDE.hpp"
#include "Basic/File.hpp"
#include "Basic/OptCst.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"

using namespace gstlrn;
/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program is meant to check the manipulation of the Db
 **
 *****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setPrefixName("test_Db-");
  int mode = 0;

  ////////////////////////////
  // Testing the selections //
  ////////////////////////////

  if (mode == 0 || mode == 1)
  {
    mestitle(0, "Testing selection management");
    // Creating the Grid Rotated Db
    DbGrid* grid = DbGrid::create({6, 4}, {1., 2.}, {10., 20.}, {10., 0.});
    auto nech    = grid->getNSample();

    // Creating the Model
    Model* model = Model::createFromParam(ECov::CUBIC, 0., 1., 1., {10., 45.},
                                          MatrixSymmetric(), {30., 0.});

    // First selection generated with Bernoulli (proba=0.6)
    VectorDouble sel1 = VH::simulateBernoulli(nech, 0.6);
    print_vector("sel1", sel1, true, true);
    grid->addSelection(sel1, "Sel1");

    // Second selection generated with Bernoulli (proba=0.4) combined with previous one
    VectorDouble sel2 = VH::simulateBernoulli(nech, 0.4);
    print_vector("sel2", sel2, true, true);
    grid->addSelection(sel2, "Sel2", "and");

    // Retrieve resulting selection for check
    VectorDouble sel3 = grid->getSelections();
    print_vector("sel1 && sel2", sel3, true, true);

    // Testing Filters on Db printout (only Statistics on the variables "Sel*")
    DbStringFormat dbfmt(FLAG_VARS | FLAG_STATS, {"Sel*"});
    grid->display(&dbfmt);

    // Creating a Selection by setting individual values
    OptCst::define(ECst::NTROW, -1);
    DbStringFormat dbfmt2(FLAG_VARS | FLAG_ARRAY);
    grid->addSelection(VectorDouble(), "mySel");
    grid->setValue("mySel", 12, 0.);
    grid->setValue("mySel", 19, 0.);
    grid->display(&dbfmt2);

    delete grid;
    delete model;
  }

  //////////////////////////////////
  // Testing getRanks() functions //
  //////////////////////////////////

  if (mode == 0 || mode == 2)
  {
    mestitle(0, "Testing Db::getRanks functions");
    Id ndim    = 2;
    Id ndat    = 15;
    auto* db   = Db::createFillRandom(ndat, ndim, 0, 0, 0, 0., 0.3);
    auto dbfmt = DbStringFormat(FLAG_ARRAY, VectorString(), VectorInt(), false);
    db->display(&dbfmt);

    print_vector("Vector of Relative ranks of active samples",
                 db->getRankAbsoluteToRelativeVec(), true, true);

    print_vector("Vector of Absolute ranks of active samples",
                 db->getRankRelativeToAbsoluteVec(), true, true);
    delete db;
  }

  ///////////////////////
  // Testing db_reduce //
  ///////////////////////

  if (mode == 0 || mode == 3)
  {
    mestitle(0, "Testing db_reduce facility");
    Id nvar    = 3;
    Id ndim    = 2;
    Id ndat    = 15;
    auto dbfmt = DbStringFormat(FLAG_ARRAY, VectorString(), VectorInt(), false);
    VectorDouble hetero(nvar, 0.1);
    auto* db = Db::createFillRandom(ndat, ndim, nvar, 0, 0, 0., 0.1, hetero);
    db->display(&dbfmt);
    VectorString names = {"z-1", "z-3"};
    VectorInt ranks    = VH::sequence(5., 3, 1);
    Db* dbaux;

    message("\n---> Reducing by:\n");
    message(" - selecting some variables (z-1 and z-3)\n");
    message(" - suppressing masked samples\n");
    dbaux = Db::createReduce(db, names);
    dbaux->display(&dbfmt);
    delete dbaux;

    message("\n---> Reducing by:\n");
    message(" - selecting some variables (z-1 and z-3)\n");
    message(" - selecting some sample ranks (5 samples starting from rank 3)\n");
    dbaux = Db::createReduce(db, names, ranks);
    dbaux->display(&dbfmt);
    delete dbaux;

    message("\n---> Reducing by:\n");
    message(" - selecting some variables (z-1 and z-3)\n");
    message(" - selecting only samples where all variables are defined\n");
    dbaux = Db::createReduce(db, names, VectorInt(), true);
    dbaux->display(&dbfmt);
    delete dbaux;

    delete db;
  }

  ///////////////////////////////
  // Testing operator overload //
  ///////////////////////////////
  if (mode == 0 || mode == 4)
  {
    mestitle(0, "Testing operator overload");
    Id nvar    = 3;
    Id ndim    = 2;
    Id ndat    = 15;
    auto dbfmt = DbStringFormat(FLAG_ARRAY, VectorString(), VectorInt(), false);
    VectorDouble hetero(nvar, 0.1);
    auto* db           = Db::createFillRandom(ndat, ndim, nvar, 0, 0, 0., 0.1, hetero);
    VectorString names = {"z-1", "z-3"};
    db->display(&dbfmt);
    double value = (*db)(3, "z-1");
    message("Initial Value = %f\n", value);
    (*db)(3, "z-1") = 12.;
    db->display(&dbfmt);

    // Testing invalid answers
    value = (*db)(3000, "z-1");
    messerr("This should be an error = %f", value);
    value = (*db)(3, "stupid");
    messerr("This should be an error = %f", value);

    (*db)(30000, "z-1") = 12.;
    messerr("Previous statement is an error");
    (*db)(3, "stupid") = 12.;
    messerr("Previous statement is an error");
    delete db;
  }
  return 0;
}
