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
#include "Anamorphosis/AnamHermite.hpp"
#include "Anamorphosis/CalcAnamTransform.hpp"
#include "Basic/File.hpp"
#include "Basic/Message.hpp"
#include "Basic/OptCst.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/VectorNumT.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Enum/ECst.hpp"
#include "Estimation/Estimations.hpp"
#include "Estimation/Vecchia.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighImage.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Simulation/Simulations.hpp"
#include "Tree/Ball.hpp"

#include "geoslib_f.h"

using namespace gstlrn;

static Db* _createDb(Id nvar = 1, Id ndat = 5, Id seed = 4243)
{
  Id ndim = 2;
  Db* db  = Db::createFillRandom(ndat, ndim, nvar, 0, 0, 0., 0.,
                                 VectorDouble(), VectorDouble(), VectorDouble(), seed);
  VectorDouble sel(ndat, 1.);
  db->addSelection(sel, "sel");
  return db;
}

static Model* _createModel(Id nvar = 1, double range = 0.2, bool verbose = false)
{
  MatrixSymmetric* sills = MatrixSymmetric::createRandomDefinitePositive(nvar);
  Model* model           = Model::createFromParam(ECov::EXPONENTIAL, range, TEST, 1., VectorDouble(), *sills);
  delete sills;
  if (verbose) model->display();
  return model;
}

static DbGrid* _createGrid(Id nx = 2)
{
  DbGrid* grid = DbGrid::create({nx, nx}, {1. / nx, 1. / nx});
  return grid;
}

static void _dumpLimit(Id mode)
{
  Id limit = (mode > 0) ? -1 : 7;
  OptCst::define(ECst::NTCOL, limit);
  OptCst::define(ECst::NTROW, limit);
}

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_Vecchia-");

  // Global parameters
  Id mode      = 0;
  Id nbVecchia = 3;
  bool verbose = true;

  if (mode == 0 || mode == 1)
  {
    mestitle(0, "Checking Vecchia Class (based on 2 Dbs with NA)");
    _dumpLimit(1);
    auto dbfmt = DbStringFormat(FLAG_VARS | FLAG_ARRAY, VectorString(), VectorInt(), false);
    int indNa;
    bool addSelection = false;

    indNa   = 3;
    Db* db1 = _createDb(1, 5, 3261);
    db1->setValue("z", indNa, TEST);
    if (addSelection) db1->setValue("sel", indNa, 0.);
    db1->display(&dbfmt);

    indNa   = 2;
    Db* db2 = _createDb(1, 6, 4204);
    db2->setValue("z", indNa, TEST);
    if (addSelection) db2->setValue("sel", indNa, 0.);
    db2->display(&dbfmt);

    Model* model = _createModel(1);
    Vecchia V(model, nbVecchia, db1, db2);
    auto Ranks = findNN(db1, db2, nbVecchia + 1, false, verbose);
    (void)V.computeLower(Ranks, verbose);
    delete db1;
    delete db2;
    delete model;
    _dumpLimit(-1);
  }

  if (mode == 0 || mode == 2)
  {
    mestitle(0, "Kriging with Vecchia approximation");
    _dumpLimit(1);
    auto dbfmt   = DbStringFormat(FLAG_STATS, {"Vecchia*"});
    Db* db       = _createDb(1, 5);
    Model* model = _createModel(1);
    DbGrid* grid = _createGrid();
    krigingVecchia(db, grid, model, nbVecchia, true);
    grid->display(&dbfmt);
    delete db;
    delete model;
    delete grid;
    _dumpLimit(-1);
  }

  if (mode == 0 || mode == 3)
  {
    mestitle(0, "Log-Likelihood");
    Db* db              = _createDb(1, 20000);
    Model* model        = _createModel(1);
    const double result = logLikelihoodVecchia(db, model, nbVecchia, false);
    message("Log-likelihood = %f\n", result);
    delete db;
    delete model;
  }

  if (mode == 0 || mode == 4)
  {
    nbVecchia = 4;
    mestitle(0, "Kriging with Vecchia approximation (nvar=2)");
    auto dbfmt   = DbStringFormat(FLAG_STATS, {"Vecchia*"});
    Db* db       = _createDb(2, 10);
    Model* model = _createModel(2);
    DbGrid* grid = _createGrid(100);
    krigingVecchia(db, grid, model, nbVecchia, false);
    grid->display(&dbfmt);
    (void)db->dumpToNF("Data.dat");
    (void)grid->dumpToNF("Grid.dat");
    delete db;
    delete model;
    delete grid;
  }

  if (mode == 0 || mode == 5)
  {
    int ndat = 6;
    mestitle(0, "Comparing Log-Likelihood with Vecchia approximation when all neighbors are used");
    Db* db = Db::createFillRandom(ndat, 2, 1, 0, 0, 1., 0., 0., {0., 0.}, {1., 1.}, 1234, false);
    db->setColumn(VectorDouble(ndat, 0.01), "verr");
    db->setLocator("verr", ELoc::V);
    auto ran = VH::sequenceVD(1, 1 + ndat - 1, 1);
    db->setColumn(ran, "ranges");
    // Add duplicate
    db->setSampleCoordinates(1, {db->getCoordinate(0, 0), db->getCoordinate(0, 1)});
    Model* model = _createModel(1);
    model->getCovAniso(0)->attachNoStatDb(db);
    model->getCovAniso(0)->makeRangeNoStatDb("ranges", 0.);
    double result = logLikelihoodVecchia(db, model, ndat, false);
    message("Log-likelihood with    Vecchia= %f\n", result);
    result = model->computeLogLikelihood(db);
    message("Log-likelihood without Vecchia= %f\n", result);
    delete db;
    delete model;
  }
  return (0);
}
