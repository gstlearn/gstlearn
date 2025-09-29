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
#include "Basic/AStringable.hpp"
#include "Basic/File.hpp"
#include "Basic/OptCst.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Estimation/CalcGlobal.hpp"
#include "Estimation/CalcImage.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Estimation/Vecchia.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighImage.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Tree/Ball.hpp"

#include "geoslib_f.h"

using namespace gstlrn;

static Db* _createDb(Id nvar = 1, Id ndat = 5, bool verbose = false)
{
  Id ndim = 2;
  Db* db  = Db::createFillRandom(ndat, ndim, nvar);
  // db->setZVariable(3, 0, TEST);
  DbStringFormat dbfmt1(FLAG_ARRAY);
  if (verbose) db->display(&dbfmt1);
  return db;
}

static Model* _createModel(Id nvar = 1, double range = 0.5, bool verbose = false)
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

  // Global parameters
  Id mode       = 0;
  Id nb_vecchia = 3;
  DbStringFormat dbfmt(FLAG_STATS, {"Vecchia*"});

  if (mode == 0 || mode == 1)
  {
    mestitle(0, "Checking Vecchia Class");
    _dumpLimit(1);
    Db* db       = _createDb(1, 5, true);
    Model* model = _createModel(1);
    Vecchia V(model, nb_vecchia, db);
    auto Ranks = findNN(db, nullptr, nb_vecchia + 1, false, true);
    (void)V.computeLower(Ranks, true);
    delete db;
    delete model;
    _dumpLimit(-1);
  }

  if (mode == 0 || mode == 2)
  {
    mestitle(0, "Kriging with Vecchia approximation");
    _dumpLimit(1);
    Db* db       = _createDb(1, 5, false);
    Model* model = _createModel(1);
    DbGrid* grid = _createGrid();
    krigingVecchia(db, grid, model, nb_vecchia, true);
    grid->display(&dbfmt);
    delete db;
    delete model;
    delete grid;
    _dumpLimit(-1);
  }

  if (mode == 0 || mode == 3)
  {
    mestitle(0, "Log-Likelihood");
    Db* db              = _createDb(1, 20000, false);
    Model* model        = _createModel(1);
    const double result = logLikelihoodVecchia(db, model, nb_vecchia, false);
    message("Log-likelihood = %f\n", result);
    delete db;
    delete model;
  }

  if (mode == 0 || mode == 4)
  {
    nb_vecchia = 2;
    mestitle(0, "Kriging with Vecchia approximation (nvar=2)");
    Db* db       = _createDb(2, 10, false);
    Model* model = _createModel(2);
    DbGrid* grid = _createGrid(100);
    krigingVecchia(db, grid, model, nb_vecchia, false);
    grid->display(&dbfmt);
    delete db;
    delete model;
    delete grid;
  }

  return (0);
}
