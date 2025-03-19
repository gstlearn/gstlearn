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
#include "geoslib_f.h"

#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Basic/File.hpp"
#include "Basic/OptCst.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighImage.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Anamorphosis/CalcAnamTransform.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Estimation/CalcImage.hpp"
#include "Estimation/CalcGlobal.hpp"
#include "Estimation/Vecchia.hpp"
#include "Matrix/MatrixT.hpp"
#include "Tree/Ball.hpp"

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
  int nb_neigh = 3;
  int mode     = 0;
  bool verbose = false;
  OptCst::define(ECst::NTCOL, -1);
  OptCst::define(ECst::NTROW, -1);

  int ndat = 4;
  Db* db = Db::createFillRandom(ndat, 2, 1);

  double range = 0.5;
  Model* model = Model::createFromParam(ECov::EXPONENTIAL, range);

  int nx = 2;
  DbGrid* grid = DbGrid::create({nx, nx}, {1. / nx, 1. / nx});

  if (mode == 0 || mode == 1)
  {
    mestitle(0, "Checking Vecchia Class");
    verbose = true;
    Vecchia V          = Vecchia(model, db);
    MatrixT<int> Ranks = findNN(db, nullptr, nb_neigh+1, false, verbose);
    (void)V.computeLower(Ranks, verbose);
  }

  if (mode == 0 || mode == 2)
  {
    mestitle(0, "Kriging with Vecchia approximation");
    krigingVecchia(db, grid, model, nb_neigh);

    // Get some statistics for check printout
    DbStringFormat dbfmt(FLAG_STATS, {"Vecchia*"});
    grid->display(&dbfmt);
  }

  // ====================== Free pointers ==================================
  delete db;
  delete grid;
  delete model;

  return (0);
}
