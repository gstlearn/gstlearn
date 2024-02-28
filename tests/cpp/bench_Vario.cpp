/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/*                                                                            */
/* This file is meant to establish Bench Mark for Variogram calculations      */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f.h"

#include "Enum/ECalcVario.hpp"
#include "Enum/ECov.hpp"

#include "Variogram/Vario.hpp"
#include "Variogram/VMap.hpp"
#include "Model/Model.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/File.hpp"
#include "Basic/Timer.hpp"
#include "Db/Db.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

/****************************************************************************/
/*!
** Main Program for bench marking the variogram calculation
** - on a set of isolated points (including code, Faulting, ...)
** - on a regular grid
**
*****************************************************************************/
int main(int argc, char *argv[])
{
  Timer timer;
  bool verbose = true;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  int ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Creating a Point Data base in the 1x1 square with 'nech' samples
  int nech = 2000;
  int nvar = 1;
  int ncode = 5;
  int seed = 143421;
  Db *db = Db::createFillRandom(nech, ndim, nvar, 0, ncode, 0., 0.,
                                VectorDouble(), VectorDouble(), VectorDouble(), seed, 1);

  // Creating a grid covering the same space
  VectorInt nx = { 200, 200 };
  VectorDouble dx = { 0.005, 0.005 };
  DbGrid* grid = DbGrid::create(nx, dx);

  // Creating a Model(s) for simulating a variable
  Model* model = Model::createFromParam(ECov::BESSEL_K,0.2);

  // Perform a non-conditional simulation on the Db and on the Grid
  (void) simtub(nullptr,db,model);
  (void) simtub(nullptr,grid,model);

  // Defining a Fault system
  Faults* faults = new Faults();
  PolyLine2D polyline({0.1, 0.6, 0.7}, {0.1, 0.4, 0.9});
  faults->addFault(polyline);

  // ===============
  // On Data samples
  // ===============

  mestitle(1, "Experimental variogram on Data Samples");
  int ndir = 4;
  int nlag = 20;
  message("- on a Db containing %d samples\n", nech);
  message("- for a variogram calculated in %d directions with %d lags\n", ndir, nlag);

  timer.reset();
  VarioParam* varioparamP = VarioParam::createMultiple(ndir, nlag, 0.5 / nlag);
  Vario* varioP = Vario::computeFromDb(*varioparamP,db,ECalcVario::VARIOGRAM);
  timer.displayIntervalMilliseconds("Variogram on Isolated Points", 2600);
  if (verbose) varioP->display();

  // =====================================
  // On Data samples (with all attributes)
  // =====================================

  mestitle(1, "Experimental variogram on Data Samples (with attributes)");
  message("- on a Db containing %d samples (with code attribute)\n", nech);
  message("- taking a Fault System into account\n");
  message("- for a variogram calculated in 1 direction with %d lags\n", nlag);

  timer.reset();
  double dlag = 0.05;
  double toldis = 0.5;
  double tolang = 45.;
  int optcode = 1;
  double tolcode = 2;
  double angle2D = 30.;
  DirParam dirparam = DirParam(nlag, dlag, toldis, tolang, optcode, 0, TEST, TEST, tolcode,
                                VectorDouble(), VectorDouble(), angle2D);
  VarioParam varioparamC = VarioParam();
  varioparamC.addDir(dirparam);
  varioparamC.addFaults(faults);
  Vario* varioC = Vario::computeFromDb(varioparamC,db,ECalcVario::VARIOGRAM);
  timer.displayIntervalMilliseconds("Variogram on Isolated Points (with attributes)", 1100);
  if (verbose) varioC->display();

  // ===============
  // On Grid samples
  // ===============

  mestitle(1, "Experimental variogram on Grid");
  message("- on a grid of %d by %d pixels\n",nx[0],nx[1]);
  message("- for a variogram calculated along main directions with %d lags\n",nlag);

  timer.reset();
  VarioParam* varioparamG = VarioParam::createMultipleFromGrid(grid, nlag);
  Vario* varioG = Vario::computeFromDb(*varioparamG, grid, ECalcVario::VARIOGRAM);
  timer.displayIntervalMilliseconds("Variogram on Regular Grid", 1500);
  if (verbose) varioG->display();

  // ==========================================
  // Calculating Variogram Map on Isolated Data
  // ==========================================

  mestitle(1, "Variogram Map on Isolated Data");
  int ncell = 50;
  message("- on a Db containing %d samples\n", nech);
  message("- for %d by %d cells (automatic dimensions)\n", ncell, ncell);

  timer.reset();
  Db* vmapP = db_vmap(db, ECalcVario::VARIOGRAM, {ncell,ncell});
  timer.displayIntervalMilliseconds("Variogram Map on Isolated Points", 2400);

  // =================================
  // Calculating Variogram Map on Grid
  // =================================

  mestitle(1, "Variogram Map on Grid");
  timer.reset();
  Db* vmapG = db_vmap(grid, ECalcVario::VARIOGRAM, {100,100});
  timer.displayIntervalMilliseconds("Variogram Map on Regular Grid", 100);

  delete db;
  delete grid;
  delete model;
  delete varioparamP;
  delete varioparamG;
  delete varioP;
  delete varioC;
  delete varioG;
  delete vmapP;
  delete vmapG;
  delete faults;

  return (0);
}
