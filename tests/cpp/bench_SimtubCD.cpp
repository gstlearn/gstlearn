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
/* This file is meant to demonstrate the process of using                     */
/* conditional simulations with Turning Band Method                           */
/*                                                                            */
/******************************************************************************/
#include "Enum/ECov.hpp"

#include "Variogram/Vario.hpp"
#include "Model/Model.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/File.hpp"
#include "Basic/Timer.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

/****************************************************************************/
/*!
** Main Program for bench marking the non-conditional simulation using Turning Bands
**
*****************************************************************************/
int main(int argc, char *argv[])
{
  Timer timer;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  int ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);
  bool verbose = false;

  // Creating the data base
  int ndat = 4;
  int nvar = 1;
  Db* data = Db::createFillRandom(ndat, ndim, nvar);
  
  // Creating a grid covering the same space
  VectorInt nx = { 200, 200 };
  VectorDouble dx = { 0.005, 0.005 };
  DbGrid* grid = DbGrid::create(nx, dx);

  // Creating a Model(s) for simulating a variable
  Model* model = Model::createFromParam(ECov::SPHERICAL,0.2);

  // Perform a non-conditional simulation on the Grid
  int nbsimu = 2;
  int nbtuba = 1000;

  if (verbose)
  {
    message("Conditional simulation(s) on grid using Turning Bands:\n");
    message("- Data (%d) samples\n", data->getNSample());
    message("- Grid (%d x %d)\n", grid->getNX(0), grid->getNX(1));
    message("- Number of Bands = %d\n", nbtuba);
    message("- Number of simulations = %d\n", nbsimu);
  }

  timer.reset();
  (void) simtub(data,grid,model,nullptr, nbsimu, 113423, nbtuba);
  timer.displayIntervalMilliseconds("Non-conditional Simulation on Grid", 6800);

  // Produce some statistics for comparison
  DbStringFormat* dbfmt = DbStringFormat::create(FLAG_STATS, {"Simu.*"});
  grid->display(dbfmt);
  delete dbfmt;

  delete data;
  delete grid;
  delete model;

  return (0);
}
