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
#include "Basic/OptCustom.hpp"
#include "Basic/OptDbg.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Neigh/NeighUnique.hpp"
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

  // Parameters
  double oldstyle = 1.;
  int ndat        = 3;
  int nvar        = 2;
  int nbsimu      = 4;
  int nbtuba      = 100;
  bool flagSK     = false;
  bool verbose    = true;
  OptCustom::define("oldStyle", oldstyle);

  // Creating the data base
  Db* data = Db::createFillRandom(ndat, ndim, nvar);
  
  // Creating a grid covering the same space
  VectorInt nx = { 20, 20 };
  VectorDouble dx = { 0.05, 0.05 };
  DbGrid* grid = DbGrid::create(nx, dx);

  // Create the Model
  int order = (flagSK) ? -1 : 0;
  Model* model =
    Model::createFillRandom(ndim, nvar, {ECov::EXPONENTIAL}, 1., order);

  // Creating the Neighborhood
  ANeigh* neigh = NeighUnique::create();

  // Define the verbose option
  if (verbose) OptDbg::setReference(1);

  message("Conditional simulation(s) on grid using Turning Bands:\n");
  message("- Data (%d) samples\n", data->getNSample());
  message("- Grid (%d x %d)\n", grid->getNX(0), grid->getNX(1));
  message("- Number of Bands = %d\n", nbtuba);
  message("- Number of simulations = %d\n", nbsimu);

  timer.reset();
  (void)simtub(data, grid, model, neigh, nbsimu, 113423, nbtuba);
  timer.displayIntervalMilliseconds("Non-conditional Simulation on Grid", 6800);

  // Produce some statistics for comparison
  DbStringFormat* dbfmt = DbStringFormat::create(FLAG_STATS, {"Simu.*"});
  grid->display(dbfmt);
  delete dbfmt;

  delete data;
  delete grid;
  delete model;
  delete neigh;

  return (0);
}
