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
  int ndat        = 4;
  int nvar        = 2;
  int nbsimu      = 2;
  int nbtuba      = 1000;
  double oldstyle = 1.;
  bool flagSK     = false;
  OptCustom::define("oldStyle", oldstyle);

  // Creating the data base
  Db* data = Db::createFillRandom(ndat, ndim, nvar);
  
  // Creating a grid covering the same space
  VectorInt nx = { 20, 20 };
  VectorDouble dx = { 0.05, 0.05 };
  DbGrid* grid = DbGrid::create(nx, dx);

  // Creating a Model(s) for simulating a variable
  double scale = 0.7;
  MatrixSquareSymmetric* sills =
    MatrixSquareSymmetric::createRandomDefinitePositive(nvar);
  Model* model = Model::createFromParam(ECov::EXPONENTIAL, scale, 0., 0., VectorDouble(),
                                        *sills, VectorDouble(), nullptr, false);
  if (flagSK)
    model->setMeans(VH::simulateGaussian(nvar));
  else
    model->setDriftIRF(0, 0);

  // Creating the Neighborhood
  ANeigh* neigh = NeighUnique::create();

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
