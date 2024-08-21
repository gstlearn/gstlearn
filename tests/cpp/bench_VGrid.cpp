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
#include "Enum/ECalcVario.hpp"

#include "Variogram/Vario.hpp"
#include "Basic/Timer.hpp"
#include "Db/DbGrid.hpp"

/****************************************************************************/
/*!
** Main Program for bench marking the variogram calculation on the Grid
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

  // Creating a regular grid
  VectorInt nx = { 600, 400 };
  VectorDouble dx = { 1.0, 1.0 };
  DbGrid* grid    = DbGrid::create(nx, dx);
  VectorDouble tab = VH::simulateGaussian(grid->getSampleNumber());
  grid->addColumns(tab, "Var", ELoc::Z);
  if (verbose) grid->display();

  // ===============
  // On Grid samples
  // ===============

  int nlag = 10;
  int ndimax = 1;
  timer.reset();
  VarioParam* varioparamG =
    VarioParam::createMultipleFromGrid(grid, nlag, 0., VectorDouble(), nullptr, ndimax);
  Vario* varioG = Vario::computeFromDb(*varioparamG, grid, ECalcVario::VARIOGRAM);
  timer.displayIntervalMilliseconds("Variogram on Regular Grid", 1500);
  if (verbose && varioG != nullptr) varioG->display();

  delete grid;
  delete varioparamG;
  delete varioG;

  return (0);
}
