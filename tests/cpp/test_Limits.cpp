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
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "API/SPDE.hpp"
#include "Model/Model.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program is meant to check the manipulation of the Db
 **
 *****************************************************************************/
int main(int argc, char *argv[])
{
  VectorString names1;
  VectorString names2;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  DbStringFormat dbfmt(FLAG_STATS);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("TestLimits-");
  int seed = 10355;
  law_set_random_seed(seed);

  // Creating the Grid Rotated Db
  DbGrid* grid = DbGrid::create({30,20}, {1.,2.}, {10.,20.});

  // Creating the Model
  Model model(1, 2);
  model.addCovFromParam(ECov::CUBIC, 0., 2., 1., {10., 45.},
                        MatrixSquareSymmetric(), {30.,0.});
  model.display();

  // Simulating a variable on the grid
  (void) simtub(nullptr, grid, &model, nullptr);
  dbfmt = DbStringFormat(FLAG_STATS, {"Simu"});
  grid->display(&dbfmt);

  // Creating a set of Limits
  Limits limits({-1.,-0.5,0.,0.5,1.});
  limits.display();

  // Other option
  grid->setLocator("Simu", ELoc::Z, 0);
  limits.toIndicator(grid,"Simu",0);
  dbfmt = DbStringFormat(FLAG_ARRAY, {"Simu", "Indicator.Simu.Mean"});
  grid->display(&dbfmt);

  // Convert into Indicators
  grid->setLocator("Simu", ELoc::Z, 0);
  limits.toIndicator(grid,"Simu",1);
  dbfmt = DbStringFormat(FLAG_ARRAY, {"Indicator.Simu.Class*"});
  grid->display(&dbfmt);

  delete grid;
  return 0;
}

