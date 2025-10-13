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

// This test is meant to demonstrate the Eden Simulation
// which also demonstrates the SKIN methodology

#include "Enum/ECov.hpp"
#include "Enum/ESpaceType.hpp"

#include "Basic/File.hpp"
#include "Basic/Law.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Simulation/CalcSimuEden.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Space/ASpaceObject.hpp"

using namespace gstlrn;

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

  ASerializable::setPrefixName("test_Eden-");

  // Global parameters
  Id ndim = 2;
  law_set_random_seed(32131);

  defineDefaultSpace(ESpaceType::RN, ndim);
  DbStringFormat dbfmt(FLAG_STATS);

  // Generate the output grid
  VectorInt nx = {301, 101};
  DbGrid* grid = DbGrid::create(nx);

  // Simulate continuous variable
  double range = 30.;
  Model* model = Model::createFromParam(ECov::CUBIC, range);
  (void)simtub(nullptr, grid, model);

  // Operate a transformation to convert into 3 (nested) facies
  Id nfacies         = 3;
  double thresh      = 1.;
  VectorDouble vmini = {TEST, -thresh, thresh};
  VectorDouble vmaxi = {-thresh, thresh, TEST};

  Limits* limits = Limits::create(vmini, vmaxi);
  limits->toCategory(grid, "Simu");
  grid->setName("Category*", "Facies");

  // Add a Fluid information
  Id nfluids = 1;
  grid->addColumnsByConstant(nfluids, TEST, "Fluid");
  (void)grid->assignGridColumn("Fluid", 0, 100, 1.);

  // Fluid propagation
  // Speed: (dir + 6 * (facies * nfluids + fluid))
  // Direction = 0: +X; 1: -X; 2: +Y; 3: -Y; 4: +Z(up); 5: -Z(down)

  Id sl            = 1;
  Id sm            = 3;
  Id sh            = 10;
  VectorInt speeds = {sm, sm, sm, sm, sl, sl,
                      sh, sh, sh, sh, sl, sl,
                      sl, sl, sl, sl, sl, sl};
  (void)fluid_propagation(grid, "Facies", "Fluid", "", "", nfacies, nfluids, 1, speeds);

  grid->display(&dbfmt);
  (void)grid->dumpToNF("Grid.NF");

  delete grid;
  delete model;

  return (0);
}
