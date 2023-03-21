/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/

// This test is meant to demonstrate the Eden Simulation
// which also demonstrates the SKIN methodology

#include "geoslib_d.h"

#include "Enum/ECov.hpp"
#include "Enum/ESpaceType.hpp"

#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/OptCustom.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Simulation/CalcSimuEden.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Eden-");

  // Global parameters
  int ndim = 2;
  law_set_random_seed(32131);

  defineDefaultSpace(ESpaceType::RN, ndim);
  DbStringFormat dbfmt(FLAG_STATS);

  // Generate the output grid
  VectorInt nx = {301,101};
  DbGrid* grid = DbGrid::create(nx);

  // Simulate continuous variable
  double range = 30.;
  Model* model = Model::createFromParam(ECov::CUBIC, range);
  (void) simtub(nullptr, grid, model);

  // Operate a transformation to convert into 3 (nested) facies
  int nfacies = 3;
  double thresh = 1.;
  VectorDouble vmini = { TEST, -thresh, thresh };
  VectorDouble vmaxi = { -thresh, thresh, TEST };

  Limits* limits = Limits::create(vmini, vmaxi);
  limits->toCategory(grid,"Simu");
  grid->setName("Category*", "Facies");

  // Add a Fluid information
  int nfluids = 1;
  grid->addColumnsByConstant(nfluids, TEST, "Fluid");
  (void) grid->assignGridColumn("Fluid", 0, 100, 1.);

  // Fluid propagation
  // Speed: (dir + 6 * (facies * nfluids + fluid))
  // Direction = 0: +X; 1: -X; 2: +Y; 3: -Y; 4: +Z(up); 5: -Z(down)

  int sl = 1;
  int sm = 3;
  int sh = 10;
  VectorInt speeds = { sm, sm, sm, sm, sl, sl,
                       sh, sh, sh, sh, sl, sl,
                       sl, sl, sl, sl, sl, sl};
  (void) fluid_propagation(grid,"Facies","Fluid","","",nfacies,nfluids,1,speeds);

  grid->display(&dbfmt);
  (void) grid->dumpToNF("Grid.ascii");

  // ====================== Free pointers ==================================
  if (grid != nullptr) delete grid;
  if (model != nullptr) delete model;

  return (0);
}

