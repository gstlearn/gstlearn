#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Basic/OptCst.hpp"
#include "Basic/ECst.hpp"
#include "Basic/FunctionalSpirale.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "API/SPDE.hpp"
#include "Model/Model.hpp"
#include "Model/NoStatArray.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Stats/PCA.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

#include "geoslib_f.h"
#include "geoslib_f_private.h"
#include "geoslib_old_f.h"

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program is meant to check the manipulation of the Db
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  VectorString names1;
  VectorString names2;

  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  DbStringFormat dbfmt(FLAG_STATS);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("TestLimits-");
  int seed = 10355;
  law_set_random_seed(seed);

  // Creating the Grid Rotated Db
  DbGrid* grid = DbGrid::create({30,20}, {1.,2.}, {10.,20.});

  // Creating the Model
  Model model(1, 2);
  model.addCovFromParam(ECov::CUBIC, 0., 2., 1., {10.,45.}, {}, {30.,0.});
  model.display();

  // Simulating a variable on the grid
  (void) simtub(nullptr, grid, &model, nullptr);
  dbfmt = DbStringFormat(FLAG_STATS, {"Simu"});
  grid->display(&dbfmt);

  // Creating a set of Limits
  Limits limits({-1.,-0.5,0.,0.5,1.});
  limits.display();

  // Other option
  grid->setLocator("Simu", ELoc::Z);
  limits.toIndicator(grid,"Simu",0);
  dbfmt = DbStringFormat(FLAG_ARRAY, {"Indicator.Simu.Mean"});
  grid->display(&dbfmt);

  // Convert into Indicators
  grid->setLocator("Simu", ELoc::Z);
  limits.toIndicator(grid,"Simu",1);
  dbfmt = DbStringFormat(FLAG_ARRAY, {"Indicator.Simu.Class*"});
  grid->display(&dbfmt);

  delete grid;
  return 0;
}

