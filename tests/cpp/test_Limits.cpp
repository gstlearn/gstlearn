#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Basic/OptCst.hpp"
#include "Basic/ECst.hpp"
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
#include "Basic/FunctionalSpirale.hpp"
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
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
//  StdoutRedirect sr(sfn.str());

  DbStringFormat dbfmt(FLAG_STATS);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("TestLimits-");
  int seed = 10355;
  law_set_random_seed(seed);

  // Creating the Grid Rotated Db
  DbGrid* grid = DbGrid::create({6,4}, {1.,2.}, {10.,20.}, {10.,0.});

  // Creating the Model
  Model model(1, 2);
  model.addCova(ECov::CUBIC, 0., 2., 1., {10.,45.}, {}, {30.,0.});
  model.display();

  // Simulating a variable on the grid
  (void) simtub(nullptr, grid, &model, nullptr);
  dbfmt = DbStringFormat(FLAG_STATS, {"Simu"});
  grid->display(&dbfmt);

  // Creating a set of Limits
  Limits limits({-1.,-0.5,0.,0.5,1.});
  limits.display();

  // Convert into Indicators
  limits.toIndicator(grid,"Simu");
  dbfmt = DbStringFormat(FLAG_ARRAY, {"Simu","Indicator*"});
  dbfmt = DbStringFormat(FLAG_STATS, {"Indicator*"});
  grid->display(&dbfmt);

  // Other option
  limits.toIndicator(grid,"Simu",0);
  dbfmt = DbStringFormat(FLAG_ARRAY, {"Simu","Indicator*"});
  grid->display(&dbfmt);

  delete grid;
  return 0;
}

