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
  double eps = 0.001;
  VectorString names1;
  VectorString names2;

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
  DbGrid* grid = DbGrid::create({30,20}, {1.,2.}, {10.,20.});

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
  int nclass = limits.getLimitNumber();

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

  // Testing PCA calculations
  grid->setLocator("Indicator.Simu.Class*", ELoc::Z);
  PCA pca = PCA();
  pca.pca_compute(grid,true);
  pca.display();

  // Back and Forth transforms
  pca.dbZ2F(grid, true);
  pca.dbF2Z(grid, true);
  names1 = generateMultipleNames("Indicator.Simu.Class" , nclass, ".");
  names2 = generateMultipleNames("F2Z.Z2F.Indicator.Simu.Class", nclass, ".");
  for (int iclass = 0; iclass < nclass; iclass++)
    (void) grid->areSame(names1[iclass], names2[iclass], eps);

  // Testing MAF calculations
  grid->setLocator("Indicator.Simu.Class*", ELoc::Z);
  DirParam dirparam = DirParam(2);
  PCA maf = PCA();
  maf.maf_compute(grid, 2., 1., dirparam, true);
  maf.display();

  // Back and Forth transforms
  grid->setLocator("Indicator.Simu.Class*", ELoc::Z);
  maf.dbZ2F(grid, true, NamingConvention("Z2MAF"));
  maf.dbF2Z(grid, true, NamingConvention("MAF2Z"));
  grid->display();
  names1 = generateMultipleNames("Indicator.Simu.Class" , nclass, ".");
  names2 = generateMultipleNames("MAF2Z.Z2MAF.Indicator.Simu.Class", nclass, ".");
  for (int iclass = 0; iclass < nclass; iclass++)
    (void) grid->areSame(names1[iclass], names2[iclass], eps);
  delete grid;
  return 0;
}

