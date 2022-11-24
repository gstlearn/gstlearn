/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_d.h"
#include "geoslib_f.h"

#include "Enum/ELoadBy.hpp"
#include "Enum/ESpaceType.hpp"

#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "Basic/VectorHelper.hpp"
#include "Boolean/ShapeEllipsoid.hpp"
#include "Boolean/ShapeParallelepiped.hpp"
#include "Boolean/ModelBoolean.hpp"
#include "Simulation/SimuBooleanParam.hpp"

static Db* createLocalDb(int nech, int ndim, int nvar,
                         bool flag_sel = false, double proba = 0.5)
{
  // Coordinates
  VectorDouble tab = VH::simulateGaussian(ndim * nech, 0., 50.);
  // Variable
  for (int ivar=0; ivar<nvar; ivar++)
  {
    VectorDouble tabvar;
    if (flag_sel)
      tabvar = VH::simulateBernoulli(nech, proba);
    else
      tabvar = VH::simulateGaussian(nech);
    tab.insert(tab.end(), tabvar.begin(), tabvar.end());
  }

  Db* data = Db::createFromSamples(nech,ELoadBy::COLUMN,tab);
  data->setNameByUID(1,"x1");
  data->setNameByUID(2,"x2");

  data->setLocatorByUID(1,ELoc::X,0);
  data->setLocatorByUID(2,ELoc::X,1);

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    data->setNameByUID(3+ivar,"Var");
    data->setLocatorByUID(3+ivar,ELoc::Z,ivar);
  }
  return data;
}

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This exercise is to demonstrate the Boolean simulation capability
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("SimBool-");

  // Global parameters
  law_set_random_seed(32131);
  int ndim = 2;
  int nvar = 1;
  int nxcell = 100;
  int nech = 100;
  VectorDouble coormin(ndim);
  VectorDouble coormax(ndim);
  defineDefaultSpace(ESpaceType::RN, ndim);
  DbStringFormat dbfmt(FLAG_STATS);

  // Generate the output grid
  VectorInt nx = {nxcell,nxcell};
  DbGrid* grid = DbGrid::create(nx);
  grid->display();

  // Generate the data base
  Db* data = createLocalDb(nech, ndim, nvar, true, 0.1);
  data->display(&dbfmt);

  // ====================== Create Shape Dictionary ===================
  message("\n<----- Creating Shape Dictionary ----->\n");
  ModelBoolean* tokens = new ModelBoolean(0.01, true);
  ShapeEllipsoid token_ellipsoid(0.4, 10., 20., 2.);
  token_ellipsoid.setFactorX2Y(1.5);
  tokens->addToken(token_ellipsoid);
  ShapeParallelepiped token_parallelepiped(0.6, 5, 7., 1.);
  tokens->addToken(token_parallelepiped);
  tokens->display();

  // ====================== Perform Boolean simulation ===================
  message("\n<----- Perform Boolean Simulation ----->\n");
  (void) simbool(nullptr, grid, tokens);

  grid->display(&dbfmt);

  (void) grid->dumpToNF("grid.ascii");

  delete grid;
  delete data;

  return (0);
}
