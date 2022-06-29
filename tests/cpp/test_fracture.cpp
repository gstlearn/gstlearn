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

// This test is meant to demonstrate the fracture Simulation

#include "geoslib_d.h"
#include "geoslib_f.h"
#include "Space/Space.hpp"
#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/OptCustom.hpp"
#include "Fractures/FracEnviron.hpp"
#include "Fractures/FracFamily.hpp"
#include "Fractures/FracFault.hpp"
#include "Fractures/FracList.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
//  StdoutRedirect sr(sfn.str());

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Fractures-");

  // Global parameters
  int ndim = 2;
  law_set_random_seed(32131);

  ASpaceObject::defineDefaultSpace(SPACE_RN, ndim);
  DbStringFormat dbfmt(FLAG_STATS);

  // Generate the output grid
  VectorInt nx = {301,101};
  DbGrid* grid = DbGrid::create(nx);
  grid->display();

  // Creating the Fracture Environment
  double xmax = grid->getExtend(0);
  double ymax = grid->getExtend(1);
  double deltax = 0.;
  double deltay = 0.;
  double mean   = 20.;
  double stdev  = 10.;
  FracEnviron env = FracEnviron(xmax, ymax, deltax, deltay, mean, stdev);

  // Creating the Fault Families
  double orient  = 0.;
  double dorient = 20.;
  double theta0  = 1.;
  double alpha   = 1.;
  double ratcst  = 1;
  double prop1   = 0.5;
  double prop2   = 0.2;
  double aterm   = 1.2;
  double bterm   = 2.4;
  double range   = 12.;
  FracFamily family = FracFamily(orient, dorient, theta0, alpha, ratcst,
                                 prop1, prop2, aterm, bterm, range);
  family.display();

  // Creating the Major Fault
  double coord   = 30.;
  double forient = 0.;
  FracFault fault = FracFault(coord, forient);
  double thetal  = 1.;
  double thetar  = 2.;
  double rangel  = 10.;
  double ranger  = 20.;
  fault.addFaultPerFamily(thetal, thetar, rangel, ranger);
  fault.display();

  // Gluing all elements within the Environment
  env.addFault(fault);
  env.addFamily(family);
  env.display();

  // Simulating fractures
  FracList flist = FracList();
  int seed = 432431;
  flist.simulate(env, true, true, seed, false, VectorDouble());
  flist.display();

  // Plunge the set of fractures on the Grid
  VectorDouble permtab = { 10., 20., 30. };
  double perm_mat   = 0.;
  double perm_bench = 5.;
  (void) flist.fractureToBlock(grid, xmax, permtab, perm_mat, perm_bench);

  MatrixRectangular layinfo = flist.layinfoExport();
  layinfo.display();
  MatrixRectangular fracinfo = flist.fractureExport();
  fracinfo.display();

  grid->display(&dbfmt);

  // Save as Neutral File

  (void) grid->dumpToNF("Grid.ascii");

  // ====================== Free pointers ==================================
  if (grid != nullptr) delete grid;

  return (0);
}
