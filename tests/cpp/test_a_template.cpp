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
#include "Db/Db.hpp"
#include "Basic/File.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Variogram/DirParam.hpp"
#include "Variogram/VarioParam.hpp"
#include "Variogram/Vario.hpp"
#include "Model/Model.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Geometry/Geometry.hpp"
#include "Calculators/CalcMigrate.hpp"

/**
 * This file is meant to perform any test that needs to be coded for a quick trial
 */
int main(int /*argc*/, char */*argv*/[])

{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  ASpaceObject::defineDefaultSpace(ESpaceType::SPACE_RN, 2);

  String filename = ASerializable::getTestData("Scotland","temperatures.ascii");
  Db* temperatures = Db::createFromNF(filename);
  temperatures->display();

  int nx = 50;
  double x0 = 0;
  double dx = 1;
  DbGrid* dbTarget = DbGrid::create(nx={nx+1,nx+1},x0={x0,x0},dx={dx,dx});
  dbTarget->display();
//  dbTarget = db.locate(dbTarget,2:3,"x")

//  nObs = 5
//  set.seed(99)
//  coordObs = cbind(sample(seq(nx),nObs,replace=F),sample(seq(nx),nObs,replace=F))
//  valObs = rnorm(nObs)
//  dbObs = db.create(cbind(coordObs,valObs))
//  dbObs = db.locate(dbObs,2:3,"x")
//  dbObs = db.locate(dbObs,4,"z")
//  plot(dbObs,pos.legend=1,title="Observations")
//
//  filename = ASerializable::getTestData("Scotland","Scotland_Elevations.csv");
//  Db* mnt = Db::createFromCSV(filename);
//  mnt->setLocators({"Longitude","Latitude"},ELoc::X);
//  mnt->display();
//
//  DbGrid* grid = DbGrid::createCoveringDb(mnt, {81,137});
//  grid->display();
//
//  (void) migrateVariables(mnt,grid,{"Elevation","inshore"});
//  grid->display();

  return (0);
}
