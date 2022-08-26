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

#include "Space/ESpaceType.hpp"
#include "Space/ASpaceObject.hpp"
#include "Basic/File.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Covariances/ECov.hpp"
#include "Simulation/SimuRefineParam.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This exercise is to demonstrate the Refinement simulation capability
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("SimRefine-");

  // Global parameters
  int ndim = 2;
  int seed = 3322;
  int nxcell = 5;
  ASpaceObject::defineDefaultSpace(ESpaceType::SPACE_RN, ndim);

  // Generate the output grid
  VectorInt nx = {nxcell,nxcell};
  DbGrid* grid = DbGrid::create(nx);
  grid->display();

  // Create the Model
  Model* model = Model::createFromParam(ECov::SPHERICAL, 10.);
  model->display();

  // ====================== Create Parameter File ===================
  message("\n<----- Creating Reference Simulation ----->\n");
  simtub(nullptr, grid, model, nullptr);
  (void) grid->dumpToNF("grid_small.ascii");

  // ====================== Create Parameter File ===================
  message("\n<----- Creating Parameter File ----->\n");
  int nmult = 5;
  SimuRefineParam param(nmult);
  param.display();

  // ====================== Perform Boolean simulation ===================
  message("\n<----- Perform Refinement Simulation ----->\n");
  DbGrid* grid2 = simfine(grid, model, param, seed);
  (void) grid2->dumpToNF("grid_large.ascii");

  delete grid;
  delete grid2;
  delete model;

  return (0);
}
