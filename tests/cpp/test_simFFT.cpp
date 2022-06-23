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
#include "Space/Space.hpp"
#include "Space/ASpaceObject.hpp"
#include "Basic/File.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Covariances/ECov.hpp"
#include "Simulation/SimuFFTParam.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This exercise is to demonstrate the Substitution simulation capability
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("SimFFT-");

  // Global parameters
  int ndim = 2;
  int seed = 3322;
  int nxcell = 100;
  ASpaceObject::defineDefaultSpace(SPACE_RN, ndim);

  // Generate the output grid
  VectorInt nx = {nxcell,nxcell};
  DbGrid* grid = DbGrid::create(nx);
  grid->display();

  // Create the Model
  Model* model = Model::createFromParam(ECov::SPHERICAL, 10.);
  model->display();

  // ====================== Create Parameter File ===================
  message("\n<----- Creating Parameter File ----->\n");
  bool flag_aliasing = true;
  double percent = 0.1;
  SimuFFTParam param(flag_aliasing, percent);
  param.display();

  // ====================== Perform Boolean simulation ===================
  message("\n<----- Perform FFT Simulation ----->\n");
  (void) simfft(grid, model, param, 1, seed, true);

  (void) grid->dumpToNF("grid.ascii");

  delete grid;
  delete model;

  return (0);
}
