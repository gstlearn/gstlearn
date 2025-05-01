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
#include "Enum/ECov.hpp"
#include "Enum/ESpaceType.hpp"

#include "Space/ASpaceObject.hpp"
#include "Basic/File.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Model/Model.hpp"
#include "Simulation/SimuFFTParam.hpp"
#include "Simulation/CalcSimuFFT.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This exercise is to demonstrate the Substitution simulation capability
 **
 *****************************************************************************/
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("SimFFT-");

  // Global parameters
  int ndim = 2;
  int seed = 3322;
  int nxcell = 100;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Generate the output grid
  VectorInt nx = {nxcell,nxcell};
  DbGrid* grid = DbGrid::create(nx);
  grid->display();

  // Create the Model
  Model* model = Model::createFromParam(ECov::SPHERICAL, 10.);
  model->display();

  message("\n<----- Creating Parameter File ----->\n");
  bool flag_aliasing = true;
  double percent = 0.1;
  SimuFFTParam param(flag_aliasing, percent);
  param.display();

  message("\n<----- Perform FFT Simulation ----->\n");
  (void) simfft(grid, model, param, 1, seed, true);

  (void) grid->dumpToNF("grid.ascii");

  delete grid;
  delete model;

  return (0);
}
