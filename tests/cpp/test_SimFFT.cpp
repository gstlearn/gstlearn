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

#include "Basic/File.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Model/Model.hpp"
#include "Simulation/CalcSimuFFT.hpp"
#include "Simulation/SimuFFTParam.hpp"
#include "Simulation/Simulations.hpp"
#include "Space/ASpaceObject.hpp"
#include "Variogram/Vario.hpp"
#include "Variogram/VarioParam.hpp"

using namespace gstlrn;

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This exercise is to demonstrate the Substitution simulation capability
 **
 *****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setPrefixName("test_SimuFFT-");
  auto param = SimuFFTParam(false, 1.);
  bool verbose = true;
  auto nbsimu = 1;
  auto seed = 243131;

  // Simulate an isotropic Spherical Model on a 2D grid
  // ==================================================
  mestitle(0, "2D case");
  defineDefaultSpace(ESpaceType::RN, 2);

  // Generate the output grid
  DbGrid* grid2D = DbGrid::create({100, 100});
  grid2D->display();

  // Create the Model
  Model* model2D = Model::createFromParam(ECov::SPHERICAL, 20.);
  model2D->display();

  // Perform the simulation
  (void)simuFFT(nullptr, grid2D, model2D, param, nbsimu, seed, verbose);

  // Compute experimental variograms in main grid directions
  auto* varioparam2D = VarioParam::createMultipleFromGrid(grid2D, 30);
  auto vario2D = Vario(*varioparam2D);
  (void)vario2D.compute(grid2D);

  // Dump the results into Neutral Files
  (void)grid2D->dumpToNF("grid2D.NF");
  (void)model2D->dumpToNF("model2D.NF");
  (void)vario2D.dumpToNF("vario2D.NF");

  delete grid2D;
  delete model2D;
  delete varioparam2D;

  // Simulate an anisotropic Exponential Model on a 3D grid
  // ======================================================
  mestitle(0, "3D case");
  defineDefaultSpace(ESpaceType::RN, 3);

  // Generate the output grid
  DbGrid* grid3D = DbGrid::create({64, 64, 64});
  grid3D->display();

  // Create the Model
  Model* model3D =
    Model::createFromParam(ECov::EXPONENTIAL, 20., 1., 1., {20., 10., 5.});
  model3D->display();

  // Perform the simulation
  (void)simuFFT(nullptr, grid3D, model3D, param, nbsimu, seed, verbose);

  // Compute experimental variograms in main grid directions
  auto* varioparam3D = VarioParam::createMultipleFromGrid(grid3D, 30);
  auto vario3D = Vario(*varioparam3D);
  (void)vario3D.compute(grid3D);

  // Dump the results into Neutral Files
  (void)grid3D->dumpToNF("grid3D.NF");
  (void)model3D->dumpToNF("model3D.NF");
  (void)vario3D.dumpToNF("vario3D.NF");

  delete grid3D;
  delete model3D;
  delete varioparam3D;

  return (0);
}
