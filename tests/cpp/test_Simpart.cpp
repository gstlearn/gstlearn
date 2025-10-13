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
#include "Enum/ESpaceType.hpp"

#include "Basic/File.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Model/Model.hpp"
#include "Simulation/CalcSimuPartition.hpp"
#include "Simulation/SimuPartitionParam.hpp"
#include "Space/ASpaceObject.hpp"

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

  ASerializable::setPrefixName("test_Simpart-");

  // Global parameters
  Id ndim   = 2;
  Id seed   = 3322;
  Id nxcell = 100;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Generate the output grid
  VectorInt nx = {nxcell, nxcell};
  DbGrid* grid = DbGrid::create(nx);
  grid->display();

  // Create the Model
  Model* model = Model::createFromParam(ECov::SPHERICAL, 10.);
  model->display();

  // ====================== Create Parameter File ===================
  message("\n<----- Creating Parameter File ----->\n");
  Id nbtuba        = 50;
  double intensity = 0.1;
  SimuPartitionParam parparam(nbtuba, intensity);
  parparam.display();

  // ====================== Perform Boolean simulation ===================
  message("\n<----- Perform Partition Simulation ----->\n");
  (void)tessellation_poisson(grid, model, parparam, seed, false);

  (void)tessellation_voronoi(grid, model, parparam, seed, false);

  (void)grid->dumpToNF("grid.NF");

  delete grid;
  delete model;

  return (0);
}
