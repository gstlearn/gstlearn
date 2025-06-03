/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/*                                                                            */
/* This file is meant to demonstrate the process of using PGS                 */
/*                                                                            */
/******************************************************************************/
#include "Enum/ECalcVario.hpp"
#include "Enum/ECov.hpp"

#include "Variogram/Vario.hpp"
#include "Variogram/VMap.hpp"
#include "Model/Model.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/File.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

/****************************************************************************/
/*!
** Main Program for testing the sparse matrix algebra
**
*****************************************************************************/
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setPrefixName("Vmap3D-");

  int ndim = 3;
  defineDefaultSpace(ESpaceType::RN, ndim);
  CovContext ctxt(1,ndim,1.); // use default space

  // Creating a grid
  VectorInt nx = { 50, 40, 20 };
  DbGrid* grid = DbGrid::create(nx);
  grid->display();

  // Creating the Model(s) of the Underlying GRF(s)
  Model* model = Model::createFromParam(ECov::SPHERICAL, 0., 1., 0., {1., 1., 0.1});
  model->display();

  // Perform a non-conditional simulation on the Db and on the Grid
  (void) simtub(nullptr,grid,model);
  grid->display();

  // =================================
  // Calculating Variogram Map on Grid
  // =================================

  DbGrid* vmap = db_vmap(grid, ECalcVario::VARIOGRAM,{10,10,3});
  DbStringFormat dbfmt(FLAG_STATS,{"VMAP*"});
  vmap->display(&dbfmt);

  (void) vmap->dumpToNF("vmap.ascii");

  delete grid;
  delete vmap;
  delete model;
  
  return 0;
}
