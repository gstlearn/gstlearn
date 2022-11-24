/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/*                                                                            */
/* This file is meant to demonstrate the process of using PGS                 */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f.h"

#include "Enum/ECalcVario.hpp"
#include "Enum/ECov.hpp"

#include "Variogram/Vario.hpp"
#include "Model/Model.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/File.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

#include <stdlib.h>

/****************************************************************************/
/*!
** Main Program for testing the sparse matrix algebra
**
*****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  ASerializable::setContainerName(true);
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

  DbGrid* vmap = db_vmap_compute(grid, ECalcVario::VARIOGRAM,{10,10,3});
  DbStringFormat dbfmt(FLAG_STATS,{"VMAP*"});
  vmap->display(&dbfmt);

  (void) vmap->dumpToNF("vmap.ascii");

  delete grid;
  delete vmap;
  return 0;
}
