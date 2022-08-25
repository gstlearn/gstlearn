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
#include "geoslib_old_f.h"
#include "Db/Db.hpp"
#include "Basic/File.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

/**
 * This file is meant to perform any test that needs to be coded for a quick trial
 */
int main(int /*argc*/, char */*argv*/[])

{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
//  StdoutRedirect sr(sfn.str());

  // Generate the output grid
  DbGrid* grid = DbGrid::create({100,100},{0,1},{1.1,1.2});
  grid->display();

  Model* model = Model::createFromParam(ECov::CUBIC, 20., 3.);
  model->display();

  (void) simtub(nullptr, grid, model, nullptr, 2);
  grid->display();

  DbGrid* grid1 = DbGrid::createFromGridExtend(*grid,{"Simu.1"},{"Simu.2"},{10},true);
  grid1->display();

  DbGrid* grid2 = DbGrid::createFromGridShrink(*grid1, {1});
  grid2->display();

  // Free the pointers
  delete grid;
  delete grid1;
  delete model;

  return (0);
}
