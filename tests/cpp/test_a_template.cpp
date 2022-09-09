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
#include "geoslib_f.h"
#include "csparse_f.h"

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
#include "Mesh/MeshETurbo.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "Stats/Classical.hpp"

/**
 * This file is meant to perform any test that needs to be coded for a quick trial
 * It will be compiled but not run (not diffed)
 */
int main(int /*argc*/, char */*argv*/[])

{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  //  StdoutRedirect sr(sfn.str());

  ASpaceObject::defineDefaultSpace(ESpaceType::SPACE_RN, 2);

  Model* model = Model::createFromParam(ECov::CUBIC,20,2);

  int nech = 30;
  Db* data = Db::createFromBox(nech, {0,0}, {100,100});
  simtub(nullptr,data, model);
  data->setName(data->getLastName(), "data");
  data->display();

  DbGrid* grid = DbGrid::create({100,100});
  grid->display();

  NeighUnique* neighU = NeighUnique::create(2);
  neighU->display();

  (void) simtub(data, grid, model, neighU, 10);
  grid->display();

  (void) grid->statistics({"Simu*"}, { EStatOption::QUANT }, true, false, false,
                          0.2, TEST, TEST);
  DbStringFormat dbfmt;
  dbfmt.setFlags(true, true, false, true, false, false, {"*QUANT"});
  grid->display(&dbfmt);

  return (0);
}
