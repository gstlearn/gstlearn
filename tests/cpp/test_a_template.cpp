/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f.h"
#include "Matrix/LinkMatrixSparse.hpp"

#include "Enum/ESPDECalcMode.hpp"

#include "Db/Db.hpp"
#include "Basic/File.hpp"
#include "Basic/CSVformat.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Variogram/DirParam.hpp"
#include "Variogram/VarioParam.hpp"
#include "Variogram/Vario.hpp"
#include "Model/Model.hpp"
#include "Model/NoStatArray.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Geometry/GeometryHelper.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Mesh/MeshEStandard.hpp"
#include "Mesh/MeshSpherical.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "LinearOp/PrecisionOpCs.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "Stats/Classical.hpp"
#include "Polygon/Polygons.hpp"
#include "API/SPDE.hpp"

/**
 * This file is meant to perform any test that needs to be coded for a quick trial
 * It will be compiled but not run nor diff'ed.
 */
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("AAA_");

  DbGrid* grid = DbGrid::createFromNF("db.ascii");
  grid->display();

  Model* model = Model::createFromNF("Model.ascii");
  NoStatArray NoStat({"A","R"},grid);
  model->addNoStat(&NoStat);
  model->display();

  model->display();

  SPDE spde(model,grid,nullptr,ESPDECalcMode::SIMUNONCOND);
  spde.compute(grid);
  (void) grid->dumpToNF("Result.ascii");

  return (0);
}
