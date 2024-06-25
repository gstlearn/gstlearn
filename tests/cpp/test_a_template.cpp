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
#include "Mesh/MeshEStandard.hpp"
#include "Mesh/MeshSpherical.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "LinearOp/PrecisionOpCs.hpp"
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

  int option = 1;

  DbGrid* dbgrid = DbGrid::create({800,800,800},VectorDouble(),VectorDouble(),
                                  VectorDouble(),ELoadBy::COLUMN,
                                  VectorDouble(),VectorString(), VectorString(),
                                  0, false);
  VectorDouble x = VectorDouble(dbgrid->getSampleNumber());
  dbgrid->display();

  if (option == 0 || option == 1)
  {
    dbgrid->addColumnsByConstant(1, 12., "val0");
    dbgrid->display();

    dbgrid->deleteColumn("val0");
    dbgrid->display();
  }

  if (option == 0 || option == 2)
  {
    dbgrid->addColumns(x, "val");
    dbgrid->display();

    dbgrid->deleteColumn("val");
    dbgrid->display();
  }

  return (0);
}
