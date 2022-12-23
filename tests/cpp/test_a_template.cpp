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
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Calculators/CalcMigrate.hpp"
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
int main(int /*argc*/, char */*argv*/[])
{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  //  StdoutRedirect sr(sfn.str());

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("TestBugSpde-");

  bool flagSel = true;
  bool flagNoStat = true;

  DbGrid* grid = DbGrid::create({ 100, 100 });
  VectorDouble x1 = grid->getColumn("x1");
  VectorDouble x2 = grid->getColumn("x2");
  int size = (int) x2.size();

  if (flagSel)
  {
    VectorDouble selvec = VectorDouble(size, 0.);
    for (int i = 0; i < size; i++)
      selvec[i] = (x2[i] > 10.);
    grid->addColumns(selvec, "sel", ELoc::UNKNOWN);
    grid->setLocator("sel", ELoc::SEL);
  }

  MeshETurbo* mesh = MeshETurbo::createFromGrid(grid);

  NoStatArray nostat;
  if (flagNoStat)
  {
    VectorDouble temp = VectorDouble(size, 0.);
    for (int i = 0; i < size; i++)
      if (x1[i] > 50.) temp[i] = 40.;
    grid->addColumns(temp, "nostat", ELoc::NOSTAT);
    nostat = NoStatArray({"A"},grid);
  }

  Model* model = Model::createFromParam(ECov::BESSEL_K,0.,1.,1.,{20.,2.});
  if (flagNoStat) (void) model->addNoStat(&nostat);

  SPDE spde;
  spde.init(model,grid,nullptr,ESPDECalcMode::SIMUNONCOND,mesh);
  spde.compute();
  (void) spde.query(grid);

  (void) grid->dumpToNF("Grid.ascii");

  return (0);
}
