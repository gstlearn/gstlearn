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
#include "../../include/Geometry/GeometryHelper.hpp"
#include "geoslib_f.h"
#include "csparse_f.h"

#include "Db/Db.hpp"
#include "Basic/File.hpp"
#include "Basic/CSVformat.hpp"
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

  defineDefaultSpace(ESpaceType::SN, 2, EARTH_RADIUS);
  String filename;

  int nech = 40;
  VectorDouble extendmin = {0,0};
  VectorDouble extendmax = {150,100};
  Db* data = Db::createFromBox(nech, extendmin, extendmax);
  data->display();

  MeshEStandard mesh1 = MeshEStandard();
  (void) mesh1.resetFromDb(data,nullptr);
  mesh1.display();
  mesh1.printMeshes(2,10);

  return (0);
}
