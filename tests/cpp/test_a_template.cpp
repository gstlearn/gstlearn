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

  double radius= 6371;
  ASpaceObject::defineDefaultSpace(ESpaceType::SPACE_SN, radius);
  String filename;

  MeshSpherical mesh = MeshSpherical();
  mesh.resetFromDb(nullptr,nullptr,"-r3",false);

  Model* model = Model::createFromParam(ECov::BESSEL_K,1500,1);
  ShiftOpCs S = ShiftOpCs(&mesh,model);
  VectorDouble whitenoise = VH::simulateGaussian(mesh.getNApices());
  VectorDouble result = VectorDouble(mesh.getNApices());
  PrecisionOpCs Q = PrecisionOpCs(&S,model->getCova(0),EPowerPT::MINUSHALF);
  Q.eval(whitenoise,result);

  Db* db = Db::create();
  VectorDouble X = mesh.getCoordinates(0);
  VectorDouble Y = mesh.getCoordinates(1);
  db->addColumns(X, "long", ELoc::X, 0);
  db->addColumns(Y, "lat", ELoc::X, 1);

  return (0);
}
