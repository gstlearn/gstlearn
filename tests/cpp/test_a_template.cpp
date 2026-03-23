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
#include "Basic/ASerializable.hpp"
#include "LinearOp/PrecisionOpMatrix.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Model/Model.hpp"
#include "geoslib_define.h"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Do not remove
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-"); // Here set the test name
  Id nx = 350;

  auto mesh = MeshETurbo({nx, nx}, {}, {}, {}, false);
  auto range = nx / 2.;
  auto sill = 1.;
  auto param = 1.;
  auto* model = Model::createFromParam(ECov::MATERN, range, sill, param);
  message("Number of apices = %d\n", mesh.getNApices());

  auto QOpCs = PrecisionOpMatrix(&mesh, model->getCovAniso(0));
  return 0;
}
