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
#include "Db/DbGrid.hpp"
#include "Enum/ECst.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "geoslib_define.h"

#include "Basic/ASerializable.hpp"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Do not remove
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-"); // Here set the test name
  //

  gstlrn::DbGrid* db_out = gstlrn::DbGrid::create({5, 5});
  db_out->addColumns(
    {
      1.,
      1.,
      1.,
      1.,
      1.,
      1.,
      0.,
      1.,
      1.,
      1.,
      1.,
      1.,
      1.,
      1.,
      1.,
      1.,
      1.,
      1.,
      1.,
      1.,
      1.,
      1.,
      1.,
      1.,
      1.,
    },
    "SEL", gstlrn::ELoc::SEL);
  gstlrn::MeshETurbo* mesh = gstlrn::MeshETurbo::createFromGrid(db_out);

  gstlrn::ProjMatrix* proj2D = gstlrn::ProjMatrix::create(db_out, mesh, -1, true);
  OptCst::define(ECst::NTCOL, -1);
  OptCst::define(ECst::NTROW, -1);
  proj2D->display();
  return 0;
}
