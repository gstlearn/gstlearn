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

#include "geoslib_define.h"

#include "Basic/ASerializable.hpp"
#include "Db/DbGrid.hpp"
#include "Variogram/Vario.hpp"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Do not remove
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-"); // Here set the test name

  auto* dbgrid = DbGrid::createFillRandom({10, 10}, 2);
  dbgrid->display();

  auto* vario1 = variogramCalculate(dbgrid, 5, 0.1, 1, {}, true);
  if (vario1 != nullptr) vario1->display();

  auto* vario2 = variogridCalculate(dbgrid, 5, true);
  if (vario2 != nullptr) vario2->display();
}
