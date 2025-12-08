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
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Variogram/VCloud.hpp"
#include "Variogram/VarioParam.hpp"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Do not remove
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-"); // Here set the test name

  auto* db         = Db::createFillRandom(100, 2, 1);
  auto* varioparam = VarioParam::createOmniDirection(10, 0.1);
  auto* grid       = vcloudGrid(db, TEST, TEST, 10, 10);
  grid->display();

  VCloud vcloud(grid, varioparam);
  vcloud.setStorage(true);
  vcloud.compute(db);
  grid->display();

  vcloud.dumpStorage();
  vcloud.dumpStorage(1, 0.2, 0., 10, 10);
  vcloud.dumpStorage(2, 0.2, 0., 10, 10);
}
