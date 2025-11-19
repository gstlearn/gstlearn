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
#include "Stats/Classical.hpp"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Do not remove
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-"); // Here set the test name
  

  auto* db                       = Db::createFillRandom(100, 2, 3);
  VectorString names             = {"x-*", "z-*"};
  std::vector<EStatOption> opers = EStatOption::fromKeys({"NUM", "MINI", "MAXI", "MEAN", "STDV"});
  dbStatisticsPrint(db, names, opers);
}
