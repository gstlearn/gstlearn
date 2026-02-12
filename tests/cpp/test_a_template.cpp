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
#include "Db/Db.hpp"
#include "Variogram/Vario.hpp"
#include "geoslib_define.h"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Do not remove
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-"); // Here set the test name

  Id ndat  = 100;
  auto* db = Db::createFillRandom(ndat, 2, 1);

  Id iech0           = 0;
  Id ilag0           = 2;
  double dlag        = 0.1;
  VectorDouble codir = {1., 0.};
  VectorInt ranks    = variogramPerPoint(db, iech0, ilag0, dlag, codir, 0.5, 20.);
  ranks.dump("Ranks");
  return 0;
}
