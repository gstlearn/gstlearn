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
#include "Basic/Timer.hpp"
#include "Db/Db.hpp"
#include "Variogram/Vario.hpp"
#include "Variogram/VarioParam.hpp"
#include "geoslib_define.h"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Unless you test a specific feature, bring it to its minimal expansion
  // e.g. the few lines below
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-");

  auto* db = Db::createFillRandom(5000, 2, 1);
  auto* varioparam = VarioParam::createOmniDirection(10, 0.05);
  Timer timer;

  timer.reset();
  auto* vario = Vario::computeFromDb(*varioparam, db);
  timer.displayIntervalMilliseconds("Variogram calculation", 5);
  vario->display();
  return (0);
}
