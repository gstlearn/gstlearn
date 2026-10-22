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
#include "Enum/ELoc.hpp"
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

  auto nech = 10;
  auto nvers = 5;
  auto* db = Db::createFillRandom(nech);
  db->addColumnsByConstant(1, nvers, 0., "myf", ELoc::SIMU);

  auto colid = ColID(4, 13);
  message("%d\n", db->hasColumn(colid));
  message("%d\n", db->hasValidVersion(colid));

  exit(0);

  auto name = db->getName(ColID("z", 1), true);
  message("Name = %s\n", name.c_str());
  auto value = db->getValue(ColID("z", 1), 0);
  message("Value = %lf\n", value);
  return 0;
}
