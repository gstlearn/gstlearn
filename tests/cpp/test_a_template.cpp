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
#include "Db/DbGrid.hpp"
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

  auto* grid = DbGrid::create({3, 4}, {0.0, 0.0}, {1.0, 1.0});
  grid->display();
  grid->addColumnsRandom(1, 4, "MyVar");
  grid->display();

  grid->dumpToNF("avoir.dat");
  return (0);
}
