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
#include "Basic/NamingConvention.hpp"
#include "DataBase/ColID.hpp"
#include "Db/Db.hpp"
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
  auto ndim = 3;
  auto nvar = 2;
  auto nfex = 4;
  auto ncode = 2;
  double varmax = 3.;
  auto* db = Db::createFillRandom(nech, ndim, nvar, nfex, ncode, varmax);
  db->display();

  return 0;
}
