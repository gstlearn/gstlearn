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
using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Do not remove
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-"); // Here set the test name
  //

  tab_printg("my title", 24.);
  return 0;
}
