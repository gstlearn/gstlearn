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

  String filename = "/home/drenard/TEST_gstlearn/test_SimGaussian-Vario_TB";
  auto filetype = ASerializable::getFileIdentity(filename);
  auto* vario = Vario::createFromNF(filename, false);
  if (!vario)
  {
    std::cerr << "Error creating Vario from file: " << filename << std::endl;
    return 1;
  }
  return 0;
}
