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
#include "Matrix/AMatrix.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Matrix/MatrixSquare.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "geoslib_define.h"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Do not remove
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-"); // Here set the test name

  auto M1    = MatrixSquare(7, 7);
  auto M2    = MatrixSymmetric(7, 7);
  AMatrix* M = MatrixFactory::prodMatMat(&M1, &M2);
  M->printConcreteClassName();

  return 0;
}
