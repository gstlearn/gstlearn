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
#include "Covariances/CorAniso.hpp"
#include "geoslib_define.h"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Do not remove
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-"); // Here set the test name

  MatrixDense M1(3, 2);
  M1.fillRandom();

  message("Initial matrix M1\n");
  M1.display();

  // Testing the function-type of call
  M1.addCstInPlace(10.);
  message("After addition of 10 to M1\n");
  M1.display();

  auto M2 = M1.addCst(10.);

  message("After addition of 10 to M1 to create M2\n");
  M2->display();

  // Testing the operators
  M1 += 10.;
  message("After addition of 10 to M1\n");
  M1.display();

  auto M3 = M1 + 10.;
  message("After addition of 10 to M1 to create M3\n");
  M3->display();

  M1 -= 10.;
  message("After subtraction of 10 to M1\n");
  M1.display();

  auto M4 = M1 - 10.;
  message("After subtraction of 10 to M1 to create M4\n");
  M4->display();
  return 0;
}
