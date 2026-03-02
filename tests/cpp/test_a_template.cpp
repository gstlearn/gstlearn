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
#include "Matrix/MatrixSparse.hpp"
#include "geoslib_define.h"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Do not remove
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-"); // Here set the test name

  MatrixDense M1(4, 3);
  M1.fillRandom();
  MatrixSparse MS1(4, 3);
  MS1.fillRandom(31366, 0.3);

  message("Initial Dense Matrix M1\n");
  M1.display();
  message("Initial Sparse Matrix MS1\n");
  MS1.display();

  AMatrix::addCstT(M1, M1, 10.);
  message("After: M1 = M1 + 10 (free method: addCstT)\n");
  M1.display();

  AMatrix::addCstT(MS1, MS1, 10.);
  message("After: MS1 = MS1 + 10. (free method: addCstT)\n");
  MS1.display();

  M1.addCstInPlace(10.);
  message("After: M1 += 10. (class method: addCstInPlace)\n");
  M1.display();

  auto M2 = M1.addCst(10.);
  message("After: M2 = M1 + 10 (class method: addCst)\n");
  M2.display();

  // Testing the operators
  M1 += 10.;
  message("After: M1 += 10. (operator +=)\n");
  M1.display();

  auto M3 = M1 + 10.;
  message("After: M3 = M1 + 10 (operator +)\n");
  M3.display();

  // Testing the operators
  MS1 += 10.;
  message("After: MS1 += 10. (operator +=)\n");
  MS1.display();

  auto MS3 = MS1 + 10.;
  message("After: MS3 = MS1 + 10 (operator +)\n");
  MS3.display();

  return 0;
}
