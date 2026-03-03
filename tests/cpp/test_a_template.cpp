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

  ////////////////////////////////
  // Testing the free functions //
  ////////////////////////////////
  mestitle(0, " Testing free functions");

  AMatrix::addCstT(M1, M1, 10.);
  message("After: M1 = M1 + 10 (free method: addCstT)\n");
  M1.display();

  AMatrix::addCstT(MS1, MS1, 10.);
  message("After: MS1 = MS1 + 10. (free method: addCstT)\n");
  MS1.display();

  AMatrix::prodCstT(M1, M1, 2.);
  message("After: M1 = M1 * 2 (free method: prodCstT)\n");
  M1.display();

  AMatrix::prodCstT(MS1, MS1, 2.);
  message("After: MS1 = MS1 * 2. (free method: prodCstT)\n");
  MS1.display();

  //////////////////////////////////
  // Testing the member functions //
  //////////////////////////////////
  mestitle(0, " Testing member functions");

  M1.addCstInPlace(10.);
  message("After: M1 += 10. (class method: addCstInPlace)\n");
  M1.display();

  auto M1add = M1.addCst(10.);
  message("After: M1add = M1 + 10 (class method: addCst)\n");
  M1add.display();

  M1.prodCstInPlace(2.);
  message("After: M1 *= 2. (class method: prodCstInPlace)\n");
  M1.display();

  auto M1prod = M1.prodCst(2.);
  message("After: M1prod = M1 * 2 (class method: prodCst)\n");
  M1prod.display();

  ///////////////////////////
  // Testing the operators //
  ///////////////////////////
  mestitle(0, " Testing the operators");

  M1 += 10.;
  message("After: M1 += 10. (operator +=)\n");
  M1.display();

  auto M1addOp = M1 + 10.;
  message("After: M1addOp = M1 + 10 (operator +)\n");
  M1addOp.display();

  MS1 += 10.;
  message("After: MS1 += 10. (operator +=)\n");
  MS1.display();

  auto MS1addOp = MS1 + 10.;
  message("After: MS1addOp = MS1 + 10 (operator +)\n");
  MS1addOp.display();

  return 0;
}
