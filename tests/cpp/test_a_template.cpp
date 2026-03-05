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
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "geoslib_define.h"

using namespace gstlrn;

void reset(MatrixDense& M1,
           MatrixDense& M2,
           MatrixSparse& MS1,
           MatrixSparse& MS2,
           bool verbose = false)
{
  MS1.fillRandom(0.3, 31366);
  MS2.fillRandom(0.3, 424672);
  MatrixDense* Mt1 = MatrixSparse::createFromSparse(MS1);
  MatrixDense* Mt2 = MatrixSparse::createFromSparse(MS2);
  M1               = *Mt1;
  M2               = *Mt2;
  delete Mt1;
  delete Mt2;

  if (verbose)
  {
    message("Initial Sparse Matrix MS1\n");
    MS1.display();
    message("Initial Sparse Matrix MS2\n");
    MS2.display();
    message("Initial Dense Matrix M1\n");
    M1.display();
    message("Initial Dense Matrix M2\n");
    M2.display();
  }
  else
  {
    message("Matrices M1, M2, MS1 and MS2 are reset to initial values\n");
  }
}

void resetSquare(MatrixSquare& MSq1,
                 MatrixSquare& MSq2,
                 bool verbose = false)
{
  MSq1.fillRandom(0., 31366);
  MSq2.fillRandom(0., 424672);

  if (verbose)
  {
    message("Initial Square Matrix MSq1\n");
    MSq1.display();
    message("Initial Square Matrix MSq2\n");
    MSq2.display();
  }
  else
  {
    message("Matrices MSq1 and MSq2 are reset to initial values\n");
  }
}

int main(int argc, char* argv[])
{
  // Do not remove
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-"); // Here set the test name

  // Creating two matrix sparse and derive the equivalent MatrixDense
  // so that comparative calculations can be run on the two types of matrix.
  MatrixSparse MS1(4, 3);
  MatrixSparse MS2(4, 3);
  MatrixDense M1;
  MatrixDense M2;
  reset(M1, M2, MS1, MS2, true);

  ////////////////////////////////
  // Testing the free functions //
  ////////////////////////////////
  mestitle(0, "Testing free functions");

  // Adding a Constant
  reset(M1, M2, MS1, MS2, false);

  AMatrix::add(M1, M1, 10.);
  message("AMatrix::add(M1, M1, 10.) -> M1\n");
  M1.display();

  AMatrix::add(MS1, MS1, 10.);
  message("AMatrix::add(MS1, MS1, 10.) -> MS1\n");
  MS1.display();

  // Adding two matrices
  reset(M1, M2, MS1, MS2, false);

  AMatrix::add(M2, M1, M2);
  message("AMatrix::add(M2, M1, M2) -> M2\n");
  M2.display();

  AMatrix::add(MS2, MS1, MS2);
  message("AMatrix::add(MS2, MS1, MS2) -> MS2\n");
  MS2.display();

  // Product of a matrix by a constant
  reset(M1, M2, MS1, MS2, false);

  AMatrix::prod(M2, M1, 1.2);
  message("AMatrix::prod(M2, M1, 1.2) -> M2\n");
  M2.display();

  AMatrix::prod(MS2, MS1, 1.2);
  message("AMatrix::prod(MS2, MS1, 1.2) -> MS2\n");
  MS2.display();

  // Hadamard Product of two matrices
  reset(M1, M2, MS1, MS2, false);

  AMatrix::prodHadamard(M2, M1, M2);
  message("AMatrix::prodHadamard(M2, M1, M2) -> M2\n");
  M2.display();

  AMatrix::prodHadamard(MS2, MS1, MS2);
  message("AMatrix::prodHadamard(MS2, MS1, MS2) -> S2\n");
  MS2.display();

  //////////////////////////////////
  // Testing the member functions //
  //////////////////////////////////
  mestitle(0, "Testing member functions");

  // Adding a Constant
  reset(M1, M2, MS1, MS2, false);

  M1.addCst(10.);
  message("M1.addCst(10.) -> M1\n");
  M1.display();

  MS1.addCst(10.);
  message("MS1.addCst(10.) -> MS1\n");
  MS1.display();

  // Adding a Constant and allocating a new class
  reset(M1, M2, MS1, MS2, false);

  auto M1add = M1.addCst(10.);
  message("M1add = M1.addCst(10.) -> M1add\n");
  M1add.display();

  auto MS1add = MS1.addCst(10.);
  message("MS1add = MS1.addCst(10.) -> MS1add\n");
  MS1add.display();

  // Product of a matrix by a constant
  reset(M1, M2, MS1, MS2, false);

  M1.prodCst(2.);
  message("M1.prodCstInPlace(2.) -> M1\n");
  M1.display();

  MS1.prodCst(2.);
  message("MS1.prodCstInPlace(2.) -> MS1\n");
  MS1.display();

  // Product of a matrix by a constant and allocating a new class
  reset(M1, M2, MS1, MS2, false);

  auto M1prod = M1.prodCst(2.);
  message("M1prod = M1.prodCst(2.) -> M1prod\n");
  M1prod.display();

  auto MS1prod = MS1.prodCst(2.);
  message("MS1prod = MS1.prodCst(2.) -> MS1prod\n");
  MS1prod.display();

  ///////////////////////////
  // Testing the operators //
  ///////////////////////////
  mestitle(0, "Testing the operators");

  // Operator += with a constant
  reset(M1, M2, MS1, MS2, false);

  M1 += 10.;
  message("M1 += 10. -> M1\n");
  M1.display();

  MS1 += 10.;
  message("MS1 += 10. -> MS1\n");
  MS1.display();

  // Operator -= with a constant
  reset(M1, M2, MS1, MS2, false);

  M1 -= 10.;
  message("M1 -= 10. -> M1\n");
  M1.display();

  MS1 -= 10.;
  message("MS1 -= 10. -> MS1\n");
  MS1.display();

  // Operator + with a constant
  reset(M1, M2, MS1, MS2, false);

  auto M1addOp = M1 + 10.;
  message("M1addOp = M1 + 10. -> M1addOp\n");
  M1addOp.display();

  auto MS1addOp = MS1 + 10.;
  message("MS1addOp = MS1 + 10. -> MS1addOp\n");
  MS1addOp.display();

  // Operator + with a constant
  reset(M1, M2, MS1, MS2, false);

  auto M1subtractOp = M1 - 10.;
  message("M1subtractOp = M1 - 10. -> M1subtractOp\n");
  M1subtractOp.display();

  auto MS1subtractOp = MS1 - 10.;
  message("MS1subtractOp = MS1 - 10. -> MS1subtractOp\n");
  MS1subtractOp.display();

  ///////////////////////////
  // Testing child classes //
  ///////////////////////////
  mestitle(0, "Testing child classes");
  MatrixSquare MSq1(3);
  MatrixSquare MSq2(3);
  resetSquare(MSq1, MSq2, true);

  // Adding a constant
  resetSquare(MSq1, MSq2, false);
  AMatrix::add(MSq2, MSq1, 10.);
  message("add(MSq2, MSq1, 10.) -> MSq2\n");
  MSq2.display();

  // Adding a constant
  resetSquare(MSq1, MSq2, false);
  AMatrix::add(MSq2, MSq1, MSq2);
  message("add(MSq2, MSq1, MSq2) -> MSq2\n");
  MSq2.display();

  return 0;
}
