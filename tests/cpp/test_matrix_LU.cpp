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
#include "geoslib_old_f.h"

#include "Matrix/MatrixFactory.hpp"
#include "Matrix/AMatrix.hpp"
#include "Matrix/MatrixInt.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixSquareDiagonal.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"

/****************************************************************************/
/*!
** Main Program for testing the new classes of matrix algebra
**
*****************************************************************************/
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  int neq = 3;
  int neq2 = neq * neq;
  MatrixSquareGeneral mat(neq);
  VectorDouble tab(neq2);

  MatrixSquareGeneral a(neq);
  a(0,0) = -1;
  a(0,1) =  0;
  a(0,2) =  3;
  a(1,0) = -2;
  a(1,1) = -2;
  a(1,2) =  7;
  a(2,0) = -5;
  a(2,1) =  0;
  a(2,2) = 20;
  a.display();
  MatrixSquareGeneral ai(a);

  // LU decomposition
  VectorDouble tl(neq2,0.);
  VectorDouble tu(neq2,0.);

  matrix_LU_decompose(neq, a.getValues().data(), tl.data(), tu.data());

  MatrixSquareGeneral atl(neq);
  atl.resetFromArray(neq, neq, tl.data());
  atl.display();

  MatrixSquareGeneral atu(neq);
  atu.resetFromArray(neq, neq, tu.data());
  atu.display();

  MatrixSquareGeneral res(neq);
  res.prodMatrix(atl, atu);
  message("\nChecking the product\n");
  res.display();
  message("compared to Initial\n");
  a.display();

  VectorDouble xtest(neq);
  VectorDouble x(neq);
  VectorDouble b = { 2., 7., 0.};
  VH::display("B",b);

  message("Inverse using LU\n");
  VectorDouble ais = a.getValues();
  (void) matrix_LU_invert(neq, ais.data());
  ai.resetFromArray(neq, neq, ais.data());
  ai.display();

  message("Inverse using invreal\n");
  ais = a.getValues();
  (void) matrix_invreal(ais.data(), neq);
  ai.resetFromArray(neq, neq, ais.data());
  ai.display();

  return(0);
}
