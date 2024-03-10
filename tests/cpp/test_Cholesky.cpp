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

#include "LinearOp/Cholesky.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Matrix/NF_Triplet.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "Matrix/LinkMatrixSparse.hpp"

/****************************************************************************/
/*!
** Main Program for testing the Linear Algebra using Cholesky decomposition
** of sparse matrix
**
*****************************************************************************/
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  int n = 10;
  double proba = 0.05;

  // We create a square matrix (not necessarily sparse)
  NF_Triplet NF_T;
  for (int icol = 0; icol < n; icol++)
    for (int irow = 0; irow < n; irow++)
    {
      double value = law_gaussian();
      double tirage = law_uniform(0., 1.);
      if (icol != irow && tirage > proba) continue;
      NF_T.add(irow, icol, value);
    }
  MatrixSparse *A = MatrixSparse::createFromTriplet(NF_T);
  // The symmetric matrix is obtained as t(A) %*% A -> M is symmetric

  MatrixSparse* At = A->transpose();
  MatrixSparse* Q = MatrixFactory::prodMatMat<MatrixSparse>(A, At);

  // Create a vector random gaussian values
  VectorDouble vecin = VH::simulateGaussian(n);
  VectorDouble vecout1(n);
  VectorDouble vecout2(n);

  // Create the corresponding Symmetric matrix

  MatrixSquareSymmetric M(n);
  for (int icol = 0; icol < n; icol++)
    for (int irow = 0; irow < n; irow++)
    {
      double value = Q->getValue(irow, icol);
      M.setValue(irow, icol, value);
    }
  // Create the Cholesky object

  Cholesky Qchol(Q);
  message("Matrix used to demonstrate Cholesky Algebra\n");

  // Checking Product
  M.prodMatVecInPlace(vecin, vecout1);
  Qchol.evalDirect(vecin, vecout2);
  if (VH::isSame(vecout1,  vecout2))
  {
    message("Product Mat * V is validated\n");
  }
  else
  {
    VH::display("Product Mat * V (by Matrix)", vecout1);
    VH::display("Product Mat * V (by Cholesky)", vecout2);
  }

  // Checking Inverse
  (void) M.solve(vecin, vecout1);
  Qchol.evalInverse(vecin, vecout2);
  if (VH::isSame(vecout1,  vecout2))
    message("Product Mat^{-1} * V is validated\n");
  else
  {
    VH::display("Product Mat^{-1} * V (by Matrix)", vecout1);
    VH::display("Product Mat^{-1} * V (by Cholesky)", vecout2);
  }
  // Checking the Estimation of the Stdev vector

  MatrixSquareSymmetric MP(M);
  (void) MP.invert();
  VectorDouble vecout1b = MP.getDiagonal();
  Qchol.stdev(vecout2);
  if (VH::isSame(vecout1b,  vecout2))
    message("Standard Deviation is validated\n");
  else
  {
    VH::display("Standard Deviation (by Matrix)", vecout1b);
    VH::display("Standard Deviation (by Cholesky)", vecout2);
  }
  // Checking the calculation of Log(Det)

  double res1 = log(M.determinant());
  double res2 = Qchol.getLogDeterminant();
  if (ABS(res1 - res2) < EPSILON10)
    message("Log(Det) is validated\n");
  else
  {
    message("Log(Det) (by Matrix) = %lf\n", res1);
    message("Log(Det) (by Cholesky) = %lf\n", res2);
  }
  // Free the pointers

  delete A;
  delete At;
  delete Q;
  return(0);
}
