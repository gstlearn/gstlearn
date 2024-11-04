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
#include "LinearOp/CholeskySparse.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Matrix/NF_Triplet.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "Basic/Utilities.hpp"
#include "Matrix/LinkMatrixSparse.hpp"

MatrixSparse* _createSparseMatrix(int n, double proba)
{
  // We create a square matrix
  NF_Triplet NF_T;
  for (int icol = 0; icol < n; icol++)
    for (int irow = 0; irow < n; irow++)
    {
      double value  = law_gaussian();
      double tirage = law_uniform(0., 1.);
      if (icol != irow && tirage > proba) continue;
      NF_T.add(irow, icol, value);
    }
  MatrixSparse* A = MatrixSparse::createFromTriplet(NF_T);

  // The symmetric matrix is obtained as t(A) %*% A -> M is symmetric
  MatrixSparse* At = A->transpose();
  MatrixSparse* Q  = MatrixFactory::prodMatMat<MatrixSparse>(A, At);

  delete A;
  delete At;

  return Q;
}

MatrixSquareSymmetric* _createDenseMatrix(int n, const MatrixSparse* Q)
{
  // Create the corresponding Symmetric matrix
  MatrixSquareSymmetric* M = new MatrixSquareSymmetric(n);
  for (int icol = 0; icol < n; icol++)
    for (int irow = 0; irow < n; irow++)
    {
      double value = Q->getValue(irow, icol);
      M->setValue(irow, icol, value);
    }
  return M;
}

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

  int size                 = 10;
  double proba             = 0.05;
  MatrixSparse* Q          = _createSparseMatrix(size, proba);
  MatrixSquareSymmetric* M = _createDenseMatrix(size, Q);

  // Create a vector random gaussian values
  VectorDouble vecin = VH::simulateGaussian(size);
  VectorDouble vecout1(size);
  VectorDouble vecout2(size);
  VectorDouble vecout(size);
  VectorDouble vecref(size);

  // Creating the Cholesky objects
  CholeskySparse cholSparse(Q);
  CholeskyDense cholDense(M);

  // Checking the Cholesky decomposition
  M->prodMatVecInPlace(vecin, vecref);

  cholSparse.LtX(vecin, vecout);
  cholSparse.LX(vecout, vecout1);
  cholDense.LtX(vecin, vecout);
  cholDense.LX(vecout, vecout2);
  if (VH::isSame(vecout1, vecout2) && VH::isSame(vecref, vecout1))
    message(">>> Function 'LLt' is validated\n");
  else
  {
    VH::display("LLt (by Matrix)", vecref);
    VH::display("LLt (by CholeskySparse)", vecout1);
    VH::display("LLt (by CholeskyDense)", vecout2);
    message(">>> Function 'LLt' is INVALID =======================\n");
  }

  // Checking method 'solve'
  (void)M->solve(vecin, vecref);

  cholSparse.solve(vecin, vecout1);
  cholDense.solve(vecin, vecout2);
  if (VH::isSame(vecout1, vecout2) && VH::isSame(vecref, vecout1))
    message(">>> Function 'solve' is validated\n");
  else
  {
    VH::display("Solve (by Matrix)", vecref);
    VH::display("Solve (by CholeskySparse)", vecout1);
    VH::display("Solve (by CholeskyDense)", vecout2);
    message(">>> Function 'solve' is INVALID =======================\n");
  }

  // Checking method 'LX' followed by 'InvLX'
  cholSparse.LX(vecin, vecout);
  cholSparse.InvLX(vecout, vecout1);
  cholDense.LX(vecin, vecout);
  cholDense.InvLX(vecout, vecout2);

  if (VH::isSame(vecout1, vecout2) && VH::isSame(vecout1, vecin))
    message(">>> Function 'InvLX(LX)' is validated\n");
  else
  {
    VH::display("Function 'InvLX(LX)' (by Matrix)", vecin);
    VH::display("Function 'InvLX(LX)' (by CholeskySparse)", vecout1);
    VH::display("Function 'InvLX(LX)' (by CholeskyDense)", vecout2);
    message(">>> Function 'InvLX(LX)' is INVALID ========================\n");
  }

  // Checking method 'InvLtX' followed by 'LtX'
  cholSparse.InvLtX(vecin, vecout);
  cholSparse.LtX(vecout, vecout1);
  cholDense.InvLtX(vecin, vecout);
  cholDense.LtX(vecout, vecout2);

  if (VH::isSame(vecout1, vecout2) && VH::isSame(vecout1, vecin))
    message(">>> Function 'LtX(InvLtX)' is validated\n");
  else
  {
    VH::display("Function 'LtX(InvLtX)' (by Matrix)", vecin);
    VH::display("Function 'LtX(InvLtX)' (by CholeskySparse)", vecout1);
    VH::display("Function 'LtX(InvLtX)' (by CholeskyDense)", vecout2);
    message(">>> Function 'LtX(InvLtX)' is INVALID ========================\n");
  }

  // Checking the Stdev vector
  MatrixSquareSymmetric MP(*M);
  (void)MP.invert();
  VectorDouble vecout1b = MP.getDiagonal();

  // We use a Tim Davis sparse matrix cs as long as Qchol
  // stdev calculation is not available with eigen underlying matrix
  MatrixSparse* M2 = MatrixSparse::createFromTriplet(
    M->getMatrixToTriplet(), M->getNRows(), M->getNCols(), 0);
  CholeskySparse Qchol(M2);
  Qchol.stdev(vecout2, false);

  if (VH::isSame(vecout1b, vecout2))
    message(">>> Function 'stdev' is validated\n");
  else
  {
    VH::display("Standard Deviation (by Matrix)", vecout1b);
    VH::display("Standard Deviation (by Cholesky)", vecout2);
    message(">>> Function 'stdev' is INVALID =============================\n");
  }

  // Checking the calculation of Log(Det)
  double res1 = log(M->determinant());
  double res2 = Qchol.computeLogDeterminant();
  if (isZero(res1 - res2))
    message(">>> Log(Det) is validated\n");
  else
  {
    message("Log(Det) (by Matrix) = %lf\n", res1);
    message("Log(Det) (by Cholesky) = %lf\n", res2);
    message(">>> Log(Det) is INVALID =============================\n");
  }

  // Free the pointers
  delete Q;
  delete M;
  delete M2;

  return(0);
}
