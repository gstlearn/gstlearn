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
#include "Basic/File.hpp"
#include "Basic/Law.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "LinearOp/CholeskySparse.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Matrix/NF_Triplet.hpp"

using namespace gstlrn;
MatrixSparse* _createSparseMatrix(Id n, double proba)
{
  // We create a square matrix
  NF_Triplet NF_T;
  for (Id icol = 0; icol < n; icol++)
    for (Id irow = 0; irow < n; irow++)
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

MatrixSymmetric* _createDenseMatrix(Id n, const MatrixSparse* Q)
{
  // Create the corresponding Symmetric matrix
  auto* M = new MatrixSymmetric(n);
  for (Id icol = 0; icol < n; icol++)
    for (Id irow = 0; irow < n; irow++)
    {
      double value = Q->getValue(irow, icol);
      M->setValue(irow, icol, value);
    }
  return M;
}

void printError(const String& name)
{
  message(">>> Function '%s' is INVALID =======================\n", name.c_str());
}

/****************************************************************************/
/*!
** Main Program for testing the Linear Algebra using Cholesky decomposition
** of sparse matrix
**
*****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  Id size            = 10;
  double proba       = 0.05;
  MatrixSparse* Q    = _createSparseMatrix(size, proba);
  MatrixSymmetric* M = _createDenseMatrix(size, Q);

  // Create a vector random gaussian values
  VectorDouble vecin = VH::simulateGaussian(size);
  VectorDouble vecout1(size);
  VectorDouble vecout2(size);
  VectorDouble vecout(size);
  VectorDouble vecref(size);

  // Creating the Cholesky objects
  CholeskySparse cholSparse(*Q);

  CholeskyDense cholDense(*M);

  // Checking the Cholesky decomposition
  AMatrix::productInPlace(vecref, *M, vecin);

  cholSparse.LtX(vecin, vecout);
  cholSparse.LX(vecout, vecout1);
  cholDense.LtX(vecin, vecout);
  cholDense.LX(vecout, vecout2);
  if (vecout1.isEqual(vecout2) && vecref.isEqual(vecout1))
    message(">>> Function 'LLt' is validated\n");
  else
  {
    printVector(vecref, "LLt (by Matrix)", true, true);
    printVector(vecout1, "LLt (by CholeskySparse)", true, true);
    printVector(vecout2, "LLt (by CholeskyDense)", true, true);
    printError("LLt");
  }

  // Checking method 'solve'
  (void)M->solve(vecin, vecref);

  cholSparse.solve(vecin, vecout1);
  cholDense.solve(vecin, vecout2);
  if (vecout1.isEqual(vecout2) && vecref.isEqual(vecout1))
    message(">>> Function 'solve' is validated\n");
  else
  {
    printVector(vecref, "Solve (by Matrix)", true, true);
    printVector(vecout1, "Solve (by CholeskySparse)", true, true);
    printVector(vecout2, "Solve (by CholeskyDense)", true, true);
    printError("solve");
  }

  // Checking method 'LX' followed by 'InvLX'
  cholSparse.LX(vecin, vecout);
  cholSparse.InvLX(vecout, vecout1);
  cholDense.LX(vecin, vecout);
  cholDense.InvLX(vecout, vecout2);

  if (vecout1.isEqual(vecout2) && vecout1.isEqual(vecin))
    message(">>> Function 'InvLX(LX)' is validated\n");
  else
  {
    printVector(vecin, "Function 'InvLX(LX)' (by Matrix)", true, true);
    printVector(vecout1, "Function 'InvLX(LX)' (by CholeskySparse)", true, true);
    printVector(vecout2, "Function 'InvLX(LX)' (by CholeskyDense)", true, true);
    printError("InvLX(LX)");
  }

  // Checking method 'InvLtX' followed by 'LtX'
  cholSparse.InvLtX(vecin, vecout);
  cholSparse.LtX(vecout, vecout1);
  cholDense.InvLtX(vecin, vecout);
  cholDense.LtX(vecout, vecout2);

  if (vecout1.isEqual(vecout2) && vecout1.isEqual(vecin))
    message(">>> Function 'LtX(InvLtX)' is validated\n");
  else
  {
    printVector(vecin, "Function 'LtX(InvLtX)' (by Matrix)", true, true);
    printVector(vecout1, "Function 'LtX(InvLtX)' (by CholeskySparse)", true, true);
    printVector(vecout2, "Function 'LtX(InvLtX)' (by CholeskyDense)", true, true);
    printError("LtX(InvLtX)");
  }

  // Checking the Stdev vector
  MatrixSymmetric MP(*M);
  (void)MP.invert();
  VectorDouble vecout1b = MP.getDiagonal();

  MatrixSparse* M2 = MatrixSparse::createFromTriplet(
    M->getMatrixToTriplet(), M->getNRows(), M->getNCols(), -1);
  CholeskySparse Qchol(*M2);
  MatrixSparse* proj = MatrixSparse::Identity(M->getNRows());
  Qchol.stdev(vecout2, proj, false);
  delete proj;

  if (vecout1b.isEqual(vecout2))
    message(">>> Function 'stdev' is validated\n");
  else
  {
    printVector(vecout1b, "Standard Deviation (by Matrix)", true, true);
    printVector(vecout2, "Standard Deviation (by Cholesky)", true, true);
    printError("stdev");
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
    printError("Log(Det)");
  }

  // Free the pointers
  delete Q;
  delete M;

  return (0);
}
