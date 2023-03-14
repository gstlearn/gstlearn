/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_old_f.h"

#include "LinearOp/Cholesky.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "csparse_f.h"

/****************************************************************************/
/*!
** Main Program for testing the Linear Algebra using Cholesky decomposition
** of sparse matrix
**
*****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
//  StdoutRedirect sr(sfn.str());
//
  int n = 10;

  // We create a square matrix (not necessarily sparse)

  double proba = 0.05;
  cs *Atriplet = cs_spalloc(0, 0, 1, 1, 1);
  for (int icol = 0; icol < n; icol++)
    for (int irow = 0; irow < n; irow++)
    {
      double value = law_gaussian();
      double tirage = law_uniform(0., 1.);
      if (icol != irow && tirage > proba) continue;
      (void) cs_entry(Atriplet, irow, icol, value);
    }
  cs *A = cs_triplet(Atriplet);
  cs_print_dim("Square Initial Matrix", A);
  Atriplet = cs_spfree(Atriplet);

  // The symmetric matrix is obtained as t(A) %*% A -> M is symmetric

  cs* At = cs_transpose(A, 1);
  cs* Q = cs_multiply(A, At);

  // Create a vector random gaussian values

  VectorDouble vecin = VH::simulateGaussian(n);
  VectorDouble vecout1(n);
  VectorDouble vecout2(n);

  // Create the corresponding Symmetric matrix

  MatrixSquareSymmetric M(n);
  for (int icol = 0; icol < n; icol++)
    for (int irow = 0; irow < n; irow++)
    {
      double value = cs_get_value(Q, irow, icol);
      M.setValue(irow, icol, value);
    }

  // Create the Cholesky object

  Cholesky Qchol(Q);
  Qchol.printout("Matrix used to demonstrate Cholesky Algebra", false);

  // Checking Product

  M.prodVector(vecin, vecout1);
  Qchol.evalDirect(vecin, vecout2);
  if (VH::isSame(vecout1,  vecout2))
    message("Product Mat %*% V is validated\n");
  else
  {
    VH::display("Product Mat %*% V (by Matrix)", vecout1);
    VH::display("Product Mat %*% V (by Cholesky)", vecout2);
  }

  // Checking Inverse

  (void) M.solve(vecin, vecout1);
  Qchol.evalInverse(vecin, vecout2);
  if (VH::isSame(vecout1,  vecout2))
    message("Product Mat^{-1} %*% V is validated\n");
  else
  {
    VH::display("Product Mat^{-1} %*% V (by Matrix)", vecout1);
    VH::display("Product Mat^{-1} %*% V (by Cholesky)", vecout2);
  }

  // Checking the Estimation of the Stedv vector

  (void) M.invert();
  VectorDouble vecout1b = M.getDiagonal();
  Qchol.stdev(vecout2);
  if (VH::isSame(vecout1b,  vecout2))
    message("Standard Deviation is validated\n");
  else
  {
    VH::display("Standard Deviation (by Matrix)", vecout1b);
    VH::display("Standard Deviation (by Cholesky)", vecout2);
  }

  // Free the pointers

  A  = cs_spfree(A);
  At = cs_spfree(At);
  Q  = cs_spfree(Q);
  return(0);
}
