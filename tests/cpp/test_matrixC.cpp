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
#include "geoslib_d.h"
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "Matrix/AMatrix.hpp"
#include "Matrix/MatrixInt.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixSquareDiagonal.hpp"
#include "Matrix/MatrixSquareDiagonalCst.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Basic/Vector.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"

void reset_to_initial_contents(AMatrix* M,
                               MatrixSquareDiagonalCst& D,
                               MatrixRectangular& MRR,
                               MatrixSquareGeneral& MSG,
                               MatrixSquareSymmetric& MSS,
                               AMatrix* MSP,
                               MatrixSquareDiagonal& MSD,
                               MatrixSquareDiagonalCst& MSC)
{
  MRR.setValues (M->getValues());
  MSG.setValues (M->getValues());
  MSS.setValues (M->getValues());
  MSP->setValues(M->getValues());

  MSD.setValues(D.getValues());
  MSC.setValues(D.getValues());
}

/****************************************************************************/
/*!
** Main Program for testing the new classes of matrix algebra
**
*****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  message("Cloning Matrix of integers\n");
  MatrixInt mati(2,3);
  mati.setValues({1, 2, 3, 4, 5, 6});
  mati.display();
  MatrixInt* mati2(mati.clone()); // dynamic_cast no more needed
  // equivalent to
  // MatrixInt* mati2 = mati.clone();
  mati2->display();

  message("Cloning Matrix of doubles\n");
  MatrixRectangular matd(2,3);
  matd.setValues({1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
  matd.display();
  AMatrix* pmat = &matd;
  MatrixRectangular* matd2(dynamic_cast<MatrixRectangular*>(pmat->clone())); // dynamic_cast cannot be avoided here
  matd2->display();

  VectorDouble V1,V2,V3,Vref;
  law_set_random_seed(32432);
  int nrow = 7; // For these tests, the matrix MUST be square (ncol = nrow)
  int ncol = 7;
  double proba = 0.4; // Probability to set values to 0 (making matrix sparse)

  // We create a square symmetrical matrix (not necessarily sparse)
  
  MatrixRectangular MR(nrow, ncol);
  for (int icol = 0; icol < ncol; icol++)
    for (int irow = 0; irow < nrow; irow++)
    {
      double value = law_gaussian();
      double tirage = law_uniform(0., 1.);
      if (tirage < proba) value = 0.;
      MR.setValue(irow, icol, value);
    }

  // The symmetric matrix is obtained as t(MR) %*% MR -> M is symmetric

  AMatrix* MRt = MR.transpose(); // Using cloneable feature
  
  // Equivalent instruction using shortcut function
  //AMatrix* MRt = transpose(MR);

  // Still equivalent but in two lines
  //AMatrix* MRt = dynamic_cast<AMAtrix*>(MR.clone());
  //MRt->transposeInPlace();

  // Still equivalent but with no more pointer
  //MatrixRectangular MRt = MR;
  //MRt.transposeInPlace();

  AMatrix* M = prodMatrix(MRt, &MR);

  // Equivalent but more lengthy (but no more pointer)
  //MatrixRectangular M;
  //M.prodMatrix(*MRt, MR);

  message("Matrix M\n");
  M->display();

  // Creating two vectors for future use

  V1.resize(ncol,0.);
  V2.resize(ncol,0.);

  // Create the different matrix formats (by conversion or extraction)
 
  // To a rectangular matrix
  MatrixRectangular MRR(nrow,ncol);
  MRR.setValues(M->getValues());
  message("Matrix MRR\n");
  MRR.display();

  // To a square general matrix
  MatrixSquareGeneral MSG(*M);
  message("Matrix MSG\n");
  MSG.display();

  // To a square symmetric matrix
  MatrixSquareSymmetric MSS(*M);
  message("Matrix MSS\n");
  MSS.display();

  // To a sparse matrix
  AMatrix* MSP = M->toSparse();
  message("Matrix MSP\n");
  MSP->display();

  // Creating a Diagonal matrix (from M)
  double cst = law_gaussian();

  MatrixSquareDiagonal MSD(nrow);
  for (int irow=0; irow<nrow; irow++)
    MSD(irow,irow) = cst;
  message("Matrix MSD\n");
  MSD.display();

  // Creating a Constant Diagonal Matrix

  MatrixSquareDiagonalCst MSC, D;
  D.reset(nrow,ncol,cst);
  MSC.reset(nrow,ncol,cst);
  message("Matrix MSC\n");
  MSC.display();

  /**
   * Adding a constant to the diagonal of a matrix
   */
  double addendum = 1.432;

  mestitle(0,"Adding a constant value to the diagonal of a matrix");
  reset_to_initial_contents(M, D, MRR, MSG, MSS, MSP, MSD, MSC);

  MRR.addScalarDiag(addendum);
  MSG.addScalarDiag(addendum);
  message("Are results for MRR and MSG similar: %d\n",MRR.isSame(MSG));
  MSS.addScalarDiag(addendum);
  message("Are results for MRR and MSS similar: %d\n",MRR.isSame(MSS));
  MSP->addScalarDiag(addendum);
//  message("Are results for MRR and MSP similar: %d\n",MRR.isSame(*MSP));

  MSD.addScalarDiag(addendum);
  MSC.addScalarDiag(addendum);
  message("Are results for MSD and MSC similar: %d\n",MSD.isSame(MSC));

  /**
   * Multiplying the matrix by a constant
   */
  double multiply = 3.2;

  mestitle(0,"Multiplying a Matrix by a constant");
  reset_to_initial_contents(M, D, MRR, MSG, MSS, MSP, MSD, MSC);

  MRR.prodScalar(multiply);
  MSG.prodScalar(multiply);
  message("Are results for MRR and MSG similar: %d\n",MRR.isSame(MSG));
  MSS.prodScalar(multiply);
  message("Are results for MRR and MSS similar: %d\n",MRR.isSame(MSS));
  MSP->prodScalar(multiply);
//  message("Are results for MRR and MSP similar: %d\n",MRR.isSame(*MSP));

  MSD.prodScalar(multiply);
  MSC.prodScalar(multiply);
  message("Are results for MSD and MSC similar: %d\n",MSD.isSame(MSC));

  /**
   * Adding a constant to a matrix
   * Note: This does not make sense for sparse or diagonal matrices
   */

  mestitle(0,"Adding a constant value to the whole matrix");
  reset_to_initial_contents(M, D, MRR, MSG, MSS, MSP, MSD, MSC);

  MRR.addScalar(addendum);
  MSG.addScalar(addendum);
  message("Are results for MRR and MSG similar: %d\n",MRR.isSame(MSG));
  MSS.addScalar(addendum);
  message("Are results for MRR and MSS similar: %d\n",MRR.isSame(MSS));

  /**
    * Linear combination
    */
  mestitle(0,"Linear combination of matrices");
  reset_to_initial_contents(M, D, MRR, MSG, MSS, MSP, MSD, MSC);

  double cx =  1.3;
  double cy = -0.3;

  MRR.linearCombination(cx,cy,MRR);
  MSG.linearCombination(cx,cy,MSG);
  message("Are results for MRR and MSG similar: %d\n",MRR.isSame(MSG));
  MSS.linearCombination(cx,cy,MSS);
  message("Are results for MRR and MSS similar: %d\n",MRR.isSame(MSS));
  MSP->linearCombination(cx,cy,*MSP);
//  message("Are results for MRR and MSP similar: %d\n",MRR.isSame(*MSP));

  MSD.linearCombination(cx,cy,MSD);
  MSC.linearCombination(cx,cy,MSC);
  message("Are results for MSD and MSC similar: %d\n",MSD.isSame(MSC));

  /**
   * Extraction of a Vector
   * All the tests are not performed on all the matrix types
   */
  mestitle(0,"Extracting Vectors from Matrix");
  reset_to_initial_contents(M, D, MRR, MSG, MSS, MSP, MSD, MSC);

  message("MRR and MSP matrices are used as Reference\n");
  MRR.display();
  Vref = MRR.getDiagonal();
  V1 = MSP->getDiagonal();
  print_vector("Main Diagonal",0,(int) Vref.size(),Vref.data());
  message("Are results for MRR and MSP similar: %d\n",ut_vector_same(Vref,V1));
  Vref = MRR.getDiagonal(1);
  V1 = MSP->getDiagonal(1);
  print_vector("Second Diagonal Below",0,(int) Vref.size(),Vref.data());
  message("Are results for MRR and MSP similar: %d\n",ut_vector_same(Vref,V1));
  Vref = MRR.getDiagonal(-2);
  V1 = MSP->getDiagonal(-2);
  print_vector("Third Diagonal Above",0,(int) Vref.size(),Vref.data());
  message("Are results for MRR and MSP similar: %d\n",ut_vector_same(Vref,V1));
  Vref = MRR.getRow(2);
  V1 = MSP->getRow(2);
  print_vector("Third Row",0,(int) Vref.size(),Vref.data());
  message("Are results for MRR and MSP similar: %d\n",ut_vector_same(Vref,V1));
  Vref = MRR.getColumn(3);
  V1 = MSP->getColumn(3);
  print_vector("Fourth Column",0,(int) Vref.size(),Vref.data());
  message("Are results for MRR and MSP similar: %d\n",ut_vector_same(Vref,V1));

  /**
   * Product of the matrix by a vector
   */

  mestitle(0,"Product of the matrix by a vector");
  reset_to_initial_contents(M, D, MRR, MSG, MSS, MSP, MSD, MSC);

  Vref.resize(nrow,0.);
  MRR.prodVector(V1, Vref);
  MSG.prodVector(V1, V2);
  message("Are results for MRR and MSG similar: %d\n",ut_vector_same(Vref,V2));
  MSS.prodVector(V1, V2);
  message("Are results for MRR and MSS similar: %d\n",ut_vector_same(Vref,V2));
  MSP->prodVector(V1, V2);
  message("Are results for MRR and MSP similar: %d\n",ut_vector_same(Vref,V2));

  MSD.prodVector(V1, Vref);
  MSC.prodVector(V1, V2);
  message("Are results for MSD and MSC similar: %d\n",ut_vector_same(Vref,V2));

  /**
   * Linear solver
   */

  mestitle(0,"Matrix Linear Solver");
  reset_to_initial_contents(M, D, MRR, MSG, MSS, MSP, MSD, MSC);
  V3.resize(nrow,0.);
  message("Solve X from A*X=B. Compute A*X and compare with B\n");

  MSS.solve(V1, V2);
  MSS.prodVector(V2, V3);
  message("Are results correct for MSS: %d\n",ut_vector_same(V1,V3));
  MSP->solve(V1, V2);
  MSP->prodVector(V2, V3);
  message("Are results correct for MSP: %d\n",ut_vector_same(V1,V3));

  MSD.solve(V1, V2);
  MSD.prodVector(V2, V3);
  message("Are results correct for MSD: %d\n",ut_vector_same(V1,V3));
  MSC.solve(V1, V2);
  MSC.prodVector(V2, V3);
  message("Are results correct for MSC: %d\n",ut_vector_same(V1,V3));

  /**
   * Inversion
   */

  mestitle(0,"Matrix Inversion");
  reset_to_initial_contents(M, D, MRR, MSG, MSS, MSP, MSD, MSC);
  message("Calculate B=A^{-1}. Compute A*B and compare to Identity\n");

  AMatrix* Res;
  MatrixSquareGeneral MSGref = MSG; // Used to perform A*A-1 and check Identity

  MSG.invert();
  Res = prodMatrix(&MSG, &MSGref);
  message("Are results correct for MSG: %d\n",Res->isIdentity());
  delete Res;

  MSS.invert();
  Res = prodMatrix(&MSS, &MSGref);
  message("Are results correct for MSS: %d\n",Res->isIdentity());
  delete Res;

  MSP->invert();
  Res = prodMatrix(MSP, &MSGref);
  message("Are results correct for MSP: %d\n",Res->isIdentity());
  delete Res;

  MatrixSquareDiagonal MSDref = MSD; // Used to perform A*A-1 and check Identity

  MSD.invert();
  Res = prodMatrix(&MSD, &MSDref);
  message("Are results correct for MSD: %d\n",Res->isIdentity());
  delete Res;

  MSC.invert();
  Res = prodMatrix(&MSC, &MSDref);
  message("Are results correct for MSC: %d\n",Res->isIdentity());
  delete Res;

  // Free the pointers

  delete M;
  delete MSP;
  return(0);
}
