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
#include "Matrix/AMatrix.hpp"
#include "Matrix/MatrixInt.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixSquare.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Matrix/NF_Triplet.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "LinearOp/CholeskySparse.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "Basic/OptCst.hpp"

void st_invgen()
{
  MatrixSymmetric aaa(4);
  MatrixSymmetric bbb(4);

  aaa.setValue(0, 0, 2.);
  aaa.setValue(1, 0, 1.);
  aaa.setValue(0, 1, 1.);
  aaa.setValue(1, 1, 4.);
  aaa.display();

  (void) aaa.computeGeneralizedInverse(bbb);

  message("Generalized Inverse\n");
  bbb.display();
}

void reset_to_initial_contents(AMatrix* M,
                               MatrixDense& MRR,
                               MatrixSquare& MSG,
                               MatrixSymmetric& MSS,
                               MatrixSparse* MSP)
{
  MRR.setValues(M->getValues());
  MSG.setValues(M->getValues());
  MSS.setValues(M->getValues());
  MSP->setValues(M->getValues());
}

/****************************************************************************/
/*!
** Main program for testing the new classes of matrix algebra
**
*****************************************************************************/
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  setMultiThread(8);

  setGlobalFlagEigen(true); // Check the Eigen version or not. Essential for first part.
  OptCst::define(ECst::NTCOL, -1);
  OptCst::define(ECst::NTROW, -1);

  // Checking the inverse generalized matrix
  message("Checking the inverse generalized matrix\n");
  st_invgen();

  message("Cloning Matrix of integers\n");
  MatrixInt mati(2,3);
  mati.setValues({1, 2, 3, 4, 5, 6});
  mati.display();
  MatrixInt* mati2(mati.clone());
  mati2->display();

  // Checking assessor for MatrixInt
  int imemo = mati(1,2);
  message("Initial value of mati(1,2) = %d\n", imemo);
  mati(1,2) = -23;
  mati.display();
  mati(1,2) = imemo; // set back to initial value

  message("Cloning Matrix of doubles\n");
  MatrixDense matd(2,3);
  matd.setValues({1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
  matd.display();
  AMatrix* pmat = &matd;
  MatrixDense* matd2(dynamic_cast<MatrixDense*>(pmat->clone())); // dynamic_cast cannot be avoided here
  matd2->display();

  VectorDouble V1,V2,V3,Vref;
  law_set_random_seed(32432);
  int nrow = 7; // For these tests, the matrix MUST be square (ncol = nrow)
  int ncol = 7;
  double proba = 0.4; // Probability to set values to 0 (making matrix sparse)

  // We create a square symmetrical matrix (not necessarily sparse)
  
  MatrixDense MR(nrow, ncol);
  for (int icol = 0; icol < ncol; icol++)
    for (int irow = 0; irow < nrow; irow++)
    {
      double value = law_gaussian();
      double tirage = law_uniform(0., 1.);
      if (tirage < proba) value = 0.;
      MR.setValue(irow, icol, value);
    }
  message("Matrix MR\n");
  MR.display();

  // Checking using the operator to modify and correct the initial matrix MR()

  double memo = MR(1,2);
  message("Initial value of M(1,2) = %lf\n", memo);

  double new_value = 111.111;
  MR(1,2) = new_value;
  message("Modifying it to new value = %lf\n", new_value);
  MR.display();

  MR(1,2) = memo; // Set back to initial value
  message("Resetting to initial value\n");
  MR.display();

  // The symmetric matrix is obtained as M = t(MR) %*% MR

  AMatrix* MRt = MR.transpose();
  MRt->display();
//  AMatrix* M = prodMatMatInPlace(MRt, &MR);
  AMatrix* M = MatrixFactory::prodMatMat(MRt, &MR);
  message("Matrix M (should be symmetric). Checking = %d\n", (int) M->isSymmetric());
  M->display();

  // Creating two vectors for future use

  V1.resize(ncol,0.);
  V2.resize(ncol,0.);

  // Create the different matrix formats (by conversion or extraction)
 
  // To a rectangular matrix
  MatrixDense MRR(nrow,ncol);
  MRR.setValues(M->getValues());
  message("Matrix MRR\n");
  MRR.display();

  // To a square general matrix
  MatrixSquare MSG(*M);
  message("Matrix MSG\n");
  MSG.display();

  // To a square symmetric matrix
  MatrixSymmetric MSS(*M);
  message("Matrix MSS\n");
  MSS.display();

  // To a sparse matrix
  MatrixSparse* MSP = createFromAnyMatrix(M);
  message("Matrix MSP\n");
  MSP->display();

  /**
   * Adding a constant value to the diagonal of a matrix
   */
  double addendum = 1.432;

  mestitle(0,"Adding a constant value to the diagonal of a matrix");
  reset_to_initial_contents(M, MRR, MSG, MSS, MSP);
  message("Reference MRR (before addition)\n");
  MRR.display();
  MRR.addScalarDiag(addendum);
  message("Reference MRR (after addition)\n");
  MRR.display();

  MSG.addScalarDiag(addendum);
  message("Are results for MRR and MSG similar: %d\n",(int) MRR.isSame(MSG));
  MSS.addScalarDiag(addendum);
  message("Are results for MRR and MSS similar: %d\n",(int) MRR.isSame(MSS));
  MSP->addScalarDiag(addendum);
  message("Are results for MRR and MSP similar: %d\n",(int) MRR.isSame(*MSP));

  /**
   * Multiplying the matrix by a constant
   */
  double multiply = 3.2;

  mestitle(0,"Multiplying a Matrix by a constant");
  reset_to_initial_contents(M, MRR, MSG, MSS, MSP);
  message("Reference MRR (before multiplication)\n");
  MRR.display();
  MRR.prodScalar(multiply);
  message("Reference MRR (after multiplication)\n");
  MRR.display();

  MSG.prodScalar(multiply);
  message("Are results for MRR and MSG similar: %d\n",(int) MRR.isSame(MSG));
  MSS.prodScalar(multiply);
  message("Are results for MRR and MSS similar: %d\n",(int) MRR.isSame(MSS));
  MSP->prodScalar(multiply);
  message("Are results for MRR and MSP similar: %d\n",(int) MRR.isSame(*MSP));

  /**
   * Adding a constant to a matrix
   * Note: This does not make sense for sparse or diagonal matrices
   */

  mestitle(0,"Adding a constant value to the whole matrix");
  reset_to_initial_contents(M, MRR, MSG, MSS, MSP);
  message("Reference MRR (before addition)\n");
  MRR.display();
  MRR.addScalar(addendum);
  message("Reference MRR (after addition)\n");
  MRR.display();

  MSG.addScalar(addendum);
  message("Are results for MRR and MSG similar: %d\n",(int) MRR.isSame(MSG));
  MSS.addScalar(addendum);
  message("Are results for MRR and MSS similar: %d\n",(int) MRR.isSame(MSS));

  /**
    * Linear combination
    */
  double cx =  1.3;
  double cy = -0.5;

  mestitle(0,"Linear combination of matrices");
  reset_to_initial_contents(M, MRR, MSG, MSS, MSP);
  message("Reference MRR (before linear combination)\n");
  MRR.display();
  MRR.addMatInPlace(MRR, cx, cy);
  message("Reference MRR (after linear combination)\n");
  MRR.display();

  MSG.addMatInPlace(MSG, cx, cy);
  message("Are results for MRR and MSG similar: %d\n",(int) MRR.isSame(MSG));
  MSS.addMatInPlace(MSS, cx, cy);
  message("Are results for MRR and MSS similar: %d\n",(int) MRR.isSame(MSS));
  MSP->addMatInPlace(*MSP, cx, cy);
  message("Are results for MRR and MSP similar: %d\n",(int) MRR.isSame(*MSP));

  /**
   * Extraction of a Vector
   * All the tests are not performed on all the matrix types
   */
  mestitle(0,"Extracting Vectors from Matrix");
  reset_to_initial_contents(M, MRR, MSG, MSS, MSP);

  message("MRR and MSP matrices are used as Reference\n");
  MRR.display();
  Vref = MRR.getDiagonal();
  VH::dump("Reference Vector", Vref);

  V1 = MSP->getDiagonal();
  print_vector("Main Diagonal",0,(int) Vref.size(),Vref.data());
  message("Are results for MRR and MSP similar: %d\n",(int) VH::isEqual(Vref,V1));
  Vref = MRR.getDiagonal(1);
  V1 = MSP->getDiagonal(1);
  print_vector("Second Diagonal Below",0,(int) Vref.size(),Vref.data());
  message("Are results for MRR and MSP similar: %d\n",(int) VH::isEqual(Vref,V1));
  Vref = MRR.getDiagonal(-2);
  V1 = MSP->getDiagonal(-2);
  print_vector("Third Diagonal Above",0,(int) Vref.size(),Vref.data());
  message("Are results for MRR and MSP similar: %d\n",(int) VH::isEqual(Vref,V1));
  Vref = MRR.getRow(2);
  V1 = MSP->getRow(2);
  print_vector("Third Row",0,(int) Vref.size(),Vref.data());
  message("Are results for MRR and MSP similar: %d\n",(int) VH::isEqual(Vref,V1));
  Vref = MRR.getColumn(3);
  V1 = MSP->getColumn(3);
  print_vector("Fourth Column",0,(int) Vref.size(),Vref.data());
  message("Are results for MRR and MSP similar: %d\n",(int) VH::isEqual(Vref,V1));

  /**
   * Product of the matrix by a vector
   */
  Vref.resize(nrow,0.);

  mestitle(0,"Product of the matrix by a vector");
  reset_to_initial_contents(M, MRR, MSG, MSS, MSP);
  message("Reference Matrix\n");
  MRR.display();
  VH::dump("Reference Input Vector",V1);
  MRR.prodMatVecInPlace(V1, Vref);
  VH::dump("Reference Output Vector",Vref);

  MSG.prodMatVecInPlace(V1, V2);
  message("Are results for MRR and MSG similar: %d\n",(int) VH::isEqual(Vref,V2));
  MSS.prodMatVecInPlace(V1, V2);
  message("Are results for MRR and MSS similar: %d\n",(int) VH::isEqual(Vref,V2));
  MSP->prodMatVecInPlace(V1, V2);
  message("Are results for MRR and MSP similar: %d\n",(int) VH::isEqual(Vref,V2));

  /**
   * Linear solver
   */

  V3.resize(nrow,0.);

  mestitle(0,"Matrix Linear Solver");
  reset_to_initial_contents(M, MRR, MSG, MSS, MSP);
  message("Solve X from A*X=B. Compute A*X and compare with B\n");

  message("Reference Matrix\n");
  MSS.display();
  VH::dump("Reference Input Vector",V1);
  MSS.solve(V1, V2);
  VH::dump("Reference Output Vector",V2);

  MSS.prodMatVecInPlace(V2, V3);
  message("Are results correct for MSS: %d\n",(int) VH::isEqual(V1,V3));
  MSP->solve(V1, V2);
  MSP->prodMatVecInPlace(V2, V3);
  message("Are results correct for MSP: %d\n",(int) VH::isEqual(V1,V3));

  /**
   * Inversion
   */

  mestitle(0,"Matrix Inversion");
  reset_to_initial_contents(M, MRR, MSG, MSS, MSP);
  message("Calculate B=A^{-1}. Compute A*B and compare to Identity\n");

  AMatrix* Res;
  MatrixSquare MSGref = MSG; // Used to perform A*A-1 and check Identity
  message("Reference Matrix\n");
  MSGref.display();
  MSG.invert();
  message("Inverse Matrix\n");
  MSG.display();

  Res = MatrixFactory::prodMatMat(&MSG, &MSGref);
  message("Are results correct for MSG: %d\n",(int) Res->isIdentity());
  delete Res;

  MSS.invert();
  Res = MatrixFactory::prodMatMat(&MSS, &MSGref);
  message("Are results correct for MSS: %d\n",(int) Res->isIdentity());
  delete Res;

  MSP->invert();
  Res = MatrixFactory::prodMatMat(MSP, &MSGref);
  message("Are results correct for MSP: %d\n",(int) Res->isIdentity());
  delete Res;

  /*
   * Auxiliary functions (virtual in AMatrix)
   */

  mestitle(0,"Setting various Elements to known values");
  reset_to_initial_contents(M, MRR, MSG, MSS, MSP);
  nrow = MSG.getNRows();
  ncol = MSG.getNCols();

  // For dense matrix: we use the square matrix as a reference in order to allow all possible functions
  message("Reference Square matrix\n");
  MSG.display();

  int icol0 = 1;
  message("Setting Column (%d) to a vector (sequence from 1 to %d)\n", icol0, nrow);
  VectorDouble myCol = VH::sequence(1., (double) nrow);
  MSG.setColumn(icol0, myCol);
  MSG.display();

  int irow0 = 2;
  message("Setting Row (%d) to a vector (sequence from 1 to %d)\n", irow0, ncol);
  VectorDouble myRow = VH::sequence(1., (double) ncol);
  MSG.setRow(irow0, myRow);
  MSG.display();

  double vadd0 = 12.;
  message("Adding constant %lf to all terms of matrix\n", vadd0);
  MSG.addScalar(vadd0);
  MSG.display();

  double vadddiag0 = -23.;
  message("Adding constant %lf to diagonal terms of matrix\n", vadddiag0);
  MSG.addScalarDiag(vadddiag0);
  MSG.display();

  double vprod0 = 1.2;
  message("Product of all terms of matrix by constant %lf\n", vprod0);
  MSG.prodScalar(vprod0);
  MSG.display();

  message("Adding the matrix to itself\n");
  MSG.addMatInPlace(MSG);
  MSG.display();

  cx = 1.2;
  cy = -2.3;
  message("Making the linear combination of the matrix (multiplied by %f) and itself (multiplied by %lf)\n", cx, cy);
  MatrixSquare MSG3(MSG);
  MSG.addMatInPlace(MSG3, cx, cy);
  MSG.display();

  message("Multiplying current matrix column-wise by a vector (sequence)\n");
  myCol = VH::sequence(1., (double) nrow);
  MSG.multiplyColumn(myCol);
  MSG.display();

  message("Dividing current matrix column-wise by a vector (sequence)\n");
  myCol = VH::sequence(1., (double) nrow);
  MSG.divideColumn(myCol);
  MSG.display();

  message("Multiplying current matrix row-wise by a vector (sequence)\n");
  myRow = VH::sequence(1., (double) ncol);
  MSG.multiplyRow(myRow);
  MSG.display();

  message("Dividing current matrix row-wise by a vector (sequence)\n");
  myRow = VH::sequence(1., (double) ncol);
  MSG.divideRow(myRow);
  MSG.display();

  VectorDouble myRowRes;
  VectorDouble myColRes;

  message("Multiplying sequence vector by matrix\n");
  myCol = VH::sequence(1., (double) nrow);
  myRowRes = MSG.prodVecMat(myCol, false);
  VH::dump("Resulting Vector", myRowRes);

  message("Multiplying matrix (transposed) by sequence vector\n");
  myCol = VH::sequence(1., (double) nrow);
  myRowRes = MSG.prodMatVec(myCol, true);
  VH::dump("Resulting Vector", myRowRes);

  message("Multiplying matrix by sequence vector\n");
  myRow = VH::sequence(1., (double) ncol);
  myColRes = MSG.prodMatVec(myRow, false);
  VH::dump("Resulting Vector", myColRes);

  message("Multiplying sequence vector by matrix (transposed)\n");
  myRow = VH::sequence(1., (double) ncol);
  myColRes = MSG.prodVecMat(myRow, true);
  VH::dump("Resulting Vector", myColRes);

  message("Making the product of the matrix by itself\n");
  MatrixSquare MSG2(MSG);
  MSG.prodMatMatInPlace(&MSG2, &MSG2);
  MSG.display();

  message("Clearing matrix and Setting Diagonal to a vector (sequence from 1 to %d)\n", ncol);
  VectorDouble myDiag = VH::sequence(1., (double) ncol);
  MSG.setDiagonal(myDiag);
  MSG.display();

  double vdiag0 = -4.;
  message("Clearing matrix and Setting Diagonal to %lf\n", vdiag0);
  MSG.setDiagonalToConstant(vdiag0);
  MSG.display();

  // For sparse matrix
  message("Reference Sparse matrix\n");
  setUpdateNonZeroValue(0); // Allow flexible update of sparse matrix
  MSP->display();

  message("Setting terms of Column (%d) to a vector (sequence from 1 to %d)\n", icol0, nrow);
  myCol = VH::sequence(1., (double) nrow);
  MSP->setColumn(icol0, myCol);
  MSP->display();

  message("Setting terms of Row (%d) to a vector (sequence from 1 to %d)\n", irow0, ncol);
  myRow = VH::sequence(1., (double) ncol);
  MSP->setRow(irow0, myRow);
  MSP->display();

  message("Adding constant %lf to all non-zero terms of matrix\n", vadd0);
  MSP->addScalar(vadd0);
  MSP->display();

  message("Adding constant %lf to diagonal non-zero terms of matrix\n", vadddiag0);
  MSP->addScalarDiag(vadddiag0);
  MSP->display();

  message("Product of all non-zero terms of matrix by constant %lf\n", vprod0);
  MSP->prodScalar(vprod0);
  MSP->display();

  message("Adding the matrix to itself\n");
  MSP->addMatInPlace(*MSP);
  MSP->display();

  message("Making the linear combination of the matrix (multiplied by %f) and itself (multiplied by %lf)\n", cx, cy);
  MatrixSparse MSP3(*MSP);
  MSP->addMatInPlace(MSP3, cx, cy);
  MSP->display();

  message("Multiplying current matrix column-wise by a vector (sequence)\n");
  myCol = VH::sequence(1., (double) nrow);
  MSP->multiplyColumn(myCol);
  MSP->display();

  message("Dividing current matrix column-wise by a vector (sequence)\n");
  myCol = VH::sequence(1., (double) nrow);
  MSP->divideColumn(myCol);
  MSP->display();

  message("Multiplying current matrix row-wise by a vector (sequence)\n");
  myRow = VH::sequence(1., (double) ncol);
  MSP->multiplyRow(myRow);
  MSP->display();

  message("Dividing current matrix row-wise by a vector (sequence)\n");
  myRow = VH::sequence(1., (double) ncol);
  MSP->divideRow(myRow);
  MSP->display();

  message("Multiplying sequence vector by matrix\n");
  myCol = VH::sequence(1., (double) nrow);
  myRowRes = MSP->prodVecMat(myCol, false);
  VH::dump("Resulting Vector", myRowRes);

  message("Multiplying matrix (transposed) by sequence vector\n");
  myCol = VH::sequence(1., (double) nrow);
  myRowRes = MSP->prodMatVec(myCol, true);
  VH::dump("Resulting Vector", myRowRes);

  message("Multiplying matrix by sequence vector\n");
  myRow = VH::sequence(1., (double) ncol);
  myColRes = MSP->prodMatVec(myRow, false);
  VH::dump("Resulting Vector", myColRes);

  message("Multiplying sequence vector by matrix (transposed)\n");
  myRow = VH::sequence(1., (double) ncol);
  myColRes = MSP->prodVecMat(myRow, true);
  VH::dump("Resulting Vector", myColRes);

  message("Making the product of the matrix by itself\n");
  MatrixSparse MSP2(*MSP);
  MSP->prodMatMatInPlace(&MSP2, &MSP2);
  MSP->display();

  message("Clearing matrix and Setting Diagonal to a vector (sequence from 1 to %d)\n", ncol);
  myDiag = VH::sequence(1., (double) ncol);
  MSP->setDiagonal(myDiag);
  MSP->display();

  message("Clearing matrix and Setting Diagonal to %lf\n", vdiag0);
  MSP->setDiagonalToConstant(vdiag0);
  MSP->display();

  /*
   * Testing LU
   */

  int neq = 3;
  int neq2 = neq * neq;
  MatrixSquare mat(neq);
  VectorDouble tab(neq2);

  MatrixSquare a(neq);
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
  MatrixSquare ai(a);

  // LU decomposition
  MatrixSquare tl(neq);
  MatrixSquare tu(neq);

  a.decomposeLU(tl, tu);

  tl.display();
  tu.display();

  MatrixSquare res(neq);
  res.prodMatMatInPlace(&tl, &tu);
  message("\nChecking the product\n");
  res.display();
  message("compared to Initial\n");
  a.display();

  VectorDouble xtest(neq);
  VectorDouble x(neq);
  VectorDouble b = { 2., 7., 0.};
  VH::dump("B",b);

  message("Inverse (using LU or invreal depending on the dimension)\n");
  (void) a.invert();
  a.display();

  // ************
  // Eigen values
  // ************

  mestitle(0,"Eigen values calculation for Dense matrices");
  reset_to_initial_contents(M, MRR, MSG, MSS, MSP);

  // Get a Dense matrix
  VectorDouble temp = MSS.getValues();
  int ntemp = MSS.getNRows();
  MatrixSymmetric* MEig   = MatrixSymmetric::createFromVD(temp);
  MEig->display();
  MatrixSymmetric* MNoEig = MatrixSymmetric::createFromVD(temp);
  MNoEig->display();

  // Extract the Eigen values and vectors (both matrix types)
  (void) MEig->computeEigen();
  VectorDouble eigVal = MEig->getEigenValues();
  VH::dump("Eigen Values (Eigen Library)", eigVal);
  const MatrixSquare* eigVec = MEig->getEigenVectors();
  eigVec->display();

  (void) MNoEig->computeEigen();
  VectorDouble eigNoVal = MNoEig->getEigenValues();
  VH::dump("Eigen Values (no Eigen Library)", eigNoVal);
  const MatrixSquare* eigNoVec = MNoEig->getEigenVectors();
  eigNoVec->display();

  // *********************
  // Cholesky calculations
  // *********************

  // Compute Cholesky factorization (only for comparison)
  CholeskyDense MEigChol(MEig);

  // Compare Cholesky Decomposition calculated using Eigen Library or not (sparse matrix only)
  mestitle(0,"Cholesky Decomposition for Sparse matrices");
  reset_to_initial_contents(M, MRR, MSG, MSS, MSP);

  // Get a Sparse matrix
  NF_Triplet triplet = MSP->getMatrixToTriplet();
  MatrixSparse* MSEig = MatrixSparse::createFromTriplet(triplet, MSP->getNRows(), MSP->getNCols(),1);
  MSEig->display();
  MatrixSparse* MSNoEig = MatrixSparse::createFromTriplet(triplet, MSP->getNRows(), MSP->getNCols(),0);
  MSNoEig->display();

  VectorDouble B = VH::simulateGaussian(ntemp);
  VH::dump("Input vector",B);

  VectorDouble XEig(ntemp);
  VectorDouble XNoEig(ntemp);

  CholeskySparse MSEigChol(MSEig);
  MSEigChol.solve(B, XEig);
  VH::dump("Cholesky Solve (Eigen Library)",XEig);
  VectorDouble resEig = MSEig->prodVecMat(XEig);
  VH::dump("Verification (Eigen Library)", resEig);
  MSEigChol.addSimulateToDest(B, XEig);
  // Simulation using Cholesky cannot be compared due to different choices in embedded permutations
  //  VH::dump("Cholesky Simulate (Eigen Library)", XEig);

  CholeskySparse MSNoEigChol(MSNoEig);
  MSNoEigChol.solve(B, XNoEig);
  VH::dump("Cholesky Solve (No Eigen Library)",XNoEig);
  VectorDouble resNoEig = MSNoEig->prodVecMat(XNoEig);
  VH::dump("Verification (no Eigen Library)",resNoEig);
  MSNoEigChol.addSimulateToDest(B, XNoEig);

  // Log Determinant
  mestitle(0,"Cholesky Log Determinant");
  message("Log Determinant Sparse (No Eigen Library)   = %lf\n",
          MSNoEigChol.computeLogDeterminant());
  message("Log Determinant Dense  (traditional method) = %lf\n",
          log(MEig->determinant()));
  message("Log Determinant Dense  (Eigen Library)      = %lf\n",
          MEigChol.computeLogDeterminant());

  // Compute Cholesky factorization (for dense matrix (Eigen library)
  mestitle(0,"Cholesky Decomposition for Dense matrices");
  message("Input Square Symmetric Matrix (Dense format)\n");
  MEig->display();

  // Solving a Linear system after Cholesky decomposition
  mestitle(0,"Solving a Linear system after Cholesky decomposition");
  VH::dump("Input Vector B =", B);
  MEigChol.solve(B, XEig);
  VH::dump("Result Vector X =", XEig);
  message("Is M * X = B: %d\n",(int) VH::isEqual(B,MEig->prodMatVec(XEig)));

  // Solving a linear system after Cholesky decomposition (matrix RHS)
  mestitle(0,"Solving a Linear system after Cholesky decomposition (matrix RHS)");
  int nrows = MEig->getNRows();
  int ncols = 5;
  MatrixDense Bmat(nrows, ncols);
  MatrixDense Bres(nrows, ncols);
  for (int icol = 0; icol < ncols; icol++)
    Bmat.setColumn(icol, B);
  message("Input Matrix B =\n");
  Bmat.display();
  (void) MEigChol.solveMatrix(Bmat, Bres);
  message("Result Matrix X =\n");
  Bres.display();

  MatrixDense* Bcheck = MatrixFactory::prodMatMat<MatrixDense>(MEig, &Bres);
  message("Is M * X = B: %d\n", (int) Bmat.isSame(*Bcheck));
  delete Bcheck;

  // Product by Diagonal built from a vector

  MatrixSparse* MSNDEig = prodNormDiagVec(MSEig, B, 1);
  message("Product by Diagonal from Vector (Eigen)\n");
  MSNDEig->display();
  delete MSNDEig;
  message("Product by Diagonal from Vector In Place with operation (Eigen)\n");
  MSEig->prodNormDiagVecInPlace(B, 2);
  MSEig->display();

  MatrixSparse* MSNDNoEig = prodNormDiagVec(MSNoEig, B, 1);
  message("Product by Diagonal from Vector (No Eigen)\n");
  MSNDNoEig->display();
  delete MSNDNoEig;
  message("Product by Diagonal from Vector In Place with operation (No Eigen)\n");
  MSNoEig->prodNormDiagVecInPlace(B, 2);
  MSNoEig->display();

  // Gluing two sparse matrices

  MatrixSparse* MSGlueEig = dynamic_cast<MatrixSparse*>
    (MatrixFactory::createGlue(MSEig, MSEig, true, true));
  MSGlueEig->display();
  MatrixSparse* MSGlueNoEig = dynamic_cast<MatrixSparse*>
    (MatrixFactory::createGlue(MSNoEig, MSNoEig, true, true));
  MSGlueNoEig->display();

  // Compare Generalized Eigen values calculated using Eigen Library or not (dense matrix only)

  mestitle(0,"Generalized Eigen values calculation for Dense matrices");

  // We use the Square symmetric matrices created in previous paragraph
  // We must construct another square symmetric matrix (B)
  VectorDouble vbh = VH::simulateGaussian(nrow * ncol);

  MatrixDense* MREig = MatrixDense::createFromVD(vbh, nrow, ncol, false);
  AMatrix* MREigt = MREig->transpose();
  MatrixSymmetric* BEig = MatrixFactory::prodMatMat<MatrixSymmetric>(MREig, MREigt);
  delete MREig;
  delete MREigt;

  message("Verify the input for Generalized Eigen calculation\n");
  MEig->display();
  BEig->display();

  // Extract the Generalized Eigen values and vectors (both matrix types)
  (void) MEig->computeGeneralizedEigen(*BEig);
  VectorDouble genEigVal = MEig->getEigenValues();
  const MatrixSquare* genEigVec = MEig->getEigenVectors();
  VH::dump("Generalized Eigen Values (Eigen Library)", genEigVal);
  genEigVec->display();
  delete BEig;

  MatrixDense* MRNoEig = MatrixDense::createFromVD(vbh, nrow, ncol, false);
  AMatrix* MRNoEigt = MRNoEig->transpose();
  MatrixSymmetric* BNoEig = MatrixFactory::prodMatMat<MatrixSymmetric>(MRNoEig, MRNoEigt);
  delete MRNoEig;
  delete MRNoEigt;

  (void) MNoEig->computeGeneralizedEigen(*BNoEig);
  VectorDouble genEigNoVal = MNoEig->getEigenValues();
  const MatrixSquare* genEigNoVec = MNoEig->getEigenVectors();
  VH::dump("Generalized Eigen Values (no Eigen Library)", genEigNoVal);
  genEigNoVec->display();
  delete BNoEig;

  // Free the pointers

  delete M;
  return(0);
}
