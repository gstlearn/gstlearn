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
#include "geoslib_d.h"
#include "geoslib_old_f.h"

#include "Basic/MathFunc.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "Matrix/LinkMatrixSparse.hpp"

/****************************************************************************/
/*!
** Main Program for testing the sparse matrix algebra
**
*****************************************************************************/
int main(int argc, char *argv[])
{
  int nrow, ncol, *rank_rows, *rank_cols;
  double *urow,*ucol,*work,*work2;
  cs *Atriplet,*Wtriplet,*Dtriplet;
  cs *A,*At,*Bl,*Bu,*B,*Diag,*Mwork,*W,*D;
  static int flag_leak = 1;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  // Initializations
  Atriplet = Wtriplet = Dtriplet = nullptr;
  A = At = Bl = Bu = B = W = D = Diag = Mwork = nullptr;
  urow = ucol = work = work2 = nullptr;
  memory_leak_set(flag_leak);
  
  // Dimension of the matrix
  nrow = 5;
  ncol = 6;
  law_set_random_seed(1234);

  // Creating the vectors 'u' and 'v'
  urow = (double *) mem_alloc(sizeof(double) * nrow, 1);
  for (int i=0; i<nrow; i++) urow[i] = law_gaussian();
  ucol = (double *) mem_alloc(sizeof(double) * ncol, 1);
  for (int i=0; i<ncol; i++) ucol[i] = law_gaussian();
  rank_rows = (int *) mem_alloc(sizeof(int) * nrow,1);
  for (int i=0; i<nrow; i++) rank_rows[i] = -1;
  rank_cols = (int *) mem_alloc(sizeof(int) * nrow,1);
  for (int i=0; i<ncol; i++) rank_cols[i] = -1;
  rank_rows[1] = 0;
  rank_rows[3] = 1;
  rank_cols[2] = 0;
  rank_cols[4] = 1;
  rank_cols[5] = 2;

  message("Creating two testing vectors\n");
  print_matrix("URow (5 elements)",0,1,1,nrow,NULL,urow);
  print_matrix("UCol (6 elements)",0,1,ncol,1,NULL,ucol);

  work  = (double *) mem_alloc(sizeof(double) * MAX(ncol,nrow), 1);
  work2 = (double *) mem_alloc(sizeof(double) * MAX(ncol,nrow), 1);

  // Create the non-symmetric sparse matrix
  Atriplet = cs_spalloc2(0, 0, 1, 1, 1);
  for (int irow=0; irow<nrow; irow++)
    for (int icol=0; icol<ncol; icol++)
      cs_entry2(Atriplet, irow, icol, law_gaussian());
  A = cs_triplet2(Atriplet);

  // Print the matrix using Tim Davies printout
  message("Printing A using Tim Davies' primitive\n");
  cs_print2(A, 0);

  // Print A complete
  mestitle(2,"Testing nice printout of A");
  cs_print_nice("cs_print_nice:",A,-1,-1);
  cs_print_dim ("cs_print_dim:",A);
  cs_print_nice("cs_print_nice (limited to nrows=2 and ncols=4):",A,2,4);

  // Transpose the matrix A
  mestitle(2,"Transposing a Matrix: t(A)");
  At = cs_transpose2(A,1);
  cs_print_nice("cs_transpose:",At,-1,-1);
  
  // Filling a Sparse matrix with several values for the same cell
  Dtriplet = cs_spfree2(Dtriplet);
  Dtriplet = cs_spalloc2(0, 0, 1, 1, 1);
  cs_entry2(Dtriplet, 2, 1, 1.);
  cs_entry2(Dtriplet, 1, 3, 2.);
  cs_entry2(Dtriplet, 2, 1, 3.);
  cs_entry2(Dtriplet, 1, 1, 4.);
  cs_entry2(Dtriplet, 1, 3, 5.);
  D = cs_triplet2(Dtriplet);
  cs_print_nice("Filled matrix with several values per cell",D, -1,-1);

  // Create a symmetric square matrix
  mestitle(2,"Creating a Symmetric Matrix: B = A %*% t(A)");
  B = cs_multiply2(A,At);
  cs_print_nice("cs_multiply:",B,-1,-1);

  // Get the Upper triangular part (diagonal included)
  mestitle(2,"Getting the Upper triangular part of a matrix: (U+D)(B)");
  Bu = cs_triangle(B,1,1);
  cs_print_nice("cs_triangle (with diagonal):",Bu,-1,-1);

  // Get the Lower triangular part (diagonal excluded)
  mestitle(2,"Getting the Lower triangular part of a matrix: L(B)");
  Bl = cs_triangle(B,0,0);
  cs_print_nice("cs_triangle (without diagonal):",Bl,-1,-1);

  // Left Product of matrix by vector
  mestitle(2,"Testing product of matrix by vector: A %*% W");
  Wtriplet = cs_spalloc2(0, 0, 1, 1, 1);
  for (int icol=0; icol<ncol; icol++)
    cs_entry2(Wtriplet, icol, 0, ucol[icol]);
  W = cs_triplet2(Wtriplet);
  Wtriplet = cs_spfree2(Wtriplet);
  Mwork = cs_multiply2(A, W);
  cs_print_nice("cs_multiply:",Mwork,-1,-1);
  cs_mulvec(A, nrow, ucol, work);
  print_matrix("Should be equal to:",0,1,nrow,1,NULL,work);
  W = cs_spfree2(W);
  Mwork = cs_spfree2(Mwork);

  // Left Product of matrix by vector
  mestitle(2,"Testing product of matrix by vector: At %*% W");
  Wtriplet = cs_spalloc2(0, 0, 1, 1, 1);
  for (int irow=0; irow<nrow; irow++)
    cs_entry2(Wtriplet, irow, 0, urow[irow]);
  W = cs_triplet2(Wtriplet);
  Wtriplet = cs_spfree2(Wtriplet);
  Mwork = cs_multiply2(At, W);
  cs_print_nice("cs_mulvec:",Mwork,-1,-1);
  cs_mulvec(At, ncol, urow, work);
  print_matrix("Should be equal to:",0,1,ncol,1,NULL,work);
  W = cs_spfree2(W);
  Mwork = cs_spfree2(Mwork);

  // Right Product of matrix by vector
  mestitle(2,"Testing product of matrix (transposed) by vector: W %*% (A)");
  Wtriplet = cs_spalloc2(0, 0, 1, 1, 1);
  for (int irow=0; irow<nrow; irow++)
    cs_entry2(Wtriplet, 0, irow, urow[irow]);
  W = cs_triplet2(Wtriplet);
  Wtriplet = cs_spfree2(Wtriplet);
  Mwork = cs_multiply2(W, A);
  cs_print_nice("cs_tmulvec:",Mwork,-1,-1);
  cs_tmulvec(A, ncol, urow, work);
  print_matrix("Should be equal to:",0,1,ncol,1,NULL,work);
  W = cs_spfree2(W);
  Mwork = cs_spfree2(Mwork);

  // Testing the product of upper triangle by vector
  mestitle(2,"Testing product of upper triangle by vector: (U+D)(B) %*% W");
  cs_mulvec_uptri(B, ncol, urow, work, 1);
  print_matrix("cs_mulvec_uptri:",0,1,1,ncol,NULL,work);
  cs_mulvec(Bu, ncol, urow, work);
  print_matrix("Should be equal to:",0,1,ncol,1,NULL,work);

  // Testing the product of lower triangle by vector
  mestitle(2,"Testing product of lower triangle by vector: L(B) %*% W");
  cs_mulvec_lowtri(B, ncol, urow, work, 0);
  print_matrix("cs_mulvec_lowtri:",0,1,1,ncol,NULL,work);
  cs_mulvec(Bl, ncol, urow, work);
  print_matrix("Should be equal to:",0,1,ncol,1,NULL,work);

  // Solving the Lower Triangular system
  mestitle(2,"Solving the Lower Triangular system: L(B) %*% (X) = W");
  Bl = cs_spfree2(Bl);
  Bl = cs_triangle(B,0,1);
  cs_mulvec(Bl, ncol, urow, work2);
  cs_lsolve_lowtri(B, work2, work);
  print_matrix("cs_lsolve_lowtri:",0,1,1,nrow,NULL,work);
  print_matrix("Should be:",0,1,nrow,1,NULL,urow);

  // Solving the Upper Triangular system
  mestitle(2,"Solving the Upper Triangular system: (U+D)(B) %*% (X) = W");
  Bu = cs_spfree2(Bu);
  Bu = cs_triangle(B,1,1);
  cs_mulvec(Bu, ncol, urow, work2);
  cs_lsolve_uptri(B, work2, work);
  print_matrix("cs_lsolve_uptri:",0,1,1,nrow,NULL,work);
  print_matrix("Should be:",0,1,nrow,1,NULL,urow);

  // Testing extraction
  mestitle(2,"Testing Extraction from A (by ranges)");
  cs_print_nice("Initial Matrix",A,-1,-1);
  Mwork = cs_extract_submatrix(A,1,2,2,3);
  cs_print_nice("Extracted rows[from 2 length 2] and lines[from 3 length 3]",
                Mwork,-1,-1);
  Mwork = cs_spfree2(Mwork);

  // Testing extraction
  mestitle(2,"Testing Extraction from A (by indices)");
  cs_print_nice("Initial Matrix",A,-1,-1);
  Mwork = cs_extract_submatrix_by_ranks(A,rank_rows,rank_cols);
  print_imatrix("Row Selection:",0,1,nrow,1,NULL,rank_rows);
  print_imatrix("Col Selection:",0,1,ncol,1,NULL,rank_cols);
  cs_print_nice("Extracted Matrix:",Mwork,-1,-1);
  Mwork = cs_spfree2(Mwork);

  // Free the structures
  A         = cs_spfree2(A);
  At        = cs_spfree2(At);
  Bl        = cs_spfree2(Bl);
  Bu        = cs_spfree2(Bu);
  B         = cs_spfree2(B);
  W         = cs_spfree2(W);
  D         = cs_spfree2(D);
  Diag      = cs_spfree2(Diag);
  Mwork     = cs_spfree2(Mwork);
  Atriplet  = cs_spfree2(Atriplet);
  Dtriplet  = cs_spfree2(Dtriplet);
  Wtriplet  = cs_spfree2(Wtriplet);
  urow      = (double *) mem_free((char *) urow);
  ucol      = (double *) mem_free((char *) ucol);
  work      = (double *) mem_free((char *) work);
  work2     = (double *) mem_free((char *) work2);
  rank_rows = (int    *) mem_free((char *) rank_rows);
  rank_cols = (int    *) mem_free((char *) rank_cols);
  memory_leak_report();
  return(0);
}
