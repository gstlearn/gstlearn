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
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"
#include "Matrix/AMatrix.hpp"
#include "Matrix/LinkMatrixSparse.hpp"

#include <Eigen/Sparse>

/**
 * Sparse Matrix
 */
class GSTLEARN_EXPORT MatrixSparse : public AMatrix {

public:
  MatrixSparse(int nrow = 0, int ncol = 0, int opt_eigen=-1);
#ifndef SWIG
  MatrixSparse(const cs* A, int opt_eigen = -1);
#endif
  MatrixSparse(const MatrixSparse &m);
  MatrixSparse& operator= (const MatrixSparse &m);
  virtual ~MatrixSparse();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// Cloneable interface
  IMPLEMENT_CLONING(MatrixSparse)

  /// Interface for AMatrix
  /*! Returns if the current matrix is Sparse */
  bool isSparse() const { return true; }

  /*! Set the contents of a Column */
  virtual void setColumn(int icol, const VectorDouble& tab) override;
  /*! Set the contents of a Row */
  virtual void setRow(int irow, const VectorDouble& tab) override;
  /*! Set the contents of the (main) Diagonal */
  virtual void setDiagonal(const VectorDouble& tab) override;
  /*! Set the contents of the (main) Diagonal to a constant value */
  virtual void setDiagonalToConstant(double value = 1.) override;
  /*! Transpose the matrix and return it as a copy*/
  virtual MatrixSparse* transpose() const override;
  /*! Add a value to each matrix component */
  virtual void addScalar(double v) override;
  /*! Add value to matrix diagonal */
  virtual void addScalarDiag(double v) override;
  /*! Multiply each matrix component by a value */
  virtual void prodScalar(double v) override;
  /*! Set a set of values simultaneously from an input array */
  void setValuesFromTriplet(const Triplet& T) override;
  /*! Set all the values of the matrix at once */
  virtual void fill(double value) override;
  /*! Multiply the matrix row-wise */
  virtual void multiplyRow(const VectorDouble& vec) override;
  /*! Multiply the matrix column-wise */
  virtual void multiplyColumn(const VectorDouble& vec) override;
  /*! Divide the matrix row-wise */
  virtual void divideRow(const VectorDouble& vec) override;
  /*! Divide the matrix column-wise */
  virtual void divideColumn(const VectorDouble& vec) override;
  /*! Perform y = x %*% 'this' */
  virtual VectorDouble prodVecMatInPlace(const VectorDouble& x, bool transpose = false) const override;
  /*! Perform y = 'this' %*% x */
  virtual VectorDouble prodMatVecInPlace(const VectorDouble& x, bool transpose = true) const override;

#ifndef SWIG
  /*! Extract the contents of the matrix */
  virtual Triplet getValuesAsTriplets() const override;
#endif

  //// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static MatrixSparse* createFromTriplet(const Triplet &T,
                                         int nrow,
                                         int ncol,
                                         int opt_eigen = -1);

  void init(int nrows, int ncols);

  /// The next functions use specific definition of matrix (to avoid dynamic_cast)
  /// rather than manipulating AMatrix. They are no more generic of AMatrix
  /*! Add a matrix (multiplied by a constant) */
  virtual void addMatrix(const MatrixSparse& y, double value = 1.);
  /*! Multiply a matrix by another and store the result in the current matrix */
  virtual void prodMatMat(const MatrixSparse& x, const MatrixSparse& y);
  /*! Linear combination of matrices */
  virtual void linearCombination(double cx, double cy, const MatrixSparse& y);

#ifndef SWIG
  /*! Returns a pointer to the Sparse storage */
  const cs* getCS() const
  {
    return _csMatrix;
  }
  void setCS(cs* cs)
  {
    _csMatrix = cs_duplicate(cs);
  }
  void freeCS()
  {
    _csMatrix = cs_spfree2(_csMatrix);
  }
  cs* getCSUnprotected() const { return _csMatrix; } // Temporary function to get the CS contents of Sparse Matrix
#endif
  Triplet getSparseToTriplet(bool flag_from_1 = false) const;

  void reset(int nrows, int ncols);
  void reset(int nrows, int ncols, double value);
  void reset(int nrows, int ncols, const double* tab, bool byCol = true);
  void reset(int nrows,
             int ncols,
             const VectorDouble &tab,
             bool byCol = true);
  void reset(const VectorVectorDouble& tab, bool byCol = true);


  /*! Dump a specific range of samples from the internal storage */
  void dumpElements(const String& title, int ifrom, int ito) const;

  /*! Set all the values of the Matrix with random values */
  void fillRandom(int seed = 432432, double zeroPercent = 0.1);

  int computeCholesky();
  int solveCholesky(const VectorDouble& b, VectorDouble& x);

protected:
  /// Interface for AMatrix
  bool    _isPhysicallyPresent(int irow, int icol) const { DECLARE_UNUSED(irow, icol); return true; }
  bool    _isCompatible(const AMatrix& m) const override { DECLARE_UNUSED(m); return (isSparse()); }
  void    _allocate() override;
  void    _deallocate() override;

  virtual double& _getValueRef(int irow, int icol) override;
  virtual int     _getMatrixPhysicalSize() const override;
  virtual void    _setValueByRank(int rank, double value) override;
  virtual void    _setValue(int irow, int icol, double value) override;
  virtual void    _setValues(const double* values, bool byCol) override;
  virtual double  _getValueByRank(int rank) const override;
  virtual double  _getValue(int irow, int icol) const override;
  virtual int     _getIndexToRank(int irow,int icol) const override;
  virtual void    _transposeInPlace() override;

  virtual void    _prodMatVec(const double *x, double *y, bool transpose = false) const override;
  virtual void    _prodVecMat(const double *x,double *y, bool transpose = false) const override;
  virtual int     _invert() override;
  virtual int     _solve(const VectorDouble& b, VectorDouble& x) const override;

  void _clear();
  bool _isElementPresent(int irow, int icol) const;

#ifndef SWIG
  String _toMatrixLocal(const String& title, bool flagOverride);
#endif

private:
  void _forbiddenForSparse(const String& func) const;

private:
  cs*  _csMatrix; // Classical storage for Sparse matrix
  Eigen::SparseMatrix<double> _eigenMatrix; // Eigen storage in Eigen Library (always stored Eigen::ColMajor)
  bool _flagDecomposeCholesky;

  css *_S; // Cholesky decomposition
  csn *_N; // Cholesky decomposition
  Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > _cholEigen;
};

/*! Transform any matrix into a Sparse format */
GSTLEARN_EXPORT MatrixSparse *createFromAnyMatrix(const AMatrix* mat);
GSTLEARN_EXPORT void setUpdateNonZeroValue(int status = 2);
GSTLEARN_EXPORT int getUpdateNonZeroValue();

// The following functions are added while converting cs into MatrixSparse
#ifndef SWIG
GSTLEARN_EXPORT const cs* _getCS(const MatrixSparse* A, bool optional=false);
GSTLEARN_EXPORT cs* _getCSUnprotected(const MatrixSparse* A, bool optional=false);
#endif

GSTLEARN_EXPORT MatrixSparse* matCS_glue(const MatrixSparse *A1,
                                         const MatrixSparse *A2,
                                         bool shiftRow,
                                         bool shiftCol);
GSTLEARN_EXPORT MatrixSparse* matCS_matvecnorm(const MatrixSparse *A,
                                               const double *x,
                                               int oper);
GSTLEARN_EXPORT MatrixSparse* matCS_prod_norm(int mode,
                                              const MatrixSparse *A,
                                              const MatrixSparse *B);
GSTLEARN_EXPORT MatrixSparse* matCS_eye_tab(int number, double *values);
GSTLEARN_EXPORT MatrixSparse* matCS_eye(int number, double value);
GSTLEARN_EXPORT MatrixSparse* matCS_triplet(const cs *T);
GSTLEARN_EXPORT MatrixSparse* matCS_prod_norm_diagonal(int mode, const MatrixSparse *B, VectorDouble diag);
GSTLEARN_EXPORT MatrixSparse* matCS_transpose(const MatrixSparse *A, int values);
GSTLEARN_EXPORT MatrixSparse* matCS_multiply(const MatrixSparse *A, const MatrixSparse *B);
GSTLEARN_EXPORT MatrixSparse* matCS_add(const MatrixSparse *A, const MatrixSparse *B, double alpha, double beta);
GSTLEARN_EXPORT int           matCS_coarsening(MatrixSparse *Q, int type, int **indCo_ret, MatrixSparse **L_ret);
GSTLEARN_EXPORT MatrixSparse* matCS_interpolate(MatrixSparse *AA, MatrixSparse *Lt, int *Co);
GSTLEARN_EXPORT MatrixSparse* matCS_extract_diag(MatrixSparse *C, int mode);
GSTLEARN_EXPORT MatrixSparse* matCS_extract_submatrix_by_ranks(MatrixSparse *C, int *rank_rows, int *rank_cols);
GSTLEARN_EXPORT int           matCS_gaxpy(const MatrixSparse *A, const double *x, double *y);
GSTLEARN_EXPORT void          matCS_matvecnorm_inplace(MatrixSparse *A, const double* x, int oper);
GSTLEARN_EXPORT double        matCS_norm(const MatrixSparse *A);
GSTLEARN_EXPORT MatrixSparse* matCS_prod_norm_single(int mode, MatrixSparse *B);
GSTLEARN_EXPORT void          matCS_add_value(const MatrixSparse *A, int row, int col, double value);
GSTLEARN_EXPORT void          matCS_set_cste(MatrixSparse *A, double value);
GSTLEARN_EXPORT MatrixSparse* matCS_diag(VectorDouble diag, double tol = EPSILON10);
GSTLEARN_EXPORT double        matCS_get_value(const MatrixSparse *A, int row, int col);
GSTLEARN_EXPORT VectorInt     matCS_color_coding(MatrixSparse *Q, int start, int *ncols);
GSTLEARN_EXPORT MatrixSparse* matCS_extract_submatrix_by_color(MatrixSparse *C,
                                                               const VectorInt &colors,
                                                               int ref_color,
                                                               int row_ok,
                                                               int col_ok);
GSTLEARN_EXPORT void          matCS_rowcol(const MatrixSparse *A,
                                           int *nrows,
                                           int *ncols,
                                           int *count,
                                           double *percent);

GSTLEARN_EXPORT VectorDouble  matCSD_extract_diag_VD(MatrixSparse *C, int mode);
GSTLEARN_EXPORT double*       matCSD_extract_diag(const MatrixSparse *C, int mode);
GSTLEARN_EXPORT int           matCS_scale(MatrixSparse *A);
GSTLEARN_EXPORT void          matCS_print_nice(const char *title, const MatrixSparse *A,
                                               int maxrow, int maxcol);

GSTLEARN_EXPORT MatrixSparse* matCS_normalize_by_diag_and_release(MatrixSparse *Q, int flag_release);
GSTLEARN_EXPORT MatrixSparse* matCS_add_and_release(MatrixSparse *b1, MatrixSparse *b2,
                                                    double alpha, double beta, int flag_release);
GSTLEARN_EXPORT MatrixSparse* matCS_multiply_and_release(MatrixSparse *b1, const MatrixSparse *b2,int flag_release);
