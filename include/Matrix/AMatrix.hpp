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
#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"

#ifndef SWIG
#include <Eigen/Core>
#include <Eigen/Dense>
#endif

#include <Eigen/src/Core/Matrix.h>

/// TODO : Transform into template for storing something else than double

class NF_Triplet;
class EOperator;

/**
 * This class is the root of the Matrix organization in gstlearn
 * A matrix is a 2-D organization: it is characterized by its number of rows
 * and its number of columns.
 * Although the user should not bother with this remark, the elements of a matrix
 * processed in 'gstlearn' are stored in a Row-major format.
 * This is to say that the internal rank of an element characterized by its row and column numbers is:
 *  (icol * getNRows() + irow)
 * 
 * Since gstlearn version v1.3:
 * - Dense Matrices storage and algebra rely on Eigen3 library only
 * - Sparse Matrices storage and algebra rely on Eigen3 or cs library (see MatrixSparse.hpp)
 */
class GSTLEARN_EXPORT AMatrix : public AStringable, public ICloneable
{
public:
  AMatrix(int nrow = 0, int ncol = 0);
  AMatrix(const AMatrix &m);
  AMatrix& operator= (const AMatrix &m);
  virtual ~AMatrix();

  virtual void reset(int nrows, int ncols);
  virtual void resetFromValue(int nrows, int ncols, double value);
  virtual void resetFromArray(int nrows, int ncols, const double* tab, bool byCol = true);
  virtual void resetFromVD(int nrows, int ncols, const VectorDouble &tab, bool byCol = true);
  virtual void resetFromVVD(const VectorVectorDouble& tab, bool byCol = true);

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface to AMatrix
  /*! Returns if the matrix belongs to the AMatrixDense class (avoids dynamic_cast) */
  virtual bool isDense() const = 0;
  /*! Returns if the current matrix is Sparse */
  virtual bool isSparse() const = 0;
  /*! Check if the matrix is (non empty) square */
  virtual bool isSquare(bool printWhyNot = false) const;
  /*! Indicate if the given indices are valid for the current matrix size */
  virtual bool isValid(int irow, int icol, bool printWhyNot = false) const;
  /*! Check if the matrix is square and Identity */
  virtual bool isIdentity(bool printWhyNot = false) const;
  /*! Check if the input matrix is (non empty and square) symmetric */
  virtual bool isSymmetric(bool printWhyNot = false, double eps = EPSILON10) const;
  /*! Say if the matrix must be symmetric */
  virtual bool mustBeSymmetric() const { return false; }

  /*! Set the contents of a Column */
  virtual void setColumn(int icol, const VectorDouble& tab, bool flagCheck=true);
  /*! Set the contents of a Row */
  virtual void setRow(int irow, const VectorDouble& tab, bool flagCheck=true);
  /*! Set the contents of the (main) Diagonal */
  virtual void setDiagonal(const VectorDouble& tab, bool flagCheck=true);
  /*! Set the contents of the (main) Diagonal to a constant value */
  virtual void setDiagonalToConstant(double value = 1.);
  /*! Transpose the matrix in place*/
  virtual void transposeInPlace();
  /*! Transpose the matrix and return it as a copy*/
  virtual AMatrix* transpose() const;
  /*! Add a value to each matrix component */
  virtual void addScalar(double v);
  /*! Add value to matrix diagonal */
  virtual void addScalarDiag(double v);
  /*! Multiply each matrix component by a value */
  virtual void prodScalar(double v);
  /*! Set all the values of the Matrix at once */
  virtual void fill(double value);
  /*! Multiply a Matrix row-wise */
  virtual void multiplyRow(const VectorDouble& vec);
  /*! Multiply a Matrix column-wise */
  virtual void multiplyColumn(const VectorDouble& vec);
  /*! Divide a Matrix row-wise */
  virtual void divideRow(const VectorDouble& vec);
  /*! Divide a Matrix column-wise */
  virtual void divideColumn(const VectorDouble& vec);
  /*! Perform 'vec' * 'this' */
  virtual VectorDouble prodVecMat(const VectorDouble& x, bool transpose = false) const;
  /*! Perform 'this' * 'vec' */
  virtual VectorDouble prodMatVec(const VectorDouble& x, bool transpose = false) const;
  /*! Extract a Row */
  virtual VectorDouble getRow(int irow) const;
  /*! Extract a Column */
  virtual VectorDouble getColumn(int icol) const;
  /*! Multiply matrix 'x' by matrix 'y' and store the result in 'this' */
  virtual void prodMatMatInPlace(const AMatrix *x,
                                 const AMatrix *y,
                                 bool transposeX = false,
                                 bool transposeY = false);
  /*! Extract the contents of the matrix */
  virtual NF_Triplet getMatrixToTriplet(int shiftRow=0, int shiftCol=0) const;
  /*! Add a matrix (multiplied by a constant) */
  void addMatInPlace(const AMatrix& y, double cx = 1., double cy = 1.);
  /*! Multiply 'this' by matrix 'y' and store in 'this'*/
  void prodMatInPlace(const AMatrix* matY, bool transposeY = false);
  /*! Product 't(A)' %*% 'M' %*% 'A' or 'A' %*% 'M' %*% 't(A)' stored in 'this'*/
  void prodNormMatMatInPlace(const AMatrix* a,
                             const AMatrix* m,
                             bool transpose = false);
  /*! Product 't(A)' %*% ['vec'] %*% 'A' or 'A' %*% ['vec'] %*% 't(A)' stored in 'this'*/
  void prodNormMatInPlace(const AMatrix &a,
                          const VectorDouble &vec = VectorDouble(),
                          bool transpose = false);

  /*! Modify the dimension of the matrix (if needed) */
  void resize(int nrows, int ncols);
  /*! Gets the value at row 'irow' and column 'icol' */
  virtual double getValue(int irow, int icol, bool flagCheck=true) const = 0;
  /*! Sets the value at row 'irow' and column 'icol' */
  virtual void setValue(int irow, int icol, double value, bool flagCheck=true) = 0;
  /*! Update the value at row 'irow' and column 'icol' */
  virtual void updValue(int irow,
                        int icol,
                        const EOperator &oper,
                        double value,
                        bool flagCheck = true) = 0;
  /*! Add a value to a matrix term */
  void addValue(int irow, int icol, double value);
  /*! Check if a matrix is the same as me (norm L1) */
  bool isSame(const AMatrix& m, double eps = EPSILON4, bool printWhyNot = false);
  /*! Check that both matrix have the same number of rows and columns */
  bool isSameSize(const AMatrix& m) const;
  /*! Returns if the current matrix is Empty */
  bool empty() const { return (_nRows == 0 || _nCols == 0); }
  /*! Returns the sum of absolute difference between argument and this */
  double compare(const AMatrix& mat) const;
  /*! Returns the number of rows */
  int getNRows() const { return _nRows; }
  /*! Returns the number of columns */
  int getNCols() const { return _nCols; }
  /*! Get the total number of elements of the (full) matrix */
  /* The name has been chosen by analogy to VectorT class */
  int size() const { return _nRows * _nCols; }

  /*! Returns the contents of the whole matrix as a VectorDouble */
  VectorDouble getValues(bool byCol = true) const;
  /*! Extract a Diagonal (main or secondary) of this */
  VectorDouble getDiagonal(int shift = 0) const;
  /*! Checks if a Column is valid (contains a non TEST value) */
  bool isColumnDefined(int icol) const;
  /*! Checks if a Row is valid (contains a non TEST value) */
  bool isRowDefined(int irow) const;
  /*! Define the number of defined columns */
  int getNumberColumnDefined() const;
  /*! Define the number of defined rows */
  int getNumberRowDefined() const;
  /*! Check if the matrix does not contain any negative element */
  bool isNonNegative(bool verbose = false) const;

  /*! Perform 'y' = 'this' * 'x' */
  void prodMatVecInPlace(const VectorDouble& x, VectorDouble& y, bool transpose = false) const;
  #ifndef SWIG
    void prodMatVecInPlace(const Eigen::VectorXd& x, Eigen::VectorXd& y, bool transpose = false) const;
  #endif
  void prodMatVecInPlacePtr(const double* x, double* y, bool transpose = false) const;
  /*! Perform 'y' = 'x' * 'this' */
  void prodVecMatInPlace(const VectorDouble& x, VectorDouble& y, bool transpose = false) const;
  void prodVecMatInPlacePtr(const double* x, double* y, bool transpose = false) const;

  /*! Perform x %*% mat %*% y */
  double quadraticMatrix(const VectorDouble& x, const VectorDouble& y);
  /*! Matrix inversion in place */
  int invert();
  /*! Solving the Matrix Linear system */
  int solve(const VectorDouble& b, VectorDouble& x) const;
  /*! Dump a specific range of samples from the internal storage */
  void dumpElements(const String& title, int ifrom, int ito) const;
  /*! Sets the matrix as Identity */
  void setIdentity(double value = 1.);
  void fillRandom(int seed = 432432, double zeroPercent = 0);
  void setValues(const VectorDouble& values, bool byCol=true);
  double getMeanByColumn(int icol) const;
  double getMinimum() const;
  double getMaximum() const;
  double getNormInf() const;
  void copyReduce(const AMatrix *x,
                  const VectorInt &validRows,
                  const VectorInt &validCols);
  void copyElements(const AMatrix &m, double factor = 1.);
  void setFlagCheckAddress(bool flagCheckAddress) { _flagCheckAddress = flagCheckAddress; }

  void makePositiveColumn();
  void linearCombination(double val1,
                         const AMatrix* mat1,
                         double val2         = 1.,
                         const AMatrix* mat2 = nullptr,
                         double val3         = 1.,
                         const AMatrix* mat3 = nullptr);

#ifndef SWIG
  /*! Get value operator override */
  double  operator()(int row, int col) const { return getValue(row, col); }
  /*! Set value operator override */
  double &operator()(int row, int col)       { return _getValueRef(row, col); }

protected:
  virtual void    _allocate() = 0;
  virtual void    _deallocate() = 0;

  /*! Say if (irow, icol) is stored physically or not */
  virtual bool    _isPhysicallyPresent(int /*irow*/, int /*icol*/) const { return true; }
  virtual double& _getValueRef(int irow, int icol);
  virtual int     _getMatrixPhysicalSize() const;
  virtual void    _setValues(const double* values, bool byCol);

  virtual void    _setValueByRank(int rank, double value) = 0;
  virtual double  _getValueByRank(int rank) const = 0;
  virtual int     _getIndexToRank(int irow,int icol) const = 0;

  virtual void    _transposeInPlace() = 0;
  virtual void    _prodMatVecInPlacePtr(const double *x,
                                        double *y,
                                        bool transpose = false) const = 0;
  virtual void    _addProdMatVecInPlaceToDestPtr(const double *x,
                                                 double *y,
                                                 bool transpose = false) const = 0;                                      
  virtual void    _prodVecMatInPlacePtr(const double *x,
                                        double *y,
                                        bool transpose = false) const = 0;
  virtual int     _invert() = 0;
  virtual int     _solve(const VectorDouble& b, VectorDouble& x) const = 0;
  virtual void    _clear();
  virtual bool    _isNumbersValid(int nrows, int ncols) const;

  void _setNCols(int ncols) { _nCols = ncols; }
  void _setNRows(int nrows) { _nRows = nrows; }
  bool _isColumnValid(int icol) const;
  bool _isRowValid(int irow) const;
  bool _isIndexValid(int irow, int icol) const;
  bool _isRowVectorConsistent(const VectorDouble& tab) const;
  bool _isColVectorConsistent(const VectorDouble& tab) const;
  bool _isVectorSizeConsistent(const VectorDouble& tab) const;
  bool _isColumnSizeConsistent(const VectorDouble &tab) const;
  bool _isRowSizeConsistent(const VectorDouble &tab) const;
  bool _isRankValid(int rank) const;
  void _fillFromVVD(const VectorVectorDouble& X);

  bool _getFlagCheckAddress() const { return _flagCheckAddress; }

  bool _checkLink(int nrow1,
                  int ncol1,
                  bool transpose1,
                  int nrow2 = 0,
                  int ncol2 = 0,
                  bool transpose2 = false,
                  int nrow3 = 0,
                  int ncol3 = 0,
                  bool transpose3 = false) const;

private:
  int  _nRows;
  int  _nCols;
  bool _flagCheckAddress;
  double _nullTerm; // Used for returning a null constant address
#endif
};

/* Shortcut functions for C style aficionados */
GSTLEARN_EXPORT void setMultiThread(int nthreads);
GSTLEARN_EXPORT int  getMultiThread();
GSTLEARN_EXPORT bool isMultiThread();
