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

/// TODO : Transform into template for storing something else than double

/**
 * Matrix
 */
class GSTLEARN_EXPORT AMatrix : public AStringable, public ICloneable
{
public:
  AMatrix(int nrow = 0, int ncol = 0);
  AMatrix(const AMatrix &m);
  AMatrix& operator= (const AMatrix &m);
  virtual ~AMatrix();

  void init(int nrows, int ncols);
  void reset(int nrows, int ncols, double value = 0.);
  void resetFromArray(int nrows, int ncols, const double* tab, bool byCol = true);
  void resetFromVD(int nrows, int ncols, const VectorDouble &tab, bool byCol = true);
  void resetFromVVD(const VectorVectorDouble& tab, bool byCol = true);

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /*! Returns if the current matrix is Sparse */
  virtual bool isSparse() const { return false; }
  /*! Check if the matrix is (non empty) square */
  virtual bool isSquare(bool printWhyNot = false) const;
  /*! Indicate if the given indices are valid for the current matrix size */
  virtual bool isValid(int irow, int icol, bool printWhyNot = false) const;
  /*! Check if the matrix is square and Identity */
  virtual bool isIdentity(bool printWhyNot = false) const;
  /*! Check if the input matrix is (non empty and square) symmetric */
  virtual bool isSymmetric(bool printWhyNot = false) const;
  /*! Check if the matrix is (non empty) diagonal */
  virtual bool isDiagonal(bool printWhyNot = false) const;
  /*! Check if the contents of the matrix is constant and diagonal */
  virtual bool isDiagCst(bool printWhyNot = false) const;
  /*! Say if the matrix must be symmetric */
  virtual bool mustBeSymmetric() const { return false; }
  /*! Say if the matrix must be diagonal */
  virtual bool mustBeDiagonal() const { return false; }
  /*! Say if the matrix must be diagonal constant */
  virtual bool mustBeDiagCst() const { return false; }

  /*! Set the contents of a Column */
  virtual void setColumn(int icol, const VectorDouble& tab);
  /*! Set the contents of a Row */
  virtual void setRow(int irow, const VectorDouble& tab);
  /*! Set the contents of the (main) Diagonal */
  virtual void setDiagonal(const VectorDouble& tab);
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
  /*! Set a set of values simultaneously from an input array */
  virtual void setValuesByArrays(const VectorInt &irows,
                                 const VectorInt &icols,
                                 const VectorDouble &values);
  /*! Add a matrix (multiplied by a constant) */
  virtual void addMatrix(const AMatrix& y, double value = 1.);
  /*! Multiply a matrix by another and store the result in the current matrix */
  virtual void prodMatrix(const AMatrix& x, const AMatrix& y);
  /*! Linear combination of matrices */
  virtual void linearCombination(double cx, double cy, const AMatrix& y);
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
  /*! Perform M * 'vec' */
  virtual VectorDouble prodVector(const VectorDouble& vec) const;
  /*! Perform 'vec'^T * M */
  virtual VectorDouble prodTVector(const VectorDouble& vec) const;

#ifndef SWIG
  /*! Extract the contents of the matrix */
  virtual void getValuesAsTriplets(VectorInt &irows,
                                   VectorInt &icols,
                                   VectorDouble &values) const;
#endif

  /*! Gets the value at row 'irow' and column 'icol' */
  double getValue(int irow, int icol) const;
  /*! Sets the value at row 'irow' and column 'icol' */
  void setValue(int irow, int icol, double value);
  /*! Add a value to a matrix term */
  void addValue(int irow, int icol, double value);
  /*! Check if a matrix is the same as me (norm L1) */
  bool isSame(const AMatrix& m, double eps = EPSILON10);
  /*! Check that both matrix have the same number of rows and columns */
  bool isSameSize(const AMatrix& m) const;
  /*! Returns if the current matrix is Empty */
  bool isEmpty() const { return (_nRows == 0 || _nCols == 0); }
  /*! Returns the sum of absolute difference between argument and this */
  double compare(const AMatrix& mat) const;
  /*! Returns the number of rows */
  int getNRows() const { return _nRows; }
  /*! Returns the number of columns */
  int getNCols() const { return _nCols; }
  /*! Get the total number of elements of the (full) matrix */
  int getNTotal() const { return _nRows * _nCols; }

  /*! Returns the contents of the whole matrix as a VectorDouble */
  VectorDouble getValues(bool byCol = true) const;
  /*! Extract a Diagonal (main or secondary) of this */
  VectorDouble getDiagonal(int shift = 0) const;
  /*! Extract a Row */
  VectorDouble getRow(int irow) const;
  /*! Extract a Column */
  VectorDouble getColumn(int icol) const;
  /*! Checks if a Column is valid (contains a non TEST value) */
  bool isColumnDefined(int icol) const;
  /*! Checks if a Row is valid (contains a non TEST value) */
  bool isRowDefined(int irow) const;
  /*! Define the number of defined columns */
  int getNumberColumnDefined() const;
  /*! Define the number of defined rows */
  int getNumberRowDefined() const;

  /*! Perform 'outv' = M * 'inv' */
  void prodVectorInPlace(const VectorDouble& inv, VectorDouble& outv) const;
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
  void fillRandom(int seed = 432432, double zeroPercent = 0.1);
  void setValues(const VectorDouble& values, bool byCol=true);
  double getMeanByColumn(int icol) const;
  double getMinimum() const;
  double getMaximum() const;
  void copyReduce(const AMatrix *x,
                  const VectorInt &activeRows,
                  const VectorInt &activeCols);
  void setFlagCheckAddress(bool flagCheckAddress) { _flagCheckAddress = flagCheckAddress; }

#ifndef SWIG
  /*! Get value operator override */
  double  operator()(int row, int col) const { return getValue(row, col); }
  /*! Set value operator override */
  double &operator()(int row, int col)       { return _getValueRef(row, col); }

protected:
  /*! Say if (irow, icol) is stored physically or not */
  virtual bool    _isPhysicallyPresent(int /*irow*/, int /*icol*/) const { return true; }
  virtual bool    _isCompatible(const AMatrix& m) const { return isSameSize(m); }
  virtual double& _getValueRef(int irow, int icol);
  virtual int     _getMatrixPhysicalSize() const;
  virtual void    _setValues(const double* values, bool byCol);
  virtual void    _clearDecoration() {};

  virtual void    _allocate() = 0;
  virtual void    _deallocate() = 0;
  virtual void    _setValueByRank(int rank, double value) = 0;
  virtual double  _getValue(int irow, int icol) const = 0;
  virtual double  _getValueByRank(int rank) const = 0;
  virtual void    _setValue(int irow, int icol, double value) = 0;
  virtual int     _getIndexToRank(int irow,int icol) const = 0;

  virtual void    _transposeInPlace() = 0;
  virtual void    _prodVectorInPlace(const double *inv,double *outv) const = 0;
  virtual int     _invert() = 0;
  virtual int     _solve(const VectorDouble& b, VectorDouble& x) const = 0;

  void _setNCols(int ncols) { _nCols = ncols; }
  void _setNRows(int nrows) { _nRows = nrows; }
  bool _isNumbersValid(int nrows,int ncols) const;
  bool _isColumnValid(int icol) const;
  bool _isRowValid(int irow) const;
  bool _isIndexValid(int irow, int icol) const;
  bool _isRowVectorConsistent(const VectorDouble& tab);
  bool _isColVectorConsistent(const VectorDouble& tab);
  bool _isVectorSizeConsistent(int nrows, int ncols, const VectorDouble& tab);
  bool _isRankValid(int rank) const;
  void _clear();
  void _fillFromVVD(const VectorVectorDouble& X);
  void _recopy(const AMatrix &m);

  bool _getFlagCheckAddress() const { return _flagCheckAddress; }

private:
  int  _nRows;
  int  _nCols;
  bool _flagCheckAddress;
  double _nullTerm; // Used for returning a null constant address
#endif
};

/* Shortcut functions for C style aficionados */
GSTLEARN_EXPORT AMatrix* prodMatrix(const AMatrix* mat1, const AMatrix* mat2);
GSTLEARN_EXPORT void prodMatrixInPlace(AMatrix* mat1, const AMatrix* mat2);
GSTLEARN_EXPORT void setFlagEigen(bool flagEigen);
GSTLEARN_EXPORT bool isFlagEigen();
GSTLEARN_EXPORT void setMultiThread(int nthreads);
GSTLEARN_EXPORT int  getMultiThread();
GSTLEARN_EXPORT bool isMultiThread();
