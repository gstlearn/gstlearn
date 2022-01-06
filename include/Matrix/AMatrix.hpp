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
#pragma once

#include "gstlearn_export.hpp"
#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/IClonable.hpp"
#include "csparse_d.h"

typedef struct {
  VectorInt rows;
  VectorInt cols;
  VectorDouble values;
} cs_Output;

/// TODO : Transform into template for storing something else from double

/**
 * Matrix
 */
class GSTLEARN_EXPORT AMatrix : public AStringable, public IClonable
{
protected:
  AMatrix(int nrow = 0, int ncol = 0, bool sparse = false);
  AMatrix(const cs* A);
  AMatrix(const AMatrix &m);
  AMatrix& operator= (const AMatrix &m);

public:
  virtual ~AMatrix();

  void init(int nrows, int ncols, bool sparse = false);

  /*! Returns the number of rows */
  int  getNRows() const { return _nRows; }
  /*! Returns the number of columns */
  int  getNCols() const { return _nCols; }
  /*! Get the total number of elements of the (full) matrix */
  int  getNTotal() const { return _nRows * _nCols; }
  /*! Gets the value at row 'irow' and column 'icol' */
  virtual double getValue(int irow, int icol) const;
  /*! Gets the value at rank 'rank' */
  virtual double getValue(int rank) const;
  /*! Gets a reference to the value at row 'irow' and column 'icol' */
  virtual double& getValueRef(int irow, int icol);
  /*! Returns the contents of the whole matrix as a VectorDouble */
  VectorDouble getValues() const;
#ifndef SWIG
  /*! Extract the contents of the matrix */
  void getValuesAsTriplets(VectorInt& irows,
                           VectorInt& icols,
                           VectorDouble& values) const;
#endif
  cs_Output getValuesAsTriplets(bool flag_from_1 = true) const;
  /*! Extract a Diagonal (main or secondary) of this */
  VectorDouble getDiagonal(int shift=0) const;
  /*! Extract a Row */
  VectorDouble getRow(int irow) const;
  /*! Extract a Column */
  VectorDouble getColumn(int icol) const;
  /*! Add a matrix to this */
  void add(const AMatrix& tab, double value = 1.);
  /*! Subtract a matrix to this */
  void subtract(const AMatrix& tab, double value = 1.);
  /*! Set the contents of a Column */
  virtual void setColumn(int icol, const VectorDouble& tab);
  /*! Set the contents of a Row */
  virtual void setRow(int irow, const VectorDouble& tab);
  /*! Set the contents of the (main) Diagonal */
  virtual void setDiagonal(const VectorDouble& tab);
  virtual void setDiagonal(double value = 1.);

#ifndef SWIG
  /*! Returns a pointer to the Sparse storage */
  const cs* getCs() const { return _csMatrix; }
#endif

  /*! Check that both matrix have the same number of rows and columns */
  bool isSameSize(const AMatrix& m) const;
  /*! Returns if the current matrix is Sparse */
  bool isSparse() const { return _sparse; }
  /*! Returns if the current matrix is Empty */
  bool isEmpty() const { return (_nRows == 0 || _nCols == 0); }
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
  /*! Check if a matrix is the same as me (norm L1) */
  bool isSame(const AMatrix& m, double eps = EPSILON10);

  void reset(int nrows, int ncols, bool sparse = false);
  void reset(int nrows, int ncols, double value, bool sparse = false);
  void reset(int nrows, int ncols, const double* tab, bool sparse = false);
  void reset(int nrows, int ncols, const VectorDouble& tab, bool sparse = false);

  /*! Returns the sum of absolute difference between argument and this */
  double compare(const AMatrix& mat) const;
  /*! Transform the current matrix (any format) in a Sparse format in place */
  void toSparseInPlace();
  /*! Transform any matrix in a Sparse format */
  AMatrix* toSparse() const;
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
  /*! Product of the Matrix by a vector (on its right) */
  void prodVector(const double *in,double *out) const;
  void prodVector(const VectorDouble& in, VectorDouble& out) const;
  /*! Multiply a Matrix row-wise */
  void multiplyRow(const VectorDouble& vec);
  /*! Multiply a Matrix column-wise */
  void multiplyColumn(const VectorDouble& vec);
  /*! Divide a Matrix row-wise */
  void divideRow(const VectorDouble& vec);
  /*! Divide a Matrix column-wise */
  void divideColumn(const VectorDouble& vec);
  /*! Add a matrix to this component by component */
  virtual void addMatrix(const AMatrix& y);
  /*! Multiply a matrix by another and store the result in the current matrix */
  virtual void prodMatrix(const AMatrix& x, const AMatrix& y);
  /*! Linear combination of matrices */
  virtual void linearCombination(double cx, double cy, const AMatrix& y);
  /*! Matrix inversion in place */
  int invert();
  /*! Solving the Matrix Linear system */
  int solve(const VectorDouble& b, VectorDouble& x) const;
  /*! Calculate the Determinant */
  double determinant() const;

  /*! Dump a specific range of samples from the internal storage */
  void dumpElements(const String& title, int ifrom, int ito) const;
  /*! Conversion to a string */
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /*! Sets the value at row 'irow' and column 'icol' */
  virtual void setValue(int irow, int icol, double value);
  /*! Sets the value at rank 'rank' */
  virtual void setValue(int rank, double value);

  /*! Sets the matrix as Identity */
  void setIdentity(double value = 1.);
  /*! Set all the values of the Matrix at once */
  void fill(double value);
  void fillRandom(int seed = 432432, double zeroPercent = 0.1);
  void setValues(const double* values, bool byCol=true);
  void setValues(const VectorDouble& values, bool byCol=true);
  void setValues(const VectorInt irows,
                 const VectorInt icols,
                 const VectorDouble values);


#ifndef SWIG
  /*! Get value operator override */
  double  operator()(int row, int col) const { return getValue(row, col); }
  /*! Set value operator override */
  double &operator()(int row, int col)       { return getValueRef(row, col); }

protected:
  /*! Say if (irow, icol) is stored physically or not */
  virtual bool    _isPhysicallyPresent(int /*irow*/, int /*icol*/) const { return true; }
  virtual bool    _isCompatible(const AMatrix& m) const = 0;
  virtual void    _allocate() = 0;
  virtual void    _deallocate() = 0;
  virtual int     _getMatrixSize() const = 0;
  virtual void    _setValue(int rank, double value) = 0;
  virtual void    _setValue(int irow, int icol, double value) = 0;
  virtual void    _setValues(const double* values, bool byCol) = 0;
  virtual double  _getValue(int irow, int icol) const = 0;
  virtual double  _getValue(int rank) const = 0;
  virtual double& _getValueRef(int irow, int icol) = 0;
  virtual void    _transposeInPlace() = 0;
  virtual void    _prodVector(const double *in,double *out) const = 0;
  virtual int     _invert() = 0;
  virtual int     _solve(const VectorDouble& b, VectorDouble& x) const = 0;
  virtual double  _determinant() const = 0;

  void _setNCols(int ncols) { _nCols = ncols; }
  void _setNRows(int nrows) { _nRows = nrows; }
  bool _isNumbersValid(int nrows,int ncols) const;
  bool _isIndexValid(int irow, int icol) const;
  bool _isRowVectorConsistent(const VectorDouble& tab);
  bool _isColVectorConsistent(const VectorDouble& tab);
  bool _isVectorSizeConsistent(int nrows, int ncols, const VectorDouble& tab);
  bool _isRankValid(int rank) const;
  void _clear();

private:
  void _setSparse(bool sparse) { _sparse = sparse; }
  void _initiateSparse();
  void _recopySparse(const cs* cs);
  void _deallocateSparse();
  void _forbiddenForSparse(const String& func) const;

private:
  int  _nRows;
  int  _nCols;
  bool _sparse;
  cs*  _csMatrix;
#endif
};

/* Shortcut functions for C style aficionados */
GSTLEARN_EXPORT AMatrix* createIdentity(int nrow, bool sparse);
GSTLEARN_EXPORT AMatrix* transpose(const AMatrix* mat);
GSTLEARN_EXPORT AMatrix* prodMatrix(const AMatrix* mat1, const AMatrix* mat2);
