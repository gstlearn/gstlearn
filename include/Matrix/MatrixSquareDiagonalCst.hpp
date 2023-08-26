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
#include "Matrix/MatrixSquareDiagonal.hpp"

/**
 * Square Diagonal Matrices with diagonal filled with a constant value
 */
class GSTLEARN_EXPORT MatrixSquareDiagonalCst : public MatrixSquareDiagonal {

public:
  MatrixSquareDiagonalCst(int nrow = 0);
  MatrixSquareDiagonalCst(const MatrixSquareDiagonalCst &m);
  MatrixSquareDiagonalCst& operator= (const MatrixSquareDiagonalCst &r);
	virtual ~MatrixSquareDiagonalCst();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// ICloneable interface
  IMPLEMENT_CLONING(MatrixSquareDiagonalCst)

	/// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static MatrixSquareDiagonalCst* createFromVVD(const VectorVectorDouble& X);

  /*! Say if the matrix must be diagonal constant */
  bool mustBeDiagCst() const override { return true; }

  /*! Transpose the matrix */
  void transposeInPlace() override;

  /*! Add a value to each matrix component */
  void addScalar(double v) override;
  /*! Add value to matrix diagonal */
  void addScalarDiag(double v) override;
  /*! Indicate if the given indices are valid for the current matrix size */
  bool isValid(int irow, int icol, bool printWhyNot = false) const override;

  /*! Set the contents of a Column */
  void setColumn(int icol, const VectorDouble& tab) override;
  /*! Set the contents of a Row */
  void setRow(int irow, const VectorDouble& tab) override;
  /*! Set the contents of the (main) Diagonal */
  void setDiagonal(const VectorDouble& tab) override;

  double determinant() const override;

protected:
#ifndef SWIG
  virtual double& _getValueRef(int irow, int icol) override;

private:
  bool   _isCompatible(const AMatrix& m) const override { return (isSameSize(m) && isDiagCst()); }
  int    _getIndexToRank(int irow,int icol) const override;
  double _getValue(int irow, int icol) const override;
  double _getValue(int irank) const override;
  void   _setValue(int irow, int icol, double value) override;
  void   _setValue(int irank, double value) override;
  void   _setValues(const double* values, bool byCol = true) override;
  int    _getMatrixSize() const override { return 1; }
  void   _allocate() override { return; } // nothing to be done there
  void   _deallocate() override { return; } // nothing to be done there
  void   _prodVector(const double *inv,double *outv) const override;
  int    _invert() override;
  int    _solve(const VectorDouble& b, VectorDouble& x) const override;
  bool   _isValidIndex(int irow, int icol) const;
  bool   _isPhysicallyPresent(int irow, int icol) const override;

private:
  double _cstDiagMatrix;
#endif
};
