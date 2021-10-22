/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* Created on: 9 avr. 2019 by N. Desassis                                     */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "Basic/Vector.hpp"
#include "Matrix/AMatrixSquare.hpp"

/**
 * Square Diagonal matrices
 */
class MatrixSDiag : public AMatrixSquare {

public:
  MatrixSDiag(int nrow = 0, bool sparse = false);
  MatrixSDiag(const MatrixSDiag &m);
  MatrixSDiag& operator= (const MatrixSDiag &r);
	virtual ~MatrixSDiag();

  /*! Clonable interface */
  virtual IClonable* clone() const override;

  virtual String toString(int level = 0) const override;

  /*! Transpose the matrix */
  void transposeInPlace() override;

  /*! Indicate if the given indices are valid for the current matrix size */
  bool isValid(int irow, int icol, bool printWhyNot = false) const override;
  /*! does the matrix is symmetrical ? */
  bool isSymmetric(bool printWhyNot = false) const override { return true; }
  /*! Check if the (non empty) matrix is diagonal */
  bool isDiagonal(bool printWhyNot = false) const override { return true; }

  /*! Add a value to each matrix component */
  void addScalar(double v) override;

  /*! Set the contents of a Column */
  void setColumn(int icol, const VectorDouble& tab) override;
  /*! Set the contents of a Row */
  void setRow(int irow, const VectorDouble& tab) override;

protected:
#ifndef SWIG
  virtual double& _getValueRef(int irow, int icol) override;

  /*! Say if the matrix must be symmetric */
  bool mustBeSymmetric() const override { return true; }
  /*! Say if the matrix must be diagonal */
  bool mustBeDiagonal() const override { return true; }
  /*! Say if the matrix must be diagonal constant */
  bool mustBeDiagCst() const override { return false; }

private:
  bool   _isCompatible(const AMatrix& m) const override { return (isSameSize(m) && isDiagonal()); }
  double _getValue(int irow, int icol) const override;
  double _getValue(int irank) const override;
  void   _setValue(int irow, int icol, double value) override;
  void   _setValue(int irank, double value) override;
  void   _transposeInPlace() override { return ; } // Nothing to do
  void   _setValues(const double* values, bool byCol = true) override;
  int    _getMatrixSize() const override;
  void   _allocate() override;
  void   _deallocate() override;
  void   _prodVector(const double *in,double *out) const override;
  int    _invert() override;
  int    _solve(const VectorDouble& b, VectorDouble& x) const override;
  double _determinant() const override;

  void   _recopy(const MatrixSDiag &r);
  bool   _isIndexValid(int irow,int icol) const;

private:
  VectorDouble _diagMatrix;
#endif
};
