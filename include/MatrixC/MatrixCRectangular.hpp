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
#include "MatrixC/AMatrixC.hpp"

/**
 * Rectangular matrices are stored by columns
 */
class MatrixCRectangular : public AMatrixC {

public:
  MatrixCRectangular(int nrow = 0, int ncol = 0, bool sparse = false);
  MatrixCRectangular(const MatrixCRectangular &m);
  MatrixCRectangular& operator= (const MatrixCRectangular &r);
	virtual ~MatrixCRectangular();

  /*! Clonable interface */
  virtual IClonable* clone() const override;

  /*! Say if the matrix must be symmetric */
  bool mustBeSymmetric() const override { return false; }
  /*! Say if the matrix must be diagonal */
  bool mustBeDiagonal() const override { return false; }
  /*! Say if the matrix must be diagonal constant */
  bool mustBeDiagCst() const override { return false; }

protected:
#ifndef SWIG
  virtual double& _getValueRef(int irow, int icol) override;

private:
  bool   _isCompatible(const AMatrixC& m) const override { return isSameSize(m); }
  double _getValue(int irow, int icol) const override;
  double _getValue(int irank) const override;
  void   _setValue(int rank, double value) override;
  void   _setValue(int irow, int icol, double value) override;
  void   _setValues(const double* values, bool byCol=true) override;
  int    _getMatrixSize() const override;
  void   _transposeInPlace() override;
  void   _allocate() override;
  void   _deallocate() override;
  void   _prodVector(const double *in,double *out) const override;
  int    _invert() override;
  int    _solve(const VectorDouble& b, VectorDouble& x) const override;
  double _determinant() const override;

  void   _recopy(const MatrixCRectangular &r);
  int    _getIndexToRank(int irow,int icol) const;

private:
  VectorDouble _rectMatrix;
#endif
};
