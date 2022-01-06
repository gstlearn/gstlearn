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

#include "gstlearn_export.hpp"
#include "Basic/Vector.hpp"
#include "Matrix/AMatrixSquare.hpp"

class AMatrix;

/**
 * Square Symmetric matrices are stored as Lower Triangular matrices stored by column
 */
class GSTLEARN_EXPORT MatrixSquareSymmetric : public AMatrixSquare {

public:
  MatrixSquareSymmetric(int nrow = 0, bool sparse = false);
  MatrixSquareSymmetric(const MatrixSquareSymmetric &m);
  MatrixSquareSymmetric(const AMatrix &m);
  MatrixSquareSymmetric& operator= (const MatrixSquareSymmetric &r);
	virtual ~MatrixSquareSymmetric();

  /// Clonable interface
  virtual IClonable* clone() const override { return new MatrixSquareSymmetric(*this); };

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /*! Say if the matrix must be symmetric */
  bool mustBeSymmetric() const override { return true; }
  /*! Say if the matrix must be diagonal */
  bool mustBeDiagonal() const override { return false; }
  /*! Say if the matrix must be diagonal constant */
  bool mustBeDiagCst() const override { return false; }

  /// TODO : isPositiveDefinite

  /// Is the matrix symmetrical ?
  bool isSymmetric(bool /*printWhyNot*/ = false) const override { return true; }

  void initMatTri(int nsize,double* tab);

protected:
#ifndef SWIG
  virtual double& _getValueRef(int irow, int icol) override;

private:
  bool _isCompatible(const AMatrix& m) const override
  {
    return (isSameSize(m) && isSymmetric());
  }
  bool   _isPhysicallyPresent(int irow, int icol) const override;
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

  void   _recopy(const MatrixSquareSymmetric &r);
  int    _getIndexToRank(int irow,int icol) const;

private:
  VectorDouble _squareSymMatrix;
#endif
};
