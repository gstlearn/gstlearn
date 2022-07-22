/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Basic/Vector.hpp"
#include "Matrix/AMatrixSquare.hpp"

/**
 * Square Matrix General
 */
class GSTLEARN_EXPORT MatrixSquareGeneral : public AMatrixSquare {

public:
  MatrixSquareGeneral(int nrow = 0, bool sparse = false);
  MatrixSquareGeneral(const MatrixSquareGeneral &m);
  MatrixSquareGeneral(const AMatrix &m);
  MatrixSquareGeneral& operator= (const MatrixSquareGeneral &r);
	virtual ~MatrixSquareGeneral();

  /*! Clonable interface */
  virtual ICloneable* clone() const override { return new MatrixSquareGeneral(*this); };

  /*! Say if the matrix must be symmetric */
  bool mustBeSymmetric() const override { return false; }
  /*! Say if the matrix must be diagonal */
  bool mustBeDiagonal() const override { return false; }
  /*! Say if the matrix must be diagonal constant */
  bool mustBeDiagCst() const override { return false; }

  /*! Check if the matrix is (non empty) square */
  bool isSquare(bool /*printWhyNot*/ = false) const override { return 1; }

protected:
#ifndef SWIG
  virtual double& _getValueRef(int irow, int icol) override;

private:
  bool   _isCompatible(const AMatrix& m) const override { return (isSameSize(m) && isSquare()); }
  double _getValue(int irow, int icol) const override;
  double _getValue(int irank) const override;
  void   _setValue(int irow, int icol, double value) override;
  void   _setValue(int irank, double value) override;
  void   _transposeInPlace() override;
  void   _setValues(const double *values, bool byCol = true) override;
  int    _getMatrixSize() const override;
  void   _allocate() override;
  void   _deallocate() override;
  void   _prodVector(const double *in,double *out) const override;
  int    _invert() override;
  int    _solve(const VectorDouble& b, VectorDouble& x) const override;

  void   _recopy(const MatrixSquareGeneral &r);
  int    _getIndexToRank(int irow,int icol) const;

private:
  VectorDouble _squareMatrix;
#endif
};
