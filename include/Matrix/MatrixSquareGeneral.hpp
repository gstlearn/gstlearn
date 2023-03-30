/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"
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

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// ICloneable interface
  IMPLEMENT_CLONING(MatrixSquareGeneral)

  /*! Say if the matrix must be symmetric */
  bool mustBeSymmetric() const override { return false; }
  /*! Say if the matrix must be diagonal */
  bool mustBeDiagonal() const override { return false; }
  /*! Say if the matrix must be diagonal constant */
  bool mustBeDiagCst() const override { return false; }

  /*! Check if the matrix is (non empty) square */
  bool isSquare(bool printWhyNot = false) const override { DECLARE_UNUSED(printWhyNot); return 1; }

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
  void   _prodVector(const double *inv,double *outv) const override;
  int    _invert() override;
  int    _solve(const VectorDouble& b, VectorDouble& x) const override;

  void   _recopy(const MatrixSquareGeneral &r);
  int    _getIndexToRank(int irow,int icol) const;

private:
  VectorDouble _squareMatrix;
#endif
};
