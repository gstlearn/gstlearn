/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Basic/VectorNumT.hpp"
#include "Matrix/AMatrix.hpp"

/**
 * Rectangular matrices are stored by columns
 */
class GSTLEARN_EXPORT MatrixRectangular : public AMatrix {

public:
  MatrixRectangular(int nrow = 0, int ncol = 0, bool sparse = false);
  MatrixRectangular(const MatrixRectangular &m);
  MatrixRectangular& operator= (const MatrixRectangular &r);
	virtual ~MatrixRectangular();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

	/// Cloneable interface
  IMPLEMENT_CLONING(MatrixRectangular)

  /*! Say if the matrix must be symmetric */
  bool mustBeSymmetric() const override { return false; }
  /*! Say if the matrix must be diagonal */
  bool mustBeDiagonal() const override { return false; }
  /*! Say if the matrix must be diagonal constant */
  bool mustBeDiagCst() const override { return false; }

  static MatrixRectangular* createFromVVD(const VectorVectorDouble& X, bool sparse = false);
  static MatrixRectangular* createFromVD(const VectorDouble &X,
                                         int nrow,
                                         int ncol,
                                         bool byCol = false,
                                         bool sparse = false);

  /*! Adding a Row or a Column (at the bottom of Rectangular Matrix) */
  void addRow(int nrow_added=1);
  void addColumn(int ncolumn_added = 1);

  MatrixRectangular* reduce(const VectorInt& validRows, const VectorInt& validCols) const;

protected:
#ifndef SWIG
  virtual double& _getValueRef(int irow, int icol) override;

private:
  bool   _isCompatible(const AMatrix& m) const override { return isSameSize(m); }
  double _getValue(int irow, int icol) const override;
  double _getValue(int irank) const override;
  void   _setValue(int rank, double value) override;
  void   _setValue(int irow, int icol, double value) override;
  void   _setValues(const double* values, bool byCol=true) override;
  int    _getMatrixSize() const override;
  void   _transposeInPlace() override;
  void   _allocate() override;
  void   _deallocate() override;
  void   _prodVector(const double *inv,double *outv) const override;
  int    _invert() override;
  int    _solve(const VectorDouble& b, VectorDouble& x) const override;

  int    _getIndexToRank(int irow,int icol) const;

private:
  VectorDouble _rectMatrix;
#endif
};
