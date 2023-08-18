/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"
#include "Matrix/AMatrixSquare.hpp"

class AMatrix;

/**
 * Square Symmetric matrices are stored as Lower Triangular matrices stored by column
 */
class GSTLEARN_EXPORT MatrixSquareSymmetric : public AMatrixSquare {

public:
  MatrixSquareSymmetric(int nrow = 0);
  MatrixSquareSymmetric(const MatrixSquareSymmetric &m);
  MatrixSquareSymmetric(const AMatrix &m);
  MatrixSquareSymmetric& operator= (const MatrixSquareSymmetric &r);
	virtual ~MatrixSquareSymmetric();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// ICloneable interface
  IMPLEMENT_CLONING(MatrixSquareSymmetric)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /*! Say if the matrix must be symmetric */
  bool mustBeSymmetric() const final { return true; }
  /*! Say if the matrix must be diagonal */
  virtual bool mustBeDiagonal() const override { return false; }
  /*! Say if the matrix must be diagonal constant */
  virtual bool mustBeDiagCst() const override { return false; }

  /// TODO : isPositiveDefinite

  /// Is the matrix symmetrical ?
  bool isSymmetric(bool printWhyNot = false) const final { DECLARE_UNUSED(printWhyNot); return true; }

  void initMatTri(int nsize,double* tab);
  void normSingleMatrix(const AMatrix& x);
  void normTSingleMatrix(const AMatrix& x);

  static MatrixSquareSymmetric* createFromVVD(const VectorVectorDouble& X);
  MatrixSquareSymmetric* reduce(const VectorInt &validRows) const;

protected:
#ifndef SWIG
  virtual double& _getValueRef(int irow, int icol) override;

private:
  virtual bool   _isCompatible(const AMatrix& m) const override { return (isSameSize(m) && isSymmetric()); }
  virtual bool   _isPhysicallyPresent(int irow, int icol) const override;
  virtual int    _getIndexToRank(int irow,int icol) const override;
  virtual double _getValue(int irow, int icol) const override;
  virtual double _getValue(int irank) const override;
  virtual void   _setValue(int irow, int icol, double value) override;
  virtual void   _setValue(int irank, double value) override;
  virtual void   _transposeInPlace() override { return ; } // Nothing to do
  virtual void   _setValues(const double* values, bool byCol = true) override;
  virtual int    _getMatrixSize() const override;
  virtual void   _allocate() override;
  virtual void   _deallocate() override;
  virtual void   _prodVector(const double *inv,double *outv) const override;
  virtual int    _invert() override;
  virtual int    _solve(const VectorDouble& b, VectorDouble& x) const override;

  void   _recopy(const MatrixSquareSymmetric &r);

private:
  VectorDouble _squareSymMatrix;
#endif
};
