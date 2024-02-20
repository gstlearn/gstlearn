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
#include "Matrix/AMatrixSquare.hpp"

class AMatrix;

/**
 * Square Symmetric matrices are stored as Lower Triangular matrices stored by column
 */
class GSTLEARN_EXPORT MatrixSquareSymmetric : public AMatrixSquare {

public:
  MatrixSquareSymmetric(int nrow = 0, int opt_eigen = -1);
  MatrixSquareSymmetric(const MatrixSquareSymmetric &m);
  MatrixSquareSymmetric(const AMatrix &m);
  MatrixSquareSymmetric& operator= (const MatrixSquareSymmetric &r);
	virtual ~MatrixSquareSymmetric();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// ICloneable interface
  IMPLEMENT_CLONING(MatrixSquareSymmetric)

  /// Interface to AMatrix
  /*! Say if the matrix must be symmetric */
  bool mustBeSymmetric() const final { return true; }
  /*! Say if the matrix must be diagonal */
  virtual bool mustBeDiagonal() const override { return false; }
  /*! Say if the matrix must be diagonal constant */
  virtual bool mustBeDiagCst() const override { return false; }
  /// Is the matrix symmetrical ?
  bool isSymmetric(bool printWhyNot = false) const final { DECLARE_UNUSED(printWhyNot); return true; }

  void normMatrix(const AMatrix& y, const AMatrixSquare& x = AMatrixSquare(), bool transpose = false);

  static MatrixSquareSymmetric* createFromVVD(const VectorVectorDouble &X, int opt_eigen = -1);
  static MatrixSquareSymmetric* createFromVD(const VectorDouble &X,
                                             int nrow,
                                             int opt_eigen = -1);
  MatrixSquareSymmetric* createReduce(const VectorInt &validRows) const;

  int computeEigen(bool optionPositive = true);
  int computeGeneralizedEigen(const MatrixSquareSymmetric& b, bool optionPositive = true);

private:
  /// Interface for AMatrix
  virtual bool   _isCompatible(const AMatrix& m) const override { return (isSameSize(m) && isSymmetric()); }
  virtual int    _getMatrixPhysicalSize() const override;

  virtual double& _getValueRef(int irow, int icol) override;
  virtual bool    _isPhysicallyPresent(int irow, int icol) const override;
  virtual int     _getIndexToRank(int irow,int icol) const override;
  virtual void    _allocate() override;
  virtual void    _deallocate() override;
  virtual double  _getValue(int irow, int icol) const override;
  virtual double  _getValueByRank(int irank) const override;
  virtual void    _setValue(int irow, int icol, double value) override;
  virtual void    _setValueByRank(int irank, double value) override;
  virtual void    _setValues(const double* values, bool byCol = true) override;

  virtual void    _transposeInPlace() override { return ; } // Nothing to do
  virtual void    _prodMatVecInPlacePtr(const double *x,double *y, bool transpose = false) const override;
  virtual void    _prodVecMatInPlacePtr(const double *x,double *y, bool transpose = false) const override;
  virtual int     _invert() override;
  virtual int     _solve(const VectorDouble& b, VectorDouble& x) const override;

private:
  // The subsequent methods rely on the specific local storage ('squareSymMatrix')
  void    _recopyLocal(const MatrixSquareSymmetric& r);
  double  _getValueLocal(int irow, int icol) const;
  double  _getValueLocal(int irank) const;
  double& _getValueRefLocal(int irow, int icol);
  void    _setValueLocal(int irow, int icol, double value);
  void    _setValueLocal(int irank, double value);
  void    _prodMatVecInPlacePtrLocal(const double *x, double *y, bool transpose = false) const;
  void    _prodVecMatInPlacePtrLocal(const double *x, double *y, bool transpose = false) const;
  void    _setValuesLocal(const double *values, bool byCol);
  int     _invertLocal();
  void    _allocateLocal();
  int     _getIndexToRankLocal(int irow, int icol) const;
  int     _getMatrixPhysicalSizeLocal() const;
  int     _solveLocal(const VectorDouble& b, VectorDouble& x) const;
  int     _computeEigenLocal(bool optionPositive = true);
  int     _computeGeneralizedEigenLocal(const MatrixSquareSymmetric& b, bool optionPositive = true);

  // Local functions (old style algebra)
  int _matrix_geigen(const double *a,
                     const double *b,
                     int neq,
                     double *value,
                     double *vector) const;
  void _matrix_triangular_product(int neq,
                                  int mode,
                                  const double *al,
                                  const double *b,
                                  double *x) const;

private:
  VectorDouble _squareSymMatrix; // Classical storage
};
