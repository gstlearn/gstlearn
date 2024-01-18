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

/**
 * Square Matrix General
 */
class GSTLEARN_EXPORT MatrixSquareGeneral : public AMatrixSquare {

public:
  MatrixSquareGeneral(int nrow = 0, int opt_eigen=1);
  MatrixSquareGeneral(const MatrixSquareGeneral &m);
  MatrixSquareGeneral(const AMatrix &m);
  MatrixSquareGeneral& operator= (const MatrixSquareGeneral &r);
	virtual ~MatrixSquareGeneral();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// ICloneable interface
  IMPLEMENT_CLONING(MatrixSquareGeneral)

  /// Interface for AMatrix
  /*! Say if the matrix must be symmetric */
  bool mustBeSymmetric() const override { return false; }
  /*! Say if the matrix must be diagonal */
  bool mustBeDiagonal() const override { return false; }
  /*! Say if the matrix must be diagonal constant */
  bool mustBeDiagCst() const override { return false; }
  /*! Check if the matrix is (non empty) square */
  bool isSquare(bool printWhyNot = false) const override { DECLARE_UNUSED(printWhyNot); return 1; }

  static MatrixSquareGeneral* createFromVVD(const VectorVectorDouble& X);
  MatrixSquareGeneral* reduce(const VectorInt &validRows) const;

private:
  /// Interface for AMatrix
  virtual bool    _isCompatible(const AMatrix& m) const override { return (isSameSize(m) && isSquare()); }
  virtual int     _getMatrixPhysicalSize() const override;

  virtual double& _getValueRef(int irow, int icol) override;
  virtual void    _allocate() override;
  virtual void    _deallocate() override;
  virtual double  _getValue(int irow, int icol) const override;
  virtual double  _getValueByRank(int irank) const override;
  virtual void    _setValue(int irow, int icol, double value) override;
  virtual void    _setValueByRank(int irank, double value) override;

  virtual void    _transposeInPlace() override;
  virtual void    _prodVectorInPlace(const double *inv,double *outv) const override;
  virtual int     _invert() override;
  virtual int     _solve(const VectorDouble& b, VectorDouble& x) const override;

private:
  /// ==========================================================================
  /// The subsequent methods rely on the specific local storage ('squareMatrix')
  /// ==========================================================================
  void    _allocateLocal();
  void    _recopyLocal(const MatrixSquareGeneral &r);
  double  _getValueLocal(int irow, int icol) const;
  double  _getValueLocal(int irank) const;
  double& _getValueRefLocal(int irow, int icol);
  void    _setValueLocal(int irow, int icol, double value);
  void    _setValueLocal(int irank, double value);
  void    _prodVectorLocal(const double *inv, double *outv) const;
  void    _transposeInPlaceLocal();
  int     _invertLocal();

private:
  VectorDouble _squareMatrix; // Classical storage
};
