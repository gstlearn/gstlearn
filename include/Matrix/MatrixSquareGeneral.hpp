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
  MatrixSquareGeneral(int nrow = 0, int opt_eigen=-1);
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

  static MatrixSquareGeneral* createFromVVD(const VectorVectorDouble& X, int opt_eigen = -1);
  static MatrixSquareGeneral* createFromVD(const VectorDouble &X,
                                           int nrow,
                                           bool byCol = false,
                                           int opt_eigen = -1,
                                           bool invertColumnOrder = false);

  int decomposeLU(MatrixSquareGeneral& tls,
                  MatrixSquareGeneral& tus,
                  double eps = EPSILON20);

  /// Interface for AMatrixDense
  void    _setValue(int irow, int icol, double value) override;
  double  _getValue(int irow, int icol) const override;
  void    _updValue(int irow, int icol, const EOperator& oper, double value) override;

private:
  /// Interface for AMatrix
  virtual bool    _isCompatible(const AMatrix& m) const override
  {
    return (isSameSize(m) && m.isSquare());
  }
  virtual int     _getMatrixPhysicalSize() const override;
  virtual double& _getValueRef(int irow, int icol) override;
  virtual void    _allocate() override;
  virtual void    _deallocate() override;
  virtual double  _getValueByRank(int irank) const override;
  virtual void    _setValueByRank(int irank, double value) override;

  virtual void    _transposeInPlace() override;
  virtual void    _prodMatVecInPlacePtr(const double *x,double *y, bool transpose = false) const override;
  virtual void    _prodVecMatInPlacePtr(const double *x,double *y, bool transpose = false) const override;
  virtual int     _invert() override;
  virtual int     _solve(const VectorDouble& b, VectorDouble& x) const override;

private:
  /// ==========================================================================
  /// The subsequent methods rely on the specific local storage ('squareMatrix')
  /// ==========================================================================
  void    _recopy(const MatrixSquareGeneral &r);
  int     _invertLU();
  int     _solveLU(const MatrixSquareGeneral& tus,
                   const MatrixSquareGeneral& tls,
                   const double *b,
                   double *x);
  int     _forwardLU(const MatrixSquareGeneral& tls, const double *b, double *x, double eps = EPSILON20);
  int     _backwardLU(const MatrixSquareGeneral& tus, const double *b, double *x, double eps = EPSILON20);
  int     _matrix_invreal(VectorDouble& mat, int neq);
  int     _matrix_cofactor(int neq, VectorDouble& a, VectorDouble& b);

private:
  VectorDouble _squareMatrix; // Classical storage
};

/*! Product 't(A)' %*% 'M' %*% 'A' or 'A' %*% 'M' %*% 't(A)' */
GSTLEARN_EXPORT MatrixSquareGeneral* prodNormMatMat(const AMatrixDense &a,
                                                    const AMatrixDense &m,
                                                    bool transpose = false);
/*! Product 't(A)' %*% 'A' or 'A' %*% 't(A)' */
GSTLEARN_EXPORT MatrixSquareGeneral* prodNormMat(const AMatrixDense &a,
                                                 const VectorDouble& vec = VectorDouble(),
                                                 bool transpose = false);
