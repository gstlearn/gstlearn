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

#include "Matrix/AMatrix.hpp"
#include "Basic/WarningMacro.hpp"

#ifndef SWIG
DISABLE_WARNING_PUSH
DISABLE_WARNING_COND_EXPR_CONSTANT
DISABLE_WARNING_UNUSED_BUT_SET_VARIABLE
DISABLE_WARNING_MAYBE_UNINITIALIZED
DISABLE_WARNING_DECLARATION_HIDE_GLOBAL
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
DISABLE_WARNING_POP
#endif

class MatrixSquareGeneral;
class MatrixSquareSymmetric;
class EOperator;

/**
 * Square Matrix
 */

class GSTLEARN_EXPORT AMatrixDense : public AMatrix {

public:
  AMatrixDense(int nrow = 0, int ncol = 0, int opt_eigen=-1);
  AMatrixDense(const AMatrixDense &m);
  AMatrixDense(const AMatrix &m);
  AMatrixDense& operator= (const AMatrixDense &r);
	virtual ~AMatrixDense();

  /// Interface for AMatrix
  /*! Returns if the matrix belongs to the MatrixSparse class (avoids dynamic_cast) */
  virtual bool isDense() const { return true; }
  /*! Returns if the current matrix is Sparse */
  virtual bool isSparse() const { return false; }

  /*! Set the contents of a Column */
  virtual void setColumn(int icol, const VectorDouble& tab) override;
  /*! Set the contents of a Row */
  virtual void setRow(int irow, const VectorDouble& tab) override;
  /*! Set the contents of the (main) Diagonal */
  virtual void setDiagonal(const VectorDouble& tab) override;
  /*! Set the contents of the (main) Diagonal to a constant value */
  virtual void setDiagonalToConstant(double value = 1.) override;
  /*! Add a value to each matrix component */
  virtual void addScalar(double v) override;
  /*! Add value to matrix diagonal */
  virtual void addScalarDiag(double v) override;
  /*! Multiply each matrix component by a value */
  virtual void prodScalar(double v) override;
  /*! Set all the values of the Matrix at once */
  virtual void fill(double value) override;
  /*! Multiply a Matrix row-wise */
  virtual void multiplyRow(const VectorDouble& vec) override;
  /*! Multiply a Matrix column-wise */
  virtual void multiplyColumn(const VectorDouble& vec) override;
  /*! Divide a Matrix row-wise */
  virtual void divideRow(const VectorDouble& vec) override;
  /*! Divide a Matrix column-wise */
  virtual void divideColumn(const VectorDouble& vec) override;
  /*! Perform 'vec' * 'this' */
  virtual VectorDouble prodVecMat(const VectorDouble& x, bool transpose = false) const override;
  /*! Perform 'this' * 'vec'*/
  virtual VectorDouble prodMatVec(const VectorDouble& x, bool transpose = false) const override;
  /*! Extract a Row */
  virtual VectorDouble getRow(int irow) const override;
  /*! Extract a Column */
  virtual VectorDouble getColumn(int icol) const override;
  /*! Multiply a matrix by another and stored in 'this' */
  virtual void prodMatMatInPlace(const AMatrix *x,
                                 const AMatrix *y,
                                 bool transposeX = false,
                                 bool transposeY = false) override;

  /// The next functions use specific definition of matrix (to avoid dynamic_cast)
  /// rather than manipulating AMatrix. They are not generic of AMatrix anymore.
  /// WARNING: output matrix should not match any of input matrices (speed up).
  /*! Add a matrix (multiplied by a constant) */
  void addMatInPlace(const AMatrixDense& y, double cx = 1., double cy = 1.);
  /*! Product 't(A)' %*% 'M' %*% 'A' or 'A' %*% 'M' %*% 't(A)' stored in 'this'*/
  virtual void prodNormMatMatInPlace(const AMatrixDense &a,
                                     const AMatrixDense &m,
                                     bool transpose = false);
  /*! Product 't(A)' %*% ['vec'] %*% 'A' or 'A' %*% ['vec'] %*% 't(A)' stored in 'this'*/
  virtual void prodNormMatInPlace(const AMatrixDense &a,
                                  const VectorDouble& vec = VectorDouble(),
                                  bool transpose = false);

  VectorDouble         getEigenValues()  const { return _eigenValues; }
  const MatrixSquareGeneral* getEigenVectors() const { return _eigenVectors; }

protected:
  virtual int     _getMatrixPhysicalSize() const override;
  virtual double& _getValueRef(int irow, int icol) override;

  virtual void    _allocate() override;
  virtual void    _deallocate() override;
  virtual double  _getValueByRank(int rank) const override;
  virtual double  _getValue(int irow, int icol) const override;
  virtual void    _setValueByRank(int rank, double value) override;
  virtual void    _setValue(int irow, int icol, double value) override;
  virtual void    _updValue(int irow, int icol, const EOperator& oper, double value) override;
  virtual int     _getIndexToRank(int irow,int icol) const override;

  virtual void    _transposeInPlace() override;
  virtual void    _prodMatVecInPlacePtr(const double *x,double *y, bool transpose = false) const override;
  virtual void    _prodVecMatInPlacePtr(const double *x,double *y, bool transpose = false) const override;
  virtual int     _invert() override;
  virtual int     _solve(const VectorDouble& b, VectorDouble& x) const override;

  bool            _isNumberValid(int nrows,int ncols) const;
  int             _computeEigen(bool optionPositive = true);
  int             _computeGeneralizedEigen(const MatrixSquareSymmetric& b, bool optionPositive = true);

private:
  /// =========================================================================
  /// The subsequent methods rely on the specific local storage ('eigenMatrix')
  /// =========================================================================
  void    _recopyLocal(const AMatrixDense &r);
  void    _allocateLocal();
  int     _solveEigen(const VectorDouble &b, VectorDouble &x) const;
  int     _invertEigen();
  void    _prodMatVecInPlacePtrLocal(const double *x, double *y, bool transpose = false) const;
  void    _prodVecMatInPlacePtrLocal(const double *x, double *y, bool transpose = false) const;
  void    _transposeInPlaceEigen();
  int     _getIndexToRankEigen(int irow, int icol) const;
  int     _getMatrixPhysicalSizeEigen() const;
  double& _getValueRefEigen(int irow, int icol);
  void    _setValueEigen(int irow, int icol, double value);
  void    _updValueLocal(int irow, int icol, const EOperator& oper, double value);
  void    _setValueLocal(int irank, double value);
  double  _getValueLocal(int irank) const;
  double  _getValueLocal(int irow, int icol) const;

  void _setColumnLocal(int icol, const VectorDouble& tab);
  void _setRowLocal(int irow, const VectorDouble& tab);
  void _setDiagonalLocal(const VectorDouble& tab);
  void _setDiagonalToConstantLocal(double value = 1.);
  void _addScalarLocal(double v);
  void _addScalarDiagLocal(double v);
  void _prodScalarLocal(double v);
  void _addMatInPlaceLocal(const AMatrixDense& y, double cx = 1., double cy = 1.);
  void _prodMatMatInPlaceLocal(const AMatrixDense *x,
                               const AMatrixDense *y,
                               bool transposeX = false,
                               bool transposeY = false);
  void _prodNormMatMatInPlaceLocal(const AMatrixDense& a, const AMatrixDense& m, bool transpose = false);
  void _prodNormMatInPlaceLocal(const AMatrixDense &a,
                                const VectorDouble &vec = VectorDouble(),
                                bool transpose = false);
  void _fillLocal(double value);
  void _multiplyRowLocal(const VectorDouble& vec);
  void _multiplyColumnLocal(const VectorDouble& vec);
  void _divideRowLocal(const VectorDouble& vec);
  void _divideColumnLocal(const VectorDouble& vec);
  VectorDouble _prodMatVecLocal(const VectorDouble& x, bool transpose = false) const;
  VectorDouble _prodVecMatLocal(const VectorDouble& x, bool transpose = false) const;
  VectorDouble _getRowLocal(int irow) const;
  VectorDouble _getColumnLocal(int icol) const;

  int _computeEigenLocal(bool optionPositive = true);
  int _computeGeneralizedEigenLocal(const MatrixSquareSymmetric &b, bool optionPositive = true);

protected:
  bool _flagEigenDecompose;
  VectorDouble         _eigenValues; // only when ! flag_eigen()
  MatrixSquareGeneral* _eigenVectors; // only when ! flag_eigen()

protected:
#ifndef SWIG
  Eigen::MatrixXd _eigenMatrix; // Eigen storage for Dense matrix in Eigen Library
#endif
};

