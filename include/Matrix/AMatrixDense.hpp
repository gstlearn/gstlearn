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

#include <Eigen/Dense>

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
  virtual bool isMatrixDense() const { return true; }

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
  /*! Perform M * 'vec' */
  virtual VectorDouble prodVector(const VectorDouble& vec) const override;
  /*! Perform 'vec'^T * M */
  virtual VectorDouble prodTVector(const VectorDouble& vec) const override;
  /*! Extract a Row */
  virtual VectorDouble getRow(int irow) const override;
  /*! Extract a Column */
  virtual VectorDouble getColumn(int icol) const override;

  /// The next functions use specific definition of matrix (to avoid dynamic_cast)
  /// rather than manipulating AMatrix. They are no more generic of AMatrix
  /// WARNING: output matrix should not match any of input matrices (speed up).
  /*! Add a matrix (multiplied by a constant) */
  virtual void addMatrix(const AMatrixDense& y, double value = 1.);
  /*! Multiply a matrix by another and store the result in the current matrix */
  virtual void prodMatrix(const AMatrixDense &x,
                          const AMatrixDense &y,
                          bool transposeX = false,
                          bool transposeY = false);
  /*! Linear combination of matrices */
  virtual void linearCombination(double cx, double cy, const AMatrixDense& y);

protected:
  virtual int     _getMatrixPhysicalSize() const override;
  virtual double& _getValueRef(int irow, int icol) override;

  virtual void    _allocate() override;
  virtual void    _deallocate() override;
  virtual double  _getValueByRank(int rank) const override;
  virtual double  _getValue(int irow, int icol) const override;
  virtual void    _setValueByRank(int rank, double value) override;
  virtual void    _setValue(int irow, int icol, double value) override;
  virtual int     _getIndexToRank(int irow,int icol) const override;

  virtual void    _transposeInPlace() override;
  virtual void    _prodVectorInPlace(const double *inv,double *outv) const override;
  virtual int     _invert() override;
  virtual int     _solve(const VectorDouble& b, VectorDouble& x) const override;

  bool            _isNumberValid(int nrows,int ncols) const;

private:
  /// =========================================================================
  /// The subsequent methods rely on the specific local storage ('eigenMatrix')
  /// =========================================================================
  void    _recopyLocal(const AMatrixDense &r);
  void    _allocateLocal();
  int     _solveLocal(const VectorDouble &b, VectorDouble &x) const;
  int     _invertLocal();
  void    _prodVectorLocal(const double *inv, double *outv) const;
  void    _transposeInPlaceLocal();
  int     _getIndexToRankLocal(int irow, int icol) const;
  int     _getMatrixPhysicalSizeLocal() const;
  double& _getValueRefLocal(int irow, int icol);
  void    _setValueLocal(int irow, int icol, double value);
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
  void _addMatrixLocal(const AMatrixDense& y, double value = 1.);
  void _prodMatrixLocal(const AMatrixDense &x,
                        const AMatrixDense &y,
                        bool transposeX = false,
                        bool transposeY = false);
  void _linearCombinationLocal(double cx, double cy, const AMatrixDense& y);
  void _fillLocal(double value);
  void _multiplyRowLocal(const VectorDouble& vec);
  void _multiplyColumnLocal(const VectorDouble& vec);
  void _divideRowLocal(const VectorDouble& vec);
  void _divideColumnLocal(const VectorDouble& vec);
  VectorDouble _prodVectorLocal(const VectorDouble& vec) const;
  VectorDouble _prodTVectorLocal(const VectorDouble& vec) const;
  VectorDouble _getRowLocal(int irow) const;
  VectorDouble _getColumnLocal(int icol) const;

public:
  Eigen::MatrixXd _eigenMatrix; // Eigen storage for Dense matrix in Eigen Library
};
