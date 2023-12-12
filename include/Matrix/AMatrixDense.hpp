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
  AMatrixDense(int nrow = 0, int ncol = 0);
  AMatrixDense(const AMatrixDense &m);
  AMatrixDense& operator= (const AMatrixDense &r);
	virtual ~AMatrixDense();

protected:
  virtual int     _getMatrixPhysicalSize() const override;
  virtual double& _getValueRef(int irow, int icol) override;

  virtual void    _allocate() override;
  virtual void    _deallocate() override;
  virtual double  _getValue(int rank) const override;
  virtual double  _getValue(int irow, int icol) const override;
  virtual void    _setValue(int rank, double value) override;
  virtual void    _setValue(int irow, int icol, double value) override;
  virtual int     _getIndexToRank(int irow,int icol) const override;

  virtual void    _transposeInPlace() override;
  virtual void    _prodVector(const double *inv,double *outv) const override;
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

public:
  Eigen::MatrixXd _eigenMatrix; // Eigen storage
};
