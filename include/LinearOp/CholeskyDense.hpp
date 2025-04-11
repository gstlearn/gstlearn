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

#include "LinearOp/ACholesky.hpp"
#include "Basic/VectorNumT.hpp"

#ifndef SWIG
  #include <Eigen/Core>
  #include <Eigen/Dense>
  #include <Eigen/src/Core/Matrix.h>
#endif

class MatrixSymmetric;
class MatrixDense;

class GSTLEARN_EXPORT CholeskyDense: public ACholesky
{
public:
  CholeskyDense(const MatrixSymmetric* mat = nullptr);
  CholeskyDense(const CholeskyDense &m);
  CholeskyDense& operator=(const CholeskyDense &m);
  virtual ~CholeskyDense();

  int setMatrix(const MatrixSymmetric* mat);
  double computeLogDeterminant() const override;

  VectorDouble getLowerTriangle() const;
  double getLowerTriangle(int i, int j) const;
  VectorDouble getUpperTriangleInverse() const;
  double getUpperTriangleInverse(int i, int j) const;
  void solveMatInPlace(const MatrixDense& mat, MatrixDense& res) const;
  int addSolveX(const constvect vecin, vect vecout) const override;
  int addInvLtX(const constvect vecin, vect vecout) const override;
  int addLtX(const constvect vecin, vect vecout) const override;
  int addLX(const constvect vecin, vect vecout) const override;
  int addInvLX(const constvect vecin, vect vecout) const override;

  void matProductInPlace(int mode,
                         const MatrixDense& a,
                         MatrixDense& x);
  void normMatInPlace(int mode,
                      int neq,
                      const MatrixSymmetric& a,
                      MatrixSymmetric& b);
  void clear();
  bool empty() const;
private:
  void _clear();
  int _prepare() const;
  int _getTriangleSize() const;
  int _computeTL() const;
  int _computeXL() const;

private:
  mutable VectorDouble _tl; // Lower triangular matrix
  mutable VectorDouble _xl; // Lower triangular matrix
  mutable Eigen::LLT<Eigen::MatrixXd> _factor; // Cholesky decomposition (Eigen format)
};
