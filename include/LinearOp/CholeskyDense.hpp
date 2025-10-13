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
#include "Matrix/MatrixSymmetric.hpp"
#include "LinearOp/ACholesky.hpp"
#include "Basic/VectorNumT.hpp"

#ifndef SWIG
#  include <Eigen/Core>
#  include <Eigen/Dense>
#  include <Eigen/src/Core/Matrix.h>
#endif

namespace gstlrn
{ 
class MatrixDense;

class GSTLEARN_EXPORT CholeskyDense: public ACholesky
{
public:
  CholeskyDense(const MatrixSymmetric& mat = MatrixSymmetric());
  CholeskyDense(const CholeskyDense& m);
  CholeskyDense& operator=(const CholeskyDense& m);
  virtual ~CholeskyDense();

  Id setMatrix(const MatrixSymmetric& mat);
  double computeLogDeterminant() const override;

  VectorDouble getLowerTriangle() const;
  double getLowerTriangle(Id i, Id j) const;
  VectorDouble getUpperTriangleInverse() const;
  double getUpperTriangleInverse(Id i, Id j) const;
  void solveMatInPlace(const MatrixDense& mat, MatrixDense& res) const;
  Id addSolveX(const constvect vecin, vect vecout) const override;
  Id addInvLtX(const constvect vecin, vect vecout) const override;
  Id addLtX(const constvect vecin, vect vecout) const override;
  Id addLX(const constvect vecin, vect vecout) const override;
  Id addInvLX(const constvect vecin, vect vecout) const override;

  void matProductInPlace(Id mode,
                         const MatrixDense& a,
                         MatrixDense& x);
  void normMatInPlace(Id mode,
                      Id neq,
                      const MatrixSymmetric& a,
                      MatrixSymmetric& b);
  MatrixDense inverse() const;
  void clear();
  bool empty() const;

private:
  void _clear();
  Id _prepare(const MatrixSymmetric& mat) const;
  Id _getTriangleSize() const;
  Id _computeTL() const;
  Id _computeXL() const;

private:
  mutable VectorDouble _tl;                    // Lower triangular matrix
  mutable VectorDouble _xl;                    // Lower triangular matrix
  mutable Eigen::LLT<Eigen::MatrixXd> _factor; // Cholesky decomposition (Eigen format)
};
}