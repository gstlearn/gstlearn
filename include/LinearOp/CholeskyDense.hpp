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

class MatrixSquareSymmetric;
class MatrixRectangular;

class GSTLEARN_EXPORT CholeskyDense: public ACholesky
{
public:
  CholeskyDense(const MatrixSquareSymmetric* mat = nullptr);
  CholeskyDense(const CholeskyDense &m) = delete;
  CholeskyDense& operator=(const CholeskyDense &m) = delete;
  virtual ~CholeskyDense();

  int setMatrix(const MatrixSquareSymmetric* mat);
  double computeLogDeterminant() const override;

  VectorDouble getCholeskyTL() const;
  double getCholeskyTL(int i, int j) const;
  double getCholeskyTL(int iad) const;
  VectorDouble getCholeskyXL() const;
  double getCholeskyXL(int i, int j) const;

  int addSolveX(const constvect vecin, vect vecout) const override;
  int addInvLtX(const constvect vecin, vect vecout) const override;
  int addLtX(const constvect vecin, vect vecout) const override;
  int addLX(const constvect vecin, vect vecout) const override;
  int addInvLX(const constvect vecin, vect vecout) const override;

  void productCholeskyInPlace(int mode,
                              const MatrixRectangular& a,
                              MatrixRectangular& x);
  void normCholeskyInPlace(int mode,
                           int neq,
                           const MatrixSquareSymmetric& a,
                           MatrixSquareSymmetric& b);

private:
  void _clear();
  int _prepare() const;
  int _getTriangleSize() const;
  int _computeTL() const;
  int _computeXL() const;

private:
  mutable VectorDouble _tl; // Lower triangular matrix
  mutable VectorDouble _xl; // Lower triangular matrix
  mutable Eigen::LLT<Eigen::MatrixXd>* _factor; // Cholesky decomposition (Eigen format)
};
