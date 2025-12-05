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
#include "Matrix/MatrixSquare.hpp"

namespace gstlrn
{
class MatrixSymmetric;

/**
 * @brief Prepare calculations for Eigen or a Generalized Eigen (if 'b' is provided)
 *
 * @param mat Square matrix used for Eigen decomposition
 * @param b  Auxiliary Square Symmetric matrix used for Generalized Eigen decomposition
 * @param optionPositive Positive flag
 *
 * @note: Eigen decomposition is valid for any Square Matrix.
 * @note: However, currently, it has only been implemented for Symmetric ones.
 * @note: Test isReady() to check if the decomposition has been correctly performed.
 */
class GSTLEARN_EXPORT EigenVectors
{
public:
  EigenVectors(const MatrixSquare& mat, const MatrixSymmetric* b = nullptr, bool optionPositive = true);
  EigenVectors(const EigenVectors& r);
  EigenVectors& operator=(const EigenVectors& r);
  virtual ~EigenVectors();

  const VectorDouble& getEigenValues() const { return _eigenValues; }
  const MatrixSquare& getEigenVectors() const { return _eigenVectors; }
  bool isReady() const { return _ready; }

private:
  void _computeEigen(const MatrixSymmetric* b, bool optionPositive = true);
  void _terminateEigen(const Eigen::VectorXd& eigenValues,
                       const Eigen::MatrixXd& eigenVectors,
                       bool optionPositive = true,
                       bool changeOrder    = false);

protected:
  bool _ready;
  VectorDouble _eigenValues;
  MatrixSquare _eigenVectors;

  const MatrixSquare& _mat; // Not to be deleted
  Id _nrows;
  Id _ncols;
};
} // namespace gstlrn