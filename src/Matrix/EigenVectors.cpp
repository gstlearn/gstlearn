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
#include "Matrix/EigenVectors.hpp"

#include "Basic/Message.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixSquare.hpp"
#include "Matrix/MatrixSymmetric.hpp"

namespace gstlrn
{

EigenVectors::EigenVectors(const MatrixSquare& mat,
                           const MatrixSymmetric* b,
                           bool optionPositive)
  : _ready(false)
  , _eigenValues()
  , _eigenVectors()
  , _mat(&mat)
  , _nrows(0)
  , _ncols(0)
{
  const auto* matss = dynamic_cast<const MatrixSymmetric*>(_mat);
  if (matss == nullptr)
  {
    // The matrix is not declared as Symmetric. Check if it is actually Symmetric
    if (!matss->isSymmetric())
    {
      messerr("The Eigen Decomposition is valid for any Square Matrix");
      messerr("It is currently implemented only for Symmetric ones");
      return;
    }
  }
  _nrows = _mat->getNRows();
  _ncols = _mat->getNCols();
  _computeEigen(b, optionPositive);
  _ready = true;
}

EigenVectors::EigenVectors(const EigenVectors& r)
  : _ready(r._ready)
  , _eigenValues(r._eigenValues)
  , _eigenVectors(r._eigenVectors)
  , _mat(r._mat)
  , _nrows(r._nrows)
  , _ncols(r._ncols)
{
}

EigenVectors& EigenVectors::operator=(const EigenVectors& r)
{
  if (this != &r)
  {
    _ready        = r._ready;
    _eigenValues  = r._eigenValues;
    _eigenVectors = r._eigenVectors;
    _mat          = r._mat;
    _nrows        = r._nrows;
    _ncols        = r._ncols;
  }
  return *this;
}

EigenVectors::~EigenVectors()
{
}

void EigenVectors::_computeEigen(const MatrixSymmetric* b, bool optionPositive)
{
  Eigen::VectorXd eigenValues;
  Eigen::MatrixXd eigenVectors;

  if (b == nullptr)
  {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(_mat->eigenMat());
    eigenValues  = solver.eigenvalues().real();
    eigenVectors = solver.eigenvectors().real();
  }
  else
  {
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(_mat->eigenMat(), b->eigenMat());
    eigenValues  = solver.eigenvalues().real();
    eigenVectors = solver.eigenvectors().real();
  }

  _terminateEigen(eigenValues, eigenVectors, optionPositive, true);
}

void EigenVectors::_terminateEigen(const Eigen::VectorXd& eigenValues,
                                   const Eigen::MatrixXd& eigenVectors,
                                   bool optionPositive,
                                   bool changeOrder)
{
  _eigenValues                                                         = VectorDouble(_nrows);
  Eigen::Map<Eigen::VectorXd>(_eigenValues.data(), eigenValues.size()) = eigenValues;

  if (changeOrder)
    std::reverse(_eigenValues.begin(), _eigenValues.end());

  VectorDouble vec(_nrows * _ncols);
  Eigen::Map<Eigen::MatrixXd>(vec.data(), _nrows, _ncols) = eigenVectors;

  auto* a       = MatrixSquare::createFromVD(vec, _nrows, false, changeOrder);
  _eigenVectors = *a;
  delete a;

  if (optionPositive)
    _eigenVectors.makePositiveColumn();
}

} // namespace gstlrn
