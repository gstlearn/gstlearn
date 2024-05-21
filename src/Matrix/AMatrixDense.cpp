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
#include "Matrix/AMatrixDense.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Matrix/AMatrix.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"

#include <math.h>

AMatrixDense::AMatrixDense(int nrow, int ncol, int opt_eigen)
  : AMatrix(nrow, ncol, opt_eigen),
    _flagEigenDecompose(false),
    _eigenValues(),
    _eigenVectors(nullptr),
    _eigenMatrix()
{
  _allocate();
}

AMatrixDense::AMatrixDense(const AMatrixDense &r)
  : AMatrix(r),
    _flagEigenDecompose(false),
    _eigenValues(),
    _eigenVectors(nullptr),
    _eigenMatrix()
{
  // operator "=" is faster than the copy constructor
  // (see https://stackoverflow.com/questions/47644021/eigen-copy-constructor-vs-operator-performance)
  _allocate();
  _recopy(r);
}

AMatrixDense::AMatrixDense(const AMatrix &m)
    : AMatrix(m),
      _flagEigenDecompose(false),
      _eigenValues(),
      _eigenVectors(nullptr),
      _eigenMatrix()
{
  if (m.empty())
  {
    messerr("The input matrix should be non-empty");
    _clear();
    return;
  }
  _allocate();
  if (isFlagEigen())
    copyElements(m);
}

AMatrixDense& AMatrixDense::operator= (const AMatrixDense &r)
{
  if (this != &r)
  {
    AMatrix::operator=(r);
    _recopy(r);
  }
  return *this;
}

AMatrixDense::~AMatrixDense()
{
  if (_eigenVectors != nullptr)
  {
    delete _eigenVectors;
    _eigenVectors = nullptr;
  }
}

bool AMatrixDense::_isNumberValid(int nrows, int ncols) const
{
  AMatrix::_isNumbersValid(nrows,ncols);
  if (nrows != ncols)
  {
    messerr("Arguments 'nrows' and 'ncols' should be equal for Square Matrices");
    return false;
  }
  return true;
}

void AMatrixDense::_allocate()
{
  if (isFlagEigen())
  {
    if (isMultiThread()) omp_set_num_threads(getMultiThread()); // TODO Move to multithread handling class
    _eigenMatrix = Eigen::MatrixXd::Constant(getNRows(),getNCols(),0.);
  }
  else
  {
    _allocate_();
  }
}

void AMatrixDense::_deallocate()
{
  if (isFlagEigen())
  {
    delete _eigenVectors;
    _eigenVectors = nullptr;
  }
  else
    _deallocate_();
}

double AMatrixDense::getValue(int irow, int icol, bool flagCheck) const
{
  if (flagCheck && ! _isIndexValid(irow, icol)) return TEST;
  if (isFlagEigen())
    return _eigenMatrix(irow, icol);
  else
    return _getValue(irow, icol);
}

double AMatrixDense::_getValueByRank(int irank) const
{
  if (isFlagEigen())
    return *(_eigenMatrix.data() + irank);
  else
    return _getValueByRank_(irank);
  return TEST;
}

void AMatrixDense::_setValueByRank(int irank, double value)
{
  if (isFlagEigen())
    *(_eigenMatrix.data() + irank) = value;
  else
    _setValueByRank_(irank, value);
}

void AMatrixDense::setValue(int irow, int icol, double value, bool flagCheck)
{
  if (flagCheck && ! _isIndexValid(irow, icol)) return;
  if (isFlagEigen())
  {
    _eigenMatrix(irow, icol) = value;
    if (mustBeSymmetric() && irow != icol) _eigenMatrix(icol, irow) = value;
  }
  else
    _setValue(irow, icol, value);
}

void AMatrixDense::updValue(int irow,
                            int icol,
                            const EOperator &oper,
                            double value,
                            bool flagCheck)
{
  if (flagCheck && ! _isIndexValid(irow, icol)) return;
  if (isFlagEigen())
  {
    _eigenMatrix(irow, icol) = modifyOperator(oper, _eigenMatrix(irow, icol), value);
    if (mustBeSymmetric() && irow != icol)
      _eigenMatrix(icol, irow) = modifyOperator(oper, _eigenMatrix(icol, irow), value);
  }
  else
    _updValue(irow, icol, oper, value);
}

double& AMatrixDense::_getValueRef(int irow, int icol)
{
  if (isFlagEigen())
    return *(_eigenMatrix.data() + _getIndexToRank(irow, icol));
  else
    return _getValueRef_(irow, icol);
}

int AMatrixDense::_getMatrixPhysicalSize() const
{
  if (isFlagEigen())
    return _eigenMatrix.size();
  else
    return _getMatrixPhysicalSize_();
  return ITEST;
}

// Default storage in Eigen is column-major (see https://eigen.tuxfamily.org/dox/group__TopicStorageOrders.html)
int AMatrixDense::_getIndexToRank(int irow, int icol) const
{
  if (isFlagEigen())
    // Default storage in Eigen is column-major (see https://eigen.tuxfamily.org/dox/group__TopicStorageOrders.html)
    return (icol * getNRows() + irow);
  else
    return _getIndexToRank_(irow, icol);
  return ITEST;
}

void AMatrixDense::_transposeInPlace()
{
  if (isFlagEigen())
    _eigenMatrix.transposeInPlace();
  else
    _transposeInPlace_();
}

void AMatrixDense::_prodMatVecInPlacePtr(const double *x, double *y, bool transpose) const
{
  if (isFlagEigen())
  {
    Eigen::Map<const Eigen::VectorXd> xm(x, getNCols());
    Eigen::Map<Eigen::VectorXd> ym(y, getNRows());
    if (transpose)
      ym.noalias() = _eigenMatrix.transpose() * xm;
    else
      ym.noalias() = _eigenMatrix * xm;
  }
  else
    _prodMatVecInPlacePtr_(x, y, transpose);
}

void AMatrixDense::_prodVecMatInPlacePtr(const double *x,double *y, bool transpose) const
{
  if (isFlagEigen())
  {
    Eigen::Map<const Eigen::VectorXd> xm(x, getNRows());
    Eigen::Map<Eigen::VectorXd> ym(y, getNCols());
    if (transpose)
      ym.noalias() = xm.transpose() * _eigenMatrix.transpose();
    else
      ym.noalias() = xm.transpose() * _eigenMatrix;
  }
  else
    _prodVecMatInPlacePtr_(x, y, transpose);
}

int AMatrixDense::_invert()
{
  if (isFlagEigen())
    _eigenMatrix = _eigenMatrix.inverse();
  else
    my_throw("_invert should never be called here");
  return 0;
}

int AMatrixDense::_solve(const VectorDouble &b, VectorDouble &x) const
{
  if (isFlagEigen())
  {
    Eigen::Map<const Eigen::VectorXd> bm(b.data(), getNCols());
    Eigen::Map<Eigen::VectorXd> xm(x.data(), getNRows());
    xm = _eigenMatrix.inverse() * bm;
  }
  else
    my_throw("_solve should never be called here");
  return 0;
}

void AMatrixDense::setColumn(int icol, const VectorDouble& tab)
{
  if (isFlagEigen())
  {
    Eigen::Map<const Eigen::VectorXd> tabm(tab.data(), getNRows());
    _eigenMatrix.col(icol) = tabm;
  }
  else
    AMatrix::setColumn(icol, tab);
}

void AMatrixDense::setRow(int irow, const VectorDouble& tab)
{
  if (isFlagEigen())
  {
    Eigen::Map<const Eigen::VectorXd> tabm(tab.data(), getNCols());
    _eigenMatrix.row(irow) = tabm;
  }
  else
    AMatrix::setRow(irow, tab);
}

void AMatrixDense::setDiagonal(const VectorDouble& tab)
{
  if (isFlagEigen())
  {
    _eigenMatrix.setZero();
    Eigen::Map<const Eigen::VectorXd> tabm(tab.data(), getNRows());
    _eigenMatrix.diagonal() = tabm;
  }
  else
    AMatrix::setDiagonal(tab);
}

void AMatrixDense::setDiagonalToConstant(double value)
{
  if (isFlagEigen())
  {
    _eigenMatrix.setZero();
    _eigenMatrix.diagonal() = Eigen::VectorXd::Constant(getNRows(),value);
  }
  else
    AMatrix::setDiagonalToConstant(value);
}

void AMatrixDense::addScalar(double v)
{
  if (isFlagEigen())
    _eigenMatrix.array() += v;
  else
    AMatrix::addScalar(v);
}

void AMatrixDense::addScalarDiag(double v)
{
  if (isFlagEigen())
    _eigenMatrix.diagonal() += Eigen::VectorXd::Constant(getNRows(),v);
  else
    AMatrix::addScalarDiag(v);
}

void AMatrixDense::prodScalar(double v)
{
  if (isFlagEigen())
    _eigenMatrix.array() *= v;
  else
    AMatrix::prodScalar(v);
}

void AMatrixDense::addMatInPlace(const AMatrixDense& y, double cx, double cy)
{
  if (isFlagEigen() && y.isFlagEigen())
    _eigenMatrix.noalias() = cx * _eigenMatrix + cy * y._eigenMatrix;
  else
    AMatrix::addMatInPlace(y, cx, cy);
}

void AMatrixDense::prodMatMatInPlace(const AMatrix* x,
                                     const AMatrix* y,
                                     bool transposeX,
                                     bool transposeY)
{
  const AMatrixDense* xm = dynamic_cast<const AMatrixDense*>(x);
  const AMatrixDense* ym = dynamic_cast<const AMatrixDense*>(y);
  if (xm != nullptr && ym != nullptr && isFlagEigen() &&
      xm->isFlagEigen() && ym->isFlagEigen())
  {
    if (transposeX)
    {
      if (transposeY)
      {
        _eigenMatrix.noalias() = xm->_eigenMatrix.transpose() * ym->_eigenMatrix.transpose();
      }
      else
      {
        _eigenMatrix.noalias() = xm->_eigenMatrix.transpose() * ym->_eigenMatrix;
      }
    }
    else
    {
      if (transposeY)
      {
        _eigenMatrix.noalias() = xm->_eigenMatrix * ym->_eigenMatrix.transpose();
      }
      else
      {
        _eigenMatrix.noalias() = xm->_eigenMatrix * ym->_eigenMatrix;
      }
    }
  }
  else
  {
    AMatrix::prodMatMatInPlace(x, y, transposeX, transposeY);
  }
}

/**
 * Product of matrices: 'a' * 'm' (possibly transposed) stored in 'this'
 *
 * @param a First input matrix
 * @param m Second input matrix
 * @param transpose True if 'm' should be transposed beforehand
 *
 * @note: 'a' and 'm' may NOT coincide with 'this'
 */
void AMatrixDense::prodNormMatMatInPlace(const AMatrixDense &a,
                                         const AMatrixDense &m,
                                         bool transpose)
{
  if (isFlagEigen() && a.isFlagEigen() && m.isFlagEigen())
  {
    if (transpose)
    {
      _eigenMatrix.noalias() = a._eigenMatrix.transpose() * m._eigenMatrix * a._eigenMatrix;
    }
    else
    {
      _eigenMatrix.noalias() = a._eigenMatrix * m._eigenMatrix * a._eigenMatrix.transpose();
    }
  }
  else
    AMatrix::prodNormMatMatInPlace(a, m, transpose);
}

/**
 * Product 't(A)' %*% ['vec'] %*% 'A' or 'A' %*% ['vec'] %*% 't(A)' stored in 'this'
 *
 * @param a Input matrix
 * @param vec Input vector
 * @param transpose When True, the input Matrix is transposed
 */
void AMatrixDense::prodNormMatInPlace(const AMatrixDense &a, const VectorDouble& vec, bool transpose)
{
  if (isFlagEigen() && a.isFlagEigen())
  {
    if (transpose)
     {
       if (vec.empty())
         _eigenMatrix.noalias() = a._eigenMatrix.transpose() * a._eigenMatrix;
       else
       {
         Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), vec.size());
         _eigenMatrix.noalias() = a._eigenMatrix.transpose() * vecm * a._eigenMatrix;
       }
     }
     else
     {
       if (vec.empty())
         _eigenMatrix.noalias() = a._eigenMatrix * a._eigenMatrix.transpose();
       else
       {
         Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), vec.size());
         _eigenMatrix.noalias() = a._eigenMatrix * vecm * a._eigenMatrix.transpose();
       }
     }
  }
  else
    AMatrix::prodNormMatInPlace(a, vec, transpose);
}

void AMatrixDense::fill(double value)
{
  if (isFlagEigen())
    _eigenMatrix.setConstant(value);
  else
    AMatrix::fill(value);
}

/*! Multiply a Matrix row-wise */
void AMatrixDense::multiplyRow(const VectorDouble& vec)
{
  if (isFlagEigen())
  {
    Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), getNCols());
    _eigenMatrix = vecm.asDiagonal() * _eigenMatrix;
  }
  else
  {
    AMatrix::multiplyRow(vec);
  }
}

/*! Multiply a Matrix column-wise */
void AMatrixDense::multiplyColumn(const VectorDouble& vec)
{
  if (isFlagEigen())
  {
    Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), getNRows());
    _eigenMatrix = _eigenMatrix * vecm.asDiagonal();
  }
  else
  {
    AMatrix::multiplyColumn(vec);
  }
}

/*! Divide a Matrix row-wise */
void AMatrixDense::divideRow(const VectorDouble& vec)
{
  if (isFlagEigen())
  {
    VectorDouble temp = VH::inverse(vec);
    Eigen::Map<const Eigen::VectorXd> vecm(temp.data(), getNCols());
    _eigenMatrix = vecm.asDiagonal() * _eigenMatrix;
  }
  else
  {
    AMatrix::divideRow(vec);
  }
}

/*! Divide a Matrix column-wise */
void AMatrixDense::divideColumn(const VectorDouble& vec)
{
  if (isFlagEigen())
  {
    VectorDouble temp = VH::inverse(vec);
    Eigen::Map<const Eigen::VectorXd> vecm(temp.data(), getNRows());
    _eigenMatrix = _eigenMatrix * vecm.asDiagonal();
  }
  else
  {
    AMatrix::divideColumn(vec);
  }
}

/*! Perform 'vec' * 'this' */
VectorDouble AMatrixDense::prodVecMat(const VectorDouble& x, bool transpose) const
{
  if (isFlagEigen())
  {
    Eigen::Map<const Eigen::VectorXd> xm(x.data(), getNRows());
    Eigen::VectorXd ym;
    if (transpose)
      ym = xm.transpose() * _eigenMatrix.transpose();
    else
      ym = xm.transpose() * _eigenMatrix;
    VectorDouble y(ym.data(), ym.data() + ym.size());
    return y;
  }
  else
    return AMatrix::prodVecMat(x, transpose);
}

/*! Perform 'this' * 'vec' */
VectorDouble AMatrixDense::prodMatVec(const VectorDouble& x, bool transpose) const
{
  if (isFlagEigen())
  {
    Eigen::Map<const Eigen::VectorXd> xm(x.data(), getNCols());
    Eigen::VectorXd ym;
    if (transpose)
      ym = _eigenMatrix.transpose() * xm;
    else
      ym = _eigenMatrix * xm;
    VectorDouble y(ym.data(), ym.data() + ym.size());
    return y;
  }
  else
    return AMatrix::prodMatVec(x, transpose);
}

/*! Extract a Row */
VectorDouble AMatrixDense::getRow(int irow) const
{
  if (isFlagEigen())
  {
    Eigen::VectorXd resm = _eigenMatrix.row(irow);
    return VectorDouble(resm.data(), resm.data() + resm.size());
  }
  else
  {
    return AMatrix::getRow(irow);
  }
}

/*! Extract a Column */
VectorDouble AMatrixDense::getColumn(int icol) const
{
  if (isFlagEigen())
  {
    Eigen::VectorXd resm = _eigenMatrix.col(icol);
    return VectorDouble(resm.data(), resm.data() + resm.size());
  }
  else
  {
    return AMatrix::getColumn(icol);
  }
}

int AMatrixDense::_terminateEigen(const Eigen::VectorXd &eigenValues,
                                  const Eigen::MatrixXd &eigenVectors,
                                  bool optionPositive,
                                  bool changeOrder)
{
  int nrows = getNRows();
  int ncols = getNCols();

  _eigenValues = VectorDouble(nrows);
  Eigen::Map<Eigen::VectorXd>(&_eigenValues[0], eigenValues.size()) = eigenValues;

  if (changeOrder)
    std::reverse(_eigenValues.begin(), _eigenValues.end());

  if (_eigenVectors != nullptr) delete _eigenVectors;

  VectorDouble vec(nrows * ncols);
  Eigen::Map<Eigen::MatrixXd>(&vec[0], nrows, ncols) = eigenVectors;

  _eigenVectors = MatrixSquareGeneral::createFromVD(vec, nrows, false, 1, changeOrder);

  if (optionPositive) _eigenVectors->makePositiveColumn();

  _flagEigenDecompose = true;

  return 0;
}

int AMatrixDense::_computeGeneralizedEigen(const MatrixSquareSymmetric &b, bool optionPositive)
{
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(_eigenMatrix, b._eigenMatrix);
  Eigen::VectorXd eigenValues  = solver.eigenvalues().real();
  Eigen::MatrixXd eigenVectors = solver.eigenvectors().real();

  return _terminateEigen(eigenValues, eigenVectors, optionPositive, true);
}

int AMatrixDense::_computeEigen(bool optionPositive)
{
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(_eigenMatrix);
  Eigen::VectorXd eigenValues  = solver.eigenvalues().real();
  Eigen::MatrixXd eigenVectors = solver.eigenvectors().real();

  return _terminateEigen(eigenValues, eigenVectors, optionPositive, true);
}

/// ====================================================================
/// The subsequent methods rely on the specific storage in Eigen Library
/// They are valid whatever the format of Dense matrix.
/// ====================================================================
void AMatrixDense::_recopy(const AMatrixDense &r)
{
  _eigenMatrix = r._eigenMatrix;
  _flagEigenDecompose = r._flagEigenDecompose;
  _eigenValues = r._eigenValues;
  delete _eigenVectors;
  _eigenVectors = nullptr;
  if (r._eigenVectors != nullptr)
  {
    _eigenVectors = r._eigenVectors->clone();
  }
}

