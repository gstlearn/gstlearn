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
#include "Basic/Utilities.hpp"
#include "geoslib_define.h"

#include <math.h>
#include <omp.h>

AMatrixDense::AMatrixDense(int nrow, int ncol)
  : AMatrix(nrow, ncol),
    _flagEigenDecompose(false),
    _eigenValues(),
    _eigenVectors(nullptr),
    _eigenMatrix(std::make_unique<Eigen::Map<Eigen::MatrixXd>>(nullptr,0, 0))
{
  _allocate();
}

AMatrixDense::AMatrixDense(const AMatrixDense &r)
  : AMatrix(r),
    _flagEigenDecompose(false),
    _eigenValues(),
    _eigenVectors(nullptr),
    _eigenMatrix(std::make_unique<Eigen::Map<Eigen::MatrixXd>>(nullptr,0, 0))
{
  _allocate();
  _recopy(r);
}

AMatrixDense::AMatrixDense(const AMatrix &r)
  : AMatrix(r),
    _flagEigenDecompose(false),
    _eigenValues(),
    _eigenVectors(nullptr),
    _eigenMatrix(std::make_unique<Eigen::Map<Eigen::MatrixXd>>(nullptr,0, 0))
{
  _allocate();
  copyElements(r);
}

AMatrixDense& AMatrixDense::operator= (const AMatrixDense &r)
{
  if (this != &r)
  {
    AMatrix::operator=(r);
    _allocate();
    _recopy(r);
  }
  return *this;
}

AMatrixDense::~AMatrixDense()
{
  _deallocate();
}

void AMatrixDense::_allocate()
{
  _maxSize = getNRows() * getNCols();
  _deallocate();
  if (isMultiThread()) omp_set_num_threads(getMultiThread()); // TODO Move to multithread handling class
  int size = getNRows() * getNCols();
  double* data = new double[size];
  _eigenMatrix = std::make_unique<Eigen::Map<Eigen::MatrixXd>>(data,getNRows(),getNCols());
  _eigenMatrix->setZero();
}

void AMatrixDense::_deallocate()
{
  if (_eigenVectors != nullptr)
  {
    delete _eigenVectors;
    _eigenVectors = nullptr;
  }
  _eigenValues.clear();
  if (_eigenMatrix != nullptr)
    if (_eigenMatrix->data() != nullptr)
    {
      delete[] _eigenMatrix->data();
      _eigenMatrix = nullptr;
    }
}

double AMatrixDense::getValue(int irow, int icol, bool flagCheck) const
{
  if (flagCheck && ! _isIndexValid(irow, icol)) return TEST;
  return (*_eigenMatrix)(irow, icol);
}

double AMatrixDense::_getValueByRank(int irank) const
{
  return *(_eigenMatrix->data() + irank);
}

constvect AMatrixDense::getColumnPtr(int icol) const
{
  const auto a = _eigenMatrix->col(icol);
  size_t n = getNRows();
  return {a.data(), n};
}
void AMatrixDense::_setValueByRank(int irank, double value)
{
  *(_eigenMatrix->data() + irank) = value;
}

void AMatrixDense::setValue(int irow, int icol, double value, bool flagCheck)
{
  if (flagCheck && ! _isIndexValid(irow, icol)) return;
  (*_eigenMatrix)(irow, icol) = value;
  if (mustBeSymmetric() && irow != icol) (*_eigenMatrix)(icol, irow) = value;
}

void AMatrixDense::updValue(int irow,
                            int icol,
                            const EOperator &oper,
                            double value,
                            bool flagCheck)
{
  if (flagCheck && ! _isIndexValid(irow, icol)) return;
  double result = modifyOperator(oper, (*_eigenMatrix)(irow, icol), value);
  (*_eigenMatrix)(irow, icol) = result;
  if (mustBeSymmetric() && irow != icol)
    (*_eigenMatrix)(icol, irow) = result;
}

double& AMatrixDense::_getValueRef(int irow, int icol)
{
  return *(_eigenMatrix->data() + _getIndexToRank(irow, icol));
}

int AMatrixDense::_getMatrixPhysicalSize() const
{
  return _eigenMatrix->size();
}

// Default storage in Eigen is column-major (see https://eigen.tuxfamily.org/dox/group__TopicStorageOrders.html)
int AMatrixDense::_getIndexToRank(int irow, int icol) const
{
  return (icol * getNRows() + irow);
}

void AMatrixDense::_transposeInPlace()
{
  int nrows = getNRows();
  int ncols = getNCols();
  _eigenMatrix->transposeInPlace();
  _setNCols(nrows);
  _setNRows(ncols);
}

void AMatrixDense::_addProdMatVecInPlaceToDestPtr(const double *x,double *y, bool transpose) const
{
  Eigen::Map<const Eigen::VectorXd> xm(x, getNCols());
  Eigen::Map<Eigen::VectorXd> ym(y, getNRows());
  if (transpose)
    ym.noalias() += _eigenMatrix->transpose() * xm;
  else
    ym.noalias() += *_eigenMatrix * xm;

}

//TODO supress this method and implement it in the virtual class AMatrix
void AMatrixDense::_prodMatVecInPlacePtr(const double *x, double *y, bool transpose) const
{
  Eigen::Map<const Eigen::VectorXd> xm(x, getNCols());
  Eigen::Map<Eigen::VectorXd> ym(y, getNRows());
  if (transpose)
    ym.noalias() = _eigenMatrix->transpose() * xm;
  else
    ym.noalias() = *_eigenMatrix * xm;
}

void AMatrixDense::_prodVecMatInPlacePtr(const double *x,double *y, bool transpose) const
{
  Eigen::Map<const Eigen::VectorXd> xm(x, getNRows());
  Eigen::Map<Eigen::VectorXd> ym(y, getNCols());
  if (transpose)
    ym.noalias() = xm.transpose() * _eigenMatrix->transpose();
  else
    ym.noalias() = xm.transpose() * *_eigenMatrix;
}

int AMatrixDense::_invert()
{
  /// TODO : check beforehand if matrix is invertible ?
  *_eigenMatrix = _eigenMatrix->inverse();
  return 0;
}

int AMatrixDense::_solve(const VectorDouble &b, VectorDouble &x) const
{
  /// TODO : check beforehand if matrix is invertible ?
  Eigen::Map<const Eigen::VectorXd> bm(b.data(), getNCols());
  Eigen::Map<Eigen::VectorXd> xm(x.data(), getNRows());
  xm = _eigenMatrix->inverse() * bm;
  return 0;
}

void AMatrixDense::setColumn(int icol, const VectorDouble& tab, bool flagCheck)
{
  if (flagCheck)
  {
    if (! _isColumnValid(icol)) return;
    if (! _isColumnSizeConsistent(tab)) return;
  }
  Eigen::Map<const Eigen::VectorXd> tabm(tab.data(), getNRows());
  _eigenMatrix->col(icol) = tabm;
}

void AMatrixDense::setRow(int irow, const VectorDouble& tab, bool flagCheck)
{
  if (flagCheck)
  {
    if (! _isRowValid(irow)) return;
    if (! _isRowSizeConsistent(tab)) return;
  }
  Eigen::Map<const Eigen::VectorXd> tabm(tab.data(), getNCols());
  _eigenMatrix->row(irow) = tabm;
}

void AMatrixDense::setDiagonal(const VectorDouble& tab, bool flagCheck)
{
  if (flagCheck)
  {
    if (! _isRowSizeConsistent(tab)) return;
  }
  _eigenMatrix->setZero();
  Eigen::Map<const Eigen::VectorXd> tabm(tab.data(), getNRows());
  _eigenMatrix->diagonal() = tabm;
}

void AMatrixDense::setDiagonalToConstant(double value)
{
  _eigenMatrix->setZero();
  _eigenMatrix->diagonal() = Eigen::VectorXd::Constant(getNRows(),value);
}

void AMatrixDense::addScalar(double v)
{
  _eigenMatrix->array() += v;
}

void AMatrixDense::addScalarDiag(double v)
{
  _eigenMatrix->diagonal() += Eigen::VectorXd::Constant(getNRows(),v);
}

void AMatrixDense::prodScalar(double v)
{
  _eigenMatrix->array() *= v;
}

void AMatrixDense::addMatInPlace(const AMatrixDense& y, double cx, double cy)
{
  _eigenMatrix->noalias() = cx * *_eigenMatrix + cy * *y._eigenMatrix;
}

void AMatrixDense::prodMatMatInPlace(const AMatrix* x,
                                     const AMatrix* y,
                                     bool transposeX,
                                     bool transposeY)
{
  const AMatrixDense* xm = dynamic_cast<const AMatrixDense*>(x);
  const AMatrixDense* ym = dynamic_cast<const AMatrixDense*>(y);
  if (xm != nullptr && ym != nullptr)
  {
    if (transposeX)
    {
      if (transposeY)
      {
        _eigenMatrix->noalias() = xm->_eigenMatrix->transpose() * ym->_eigenMatrix->transpose();
      }
      else
      {
        _eigenMatrix->noalias() = xm->_eigenMatrix->transpose() * *ym->_eigenMatrix;
      }
    }
    else
    {
      if (transposeY)
      {
        _eigenMatrix->noalias() = *xm->_eigenMatrix * ym->_eigenMatrix->transpose();
      }
      else
      {
        _eigenMatrix->noalias() = *xm->_eigenMatrix * *ym->_eigenMatrix;
      }
    }
  }
  else
  {
    AMatrix::prodMatMatInPlace(x, y, transposeX, transposeY);
  }
}

/**
 * Product of matrices, stored in 'this'
 * - transpose = true: t('a') * 'm' * 'a'
 * - transpose = false:  'a' * 'm' * t('a')
 *
 * @param a First input matrix
 * @param m Second input matrix
 * @param transpose True if 'a' should be transposed beforehand
 *
 * @note: 'a' and 'm' may NOT coincide with 'this'
 */
void AMatrixDense::prodNormMatMatInPlace(const AMatrixDense* a,
                                         const AMatrixDense* m,
                                         bool transpose)
{
  if (transpose)
  {
    _eigenMatrix->noalias() = a->_eigenMatrix->transpose() * *m->_eigenMatrix * *a->_eigenMatrix;
  }
  else
  {
    _eigenMatrix->noalias() = *a->_eigenMatrix * *m->_eigenMatrix * a->_eigenMatrix->transpose();
  }
}

/**
 * Product 't(A)' %*% ['vec'] %*% 'A' or 'A' %*% ['vec'] %*% 't(A)' stored in 'this'
 *
 * @param a Input matrix
 * @param vec Input vector
 * @param transpose When True, the input Matrix is transposed
 */
void AMatrixDense::prodNormMatVecInPlace(const AMatrixDense &a, const VectorDouble& vec, bool transpose)
{
  if (transpose)
  {
    if (vec.empty())
      _eigenMatrix->noalias() = a._eigenMatrix->transpose() * *a._eigenMatrix;
    else
    {
      Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), vec.size());
      _eigenMatrix->noalias() = a._eigenMatrix->transpose() * vecm * *a._eigenMatrix;
    }
  }
  else
  {
    if (vec.empty())
      _eigenMatrix->noalias() = *a._eigenMatrix * a._eigenMatrix->transpose();
    else
    {
      Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), vec.size());
      _eigenMatrix->noalias() = *a._eigenMatrix * vecm * a._eigenMatrix->transpose();
    }
  }
}

void AMatrixDense::fill(double value)
{
  _eigenMatrix->setConstant(value);
}

/*! Multiply a Matrix row-wise */
void AMatrixDense::multiplyRow(const VectorDouble& vec)
{
  Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), getNCols());
  *_eigenMatrix = vecm.asDiagonal() * *_eigenMatrix;
}

/*! Multiply a Matrix column-wise */
void AMatrixDense::multiplyColumn(const VectorDouble& vec)
{
  Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), getNRows());
  *_eigenMatrix = *_eigenMatrix * vecm.asDiagonal();
}

/*! Divide a Matrix row-wise */
void AMatrixDense::divideRow(const VectorDouble& vec)
{
  VectorDouble temp = VH::inverse(vec);
  Eigen::Map<const Eigen::VectorXd> vecm(temp.data(), getNCols());
  *_eigenMatrix = vecm.asDiagonal() * *_eigenMatrix;
}

/*! Divide a Matrix column-wise */
void AMatrixDense::divideColumn(const VectorDouble& vec)
{
  VectorDouble temp = VH::inverse(vec);
  Eigen::Map<const Eigen::VectorXd> vecm(temp.data(), getNRows());
  *_eigenMatrix = *_eigenMatrix * vecm.asDiagonal();
}

/*! Perform 'vec' * 'this' */
VectorDouble AMatrixDense::prodVecMat(const VectorDouble& x, bool transpose) const
{
  Eigen::Map<const Eigen::VectorXd> xm(x.data(), x.size());
  VectorDouble y(transpose ? getNRows() : getNCols());
  Eigen::Map<Eigen::VectorXd> ym(y.data(), y.size());
  if (transpose)
    ym = xm.transpose() * _eigenMatrix->transpose();
  else
    ym = xm.transpose() * *_eigenMatrix;
  return y;
}

/*! Perform 'this' * 'vec' */
VectorDouble AMatrixDense::prodMatVec(const VectorDouble& x, bool transpose) const
{
  Eigen::Map<const Eigen::VectorXd> xm(x.data(), x.size());
  VectorDouble y(transpose ? getNCols() : getNRows());
  Eigen::Map<Eigen::VectorXd> ym(y.data(), y.size());
  if (transpose)
    ym = _eigenMatrix->transpose() * xm;
  else
    ym = *_eigenMatrix * xm;
  return y;
}

/*! Extract a Row */
VectorDouble AMatrixDense::getRow(int irow) const
{
  VectorDouble res(getNCols());
  for (size_t i = 0; i < res.size(); ++i)
  {
    res[i] = _eigenMatrix->row(irow)[i];
  }
  return res;
}

/*! Extract a Column */
VectorDouble AMatrixDense::getColumn(int icol) const
{
  VectorDouble res(getNRows());
  for (size_t i = 0; i < res.size(); ++i)
  {
    res[i] = _eigenMatrix->col(icol)[i];
  }
  return res;
}

constvect AMatrixDense::getViewOnColumn(int icol) const
{
  constvect res(_eigenMatrix->col(icol).data(),getNRows());
  return res;
}

vect AMatrixDense::getViewOnColumnModify(int icol) const
{
  vect res(_eigenMatrix->col(icol).data(),getNRows());
  return res;
}
int AMatrixDense::_terminateEigen(const Eigen::VectorXd &eigenValues,
                                  const Eigen::MatrixXd &eigenVectors,
                                  bool optionPositive,
                                  bool changeOrder)
{
  int nrows = getNRows();
  int ncols = getNCols();

  _eigenValues = VectorDouble(nrows);
  Eigen::Map<Eigen::VectorXd>(_eigenValues.data(), eigenValues.size()) = eigenValues;

  if (changeOrder)
    std::reverse(_eigenValues.begin(), _eigenValues.end());

  delete _eigenVectors;

  VectorDouble vec(nrows * ncols);
  Eigen::Map<Eigen::MatrixXd>(vec.data(), nrows, ncols) = eigenVectors;

  _eigenVectors = MatrixSquareGeneral::createFromVD(vec, nrows, false, changeOrder);

  if (optionPositive) _eigenVectors->makePositiveColumn();

  _flagEigenDecompose = true;

  return 0;
}

int AMatrixDense::_computeGeneralizedEigen(const MatrixSquareSymmetric &b, bool optionPositive)
{
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(*_eigenMatrix, *b._eigenMatrix);
  Eigen::VectorXd eigenValues  = solver.eigenvalues().real();
  Eigen::MatrixXd eigenVectors = solver.eigenvectors().real();

  return _terminateEigen(eigenValues, eigenVectors, optionPositive, true);
}

int AMatrixDense::_computeEigen(bool optionPositive)
{
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(*_eigenMatrix);
  Eigen::VectorXd eigenValues  = solver.eigenvalues().real();
  Eigen::MatrixXd eigenVectors = solver.eigenvectors().real();

  return _terminateEigen(eigenValues, eigenVectors, optionPositive, true);
}

bool AMatrixDense::_needToReset(int nrows, int ncols) 
{
  int newsize = nrows * ncols;
  
  if (newsize > _maxSize)
  {
    _maxSize = newsize;
    return true;
  }
  return false;
}

void AMatrixDense::_recopy(const AMatrixDense &r)
{
  *_eigenMatrix = *r._eigenMatrix;
  _flagEigenDecompose = r._flagEigenDecompose;
  _eigenValues = r._eigenValues;
  delete _eigenVectors;
  _eigenVectors = nullptr;
  if (r._eigenVectors != nullptr)
  {
    _eigenVectors = r._eigenVectors->clone();
  }
}

