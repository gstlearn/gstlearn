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
  _allocate();
  _recopy(r);
}

AMatrixDense::AMatrixDense(const AMatrix &r)
  : AMatrix(r),
    _flagEigenDecompose(false),
    _eigenValues(),
    _eigenVectors(nullptr),
    _eigenMatrix()
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
  _eigenMatrix.resize(getNRows() * getNCols());
}

void AMatrixDense::_deallocate()
{
  if (_eigenVectors != nullptr)
  {
    delete _eigenVectors;
    _eigenVectors = nullptr;
  }
  _eigenValues.clear();
  _eigenMatrix.clear();
}

double AMatrixDense::getValue(int irow, int icol, bool flagCheck) const
{
  if (flagCheck && ! _isIndexValid(irow, icol)) return TEST;
  return getEigenMat()(irow, icol);
}

double AMatrixDense::_getValueByRank(int irank) const
{
  return *(getEigenMat().data() + irank);
}

constvect AMatrixDense::getColumnPtr(int icol) const
{
  const auto a = getEigenMat().col(icol);
  size_t n = getNRows();
  return {a.data(), n};
}
void AMatrixDense::_setValueByRank(int irank, double value)
{
  *(getEigenMat().data() + irank) = value;
}

void AMatrixDense::setValue(int irow, int icol, double value, bool flagCheck)
{
  if (flagCheck && ! _isIndexValid(irow, icol)) return;
  getEigenMat()(irow, icol) = value;
  if (mustBeSymmetric() && irow != icol) getEigenMat()(icol, irow) = value;
}

void AMatrixDense::updValue(int irow,
                            int icol,
                            const EOperator &oper,
                            double value,
                            bool flagCheck)
{
  if (flagCheck && ! _isIndexValid(irow, icol)) return;
  double result = modifyOperator(oper, getEigenMat()(irow, icol), value);
  getEigenMat()(irow, icol) = result;
  if (mustBeSymmetric() && irow != icol)
    getEigenMat()(icol, irow) = result;
}

double& AMatrixDense::_getValueRef(int irow, int icol)
{
  return *(getEigenMat().data() + _getIndexToRank(irow, icol));
}

int AMatrixDense::_getMatrixPhysicalSize() const
{
  return getEigenMat().size();
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
  getEigenMat().transposeInPlace();
  _setNCols(nrows);
  _setNRows(ncols);
}

void AMatrixDense::_addProdMatVecInPlaceToDestPtr(const double *x,double *y, bool transpose) const
{
  Eigen::Map<const Eigen::VectorXd> xm(x, getNCols());
  Eigen::Map<Eigen::VectorXd> ym(y, getNRows());
  if (transpose)
    ym.noalias() += getEigenMat().transpose() * xm;
  else
    ym.noalias() += getEigenMat() * xm;

}

//TODO supress this method and implement it in the virtual class AMatrix
void AMatrixDense::_prodMatVecInPlacePtr(const double *x, double *y, bool transpose) const
{
  int nc = transpose? getNRows() : getNCols();
  int nr = transpose? getNCols() : getNRows();
  Eigen::Map<const Eigen::VectorXd> xm(x, nc);
  Eigen::Map<Eigen::VectorXd> ym(y, nr);
  if (transpose)
    ym.noalias() = getEigenMat().transpose() * xm;
  else
    ym.noalias() = getEigenMat() * xm;
}

void AMatrixDense::_prodVecMatInPlacePtr(const double *x,double *y, bool transpose) const
{
  Eigen::Map<const Eigen::VectorXd> xm(x, getNRows());
  Eigen::Map<Eigen::VectorXd> ym(y, getNCols());
  if (transpose)
    ym.noalias() = xm.transpose() * getEigenMat().transpose();
  else
    ym.noalias() = xm.transpose() * getEigenMat();
}

int AMatrixDense::_invert()
{
  /// TODO : check beforehand if matrix is invertible ?
  getEigenMat() = getEigenMat().inverse();
  return 0;
}

int AMatrixDense::invert2(AMatrixDense& res) const
{
  /// TODO : check beforehand if matrix is invertible ?
  res.getEigenMat() = getEigenMat().inverse();
  return 0;
}
int AMatrixDense::_solve(const VectorDouble &b, VectorDouble &x) const
{
  /// TODO : check beforehand if matrix is invertible ?
  Eigen::Map<const Eigen::VectorXd> bm(b.data(), getNCols());
  Eigen::Map<Eigen::VectorXd> xm(x.data(), getNRows());
  xm = getEigenMat().inverse() * bm;
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
  getEigenMat().col(icol) = tabm;
}

void AMatrixDense::setRow(int irow, const VectorDouble& tab, bool flagCheck)
{
  if (flagCheck)
  {
    if (! _isRowValid(irow)) return;
    if (! _isRowSizeConsistent(tab)) return;
  }
  Eigen::Map<const Eigen::VectorXd> tabm(tab.data(), getNCols());
  getEigenMat().row(irow) = tabm;
}

void AMatrixDense::setDiagonal(const VectorDouble& tab, bool flagCheck)
{
  if (flagCheck)
  {
    if (! _isRowSizeConsistent(tab)) return;
  }
  getEigenMat().setZero();
  Eigen::Map<const Eigen::VectorXd> tabm(tab.data(), getNRows());
  getEigenMat().diagonal() = tabm;
}

void AMatrixDense::setDiagonalToConstant(double value)
{
  getEigenMat().setZero();
  getEigenMat().diagonal() = Eigen::VectorXd::Constant(getNRows(),value);
}

void AMatrixDense::addScalar(double v)
{
  getEigenMat().array() += v;
}

void AMatrixDense::addScalarDiag(double v)
{
  getEigenMat().diagonal() += Eigen::VectorXd::Constant(getNRows(),v);
}

void AMatrixDense::prodScalar(double v)
{
  getEigenMat().array() *= v;
}

void AMatrixDense::addMatInPlace(const AMatrixDense& y, double cx, double cy)
{
  getEigenMat().noalias() = cx * getEigenMat() + cy * y.getEigenMat();
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
        getEigenMat().noalias() = xm->getEigenMat().transpose() * ym->getEigenMat().transpose();
      }
      else
      {
        getEigenMat().noalias() = xm->getEigenMat().transpose() * ym->getEigenMat();
      }
    }
    else
    {
      if (transposeY)
      {
        getEigenMat().noalias() = xm->getEigenMat() * ym->getEigenMat().transpose();
      }
      else
      {
        getEigenMat().noalias() = xm->getEigenMat() * ym->getEigenMat();
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
    getEigenMat().noalias() = a->getEigenMat().transpose() * m->getEigenMat() * a->getEigenMat();
  }
  else
  {
    getEigenMat().noalias() = a->getEigenMat() * m->getEigenMat() * a->getEigenMat().transpose();
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
      getEigenMat().noalias() = a.getEigenMat().transpose() * a.getEigenMat();
    else
    {
      Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), vec.size());
      getEigenMat().noalias() = a.getEigenMat().transpose() * vecm * a.getEigenMat();
    }
  }
  else
  {
    if (vec.empty())
      getEigenMat().noalias() = a.getEigenMat() * a.getEigenMat().transpose();
    else
    {
      Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), vec.size());
      getEigenMat().noalias() = a.getEigenMat() * vecm * a.getEigenMat().transpose();
    }
  }
}

void AMatrixDense::fill(double value)
{
  getEigenMat().setConstant(value);
}

/*! Multiply a Matrix row-wise */
void AMatrixDense::multiplyRow(const VectorDouble& vec)
{
  Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), getNCols());
  getEigenMat() = vecm.asDiagonal() * getEigenMat();
}

/*! Multiply a Matrix column-wise */
void AMatrixDense::multiplyColumn(const VectorDouble& vec)
{
  Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), getNRows());
  getEigenMat() = getEigenMat() * vecm.asDiagonal();
}

/*! Divide a Matrix row-wise */
void AMatrixDense::divideRow(const VectorDouble& vec)
{
  VectorDouble temp = VH::inverse(vec);
  Eigen::Map<const Eigen::VectorXd> vecm(temp.data(), getNCols());
 getEigenMat() = vecm.asDiagonal() * getEigenMat();
}

/*! Divide a Matrix column-wise */
void AMatrixDense::divideColumn(const VectorDouble& vec)
{
  VectorDouble temp = VH::inverse(vec);
  Eigen::Map<const Eigen::VectorXd> vecm(temp.data(), getNRows());
  getEigenMat() = getEigenMat() * vecm.asDiagonal();
}

/*! Perform 'vec' * 'this' */
VectorDouble AMatrixDense::prodVecMat(const VectorDouble& x, bool transpose) const
{
  Eigen::Map<const Eigen::VectorXd> xm(x.data(), x.size());
  VectorDouble y(transpose ? getNRows() : getNCols());
  Eigen::Map<Eigen::VectorXd> ym(y.data(), y.size());
  if (transpose)
    ym = xm.transpose() * getEigenMat().transpose();
  else
    ym = xm.transpose() * getEigenMat();
  return y;
}

/*! Perform 'this' * 'vec' */
VectorDouble AMatrixDense::prodMatVec(const VectorDouble& x, bool transpose) const
{
  Eigen::Map<const Eigen::VectorXd> xm(x.data(), x.size());
  VectorDouble y(transpose ? getNCols() : getNRows());
  Eigen::Map<Eigen::VectorXd> ym(y.data(), y.size());
  if (transpose)
    ym = getEigenMat().transpose() * xm;
  else
    ym = getEigenMat() * xm;
  return y;
}

/*! Extract a Row */
VectorDouble AMatrixDense::getRow(int irow) const
{
  VectorDouble res(getNCols());
  for (size_t i = 0; i < res.size(); ++i)
  {
    res[i] = getEigenMat().row(irow)[i];
  }
  return res;
}

/*! Extract a Column */
VectorDouble AMatrixDense::getColumn(int icol) const
{
  VectorDouble res(getNRows());
  for (size_t i = 0; i < res.size(); ++i)
  {
    res[i] = getEigenMat().col(icol)[i];
  }
  return res;
}

constvect AMatrixDense::getViewOnColumn(int icol) const
{
  constvect res(getEigenMat().col(icol).data(),getNRows());
  return res;
}

vect AMatrixDense::getViewOnColumnModify(int icol)
{
  vect res(getEigenMat().col(icol).data(),getNRows());
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
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(getEigenMat(), b.getEigenMat());
  Eigen::VectorXd eigenValues  = solver.eigenvalues().real();
  Eigen::MatrixXd eigenVectors = solver.eigenvectors().real();

  return _terminateEigen(eigenValues, eigenVectors, optionPositive, true);
}

int AMatrixDense::_computeEigen(bool optionPositive)
{
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(getEigenMat());
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

