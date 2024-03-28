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
  _recopyLocal(r);
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
    _recopyLocal(r);
  }
  return *this;
}

AMatrixDense::~AMatrixDense()
{
  delete _eigenVectors;
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
    _allocateLocal();
  }
  else
  {
    // Nothing to be done as the child will perform necessary duty
  }
}

void AMatrixDense::_deallocate()
{
  if (isFlagEigen())
  {
    // Possible code for Eigen should be placed here
  }
  else
  {
    // Nothing to be done as the child will perform necessary duty
  }
}

double AMatrixDense::_getValue(int irow, int icol) const
{
  if (isFlagEigen())
    return _getValueLocal(irow, icol);
  else
    my_throw("_getValue should never be called here");
  return TEST;
}

double AMatrixDense::_getValueByRank(int irank) const
{
  if (isFlagEigen())
    return _getValueLocal(irank);
  else
    my_throw("_getValue should never be called here");
  return TEST;
}

void AMatrixDense::_setValueByRank(int irank, double value)
{
  if (isFlagEigen())
    _setValueLocal(irank, value);
  else
    my_throw("_setValue should never be called here");
}

void AMatrixDense::_setValue(int irow, int icol, double value)
{
  if (isFlagEigen())
    _setValueLocal(irow, icol, value);
  else
    my_throw("_setValue should never be called here");
}

void AMatrixDense::_updValue(int irow, int icol, const EOperator& oper, double value)
{
  if (isFlagEigen())
    _updValueLocal(irow, icol, oper, value);
  else
    my_throw("_setValue should never be called here");
}

double& AMatrixDense::_getValueRef(int irow, int icol)
{
  if (isFlagEigen())
    return _getValueRefLocal(irow, icol);
  else
    my_throw("_getValueRef should never be called here");
  return AMatrix::_getValueRef(irow, icol);
}

int AMatrixDense::_getMatrixPhysicalSize() const
{
  if (isFlagEigen())
    return _getMatrixPhysicalSizeLocal();
  else
    my_throw("_getMatrixPhysicalSize should never be called here");
  return ITEST;
}

// Default storage in Eigen is column-major (see https://eigen.tuxfamily.org/dox/group__TopicStorageOrders.html)
int AMatrixDense::_getIndexToRank(int irow, int icol) const
{
  if (isFlagEigen())
    return _getIndexToRankLocal(irow, icol);
  else
    my_throw("_getIndexToRank should never be called here");
  return ITEST;
}

void AMatrixDense::_transposeInPlace()
{
  if (isFlagEigen())
    _transposeInPlaceLocal();
  else
    my_throw("_transposeInPlace should never be called here");
}

void AMatrixDense::_prodMatVecInPlacePtr(const double *x, double *y, bool transpose) const
{
  if (isFlagEigen())
    _prodMatVecInPlacePtrLocal(x, y, transpose);
  else
    my_throw("_prodMatVec should never be called here");
}

void AMatrixDense::_prodVecMatInPlacePtr(const double *x,double *y, bool transpose) const
{
  if (isFlagEigen())
    _prodVecMatInPlacePtrLocal(x, y, transpose);
  else
    my_throw("_prodVecMat should never be called here");
}

int AMatrixDense::_invert()
{
  if (isFlagEigen())
    return _invertLocal();
  else
    my_throw("_invert should never be called here");
  return 1;
}

int AMatrixDense::_solve(const VectorDouble &b, VectorDouble &x) const
{
  if (isFlagEigen())
    return _solveLocal(b, x);
  else
    my_throw("_solve should never be called here");
  return 1;
}

void AMatrixDense::setColumn(int icol, const VectorDouble& tab)
{
  if (isFlagEigen())
    _setColumnLocal(icol, tab);
  else
    AMatrix::setColumn(icol, tab);
}

void AMatrixDense::setRow(int irow, const VectorDouble& tab)
{
  if (isFlagEigen())
    _setRowLocal(irow, tab);
  else
    AMatrix::setRow(irow, tab);
}

void AMatrixDense::setDiagonal(const VectorDouble& tab)
{
  if (isFlagEigen())
    _setDiagonalLocal(tab);
  else
    AMatrix::setDiagonal(tab);
}

void AMatrixDense::setDiagonalToConstant(double value)
{
  if (isFlagEigen())
    _setDiagonalToConstantLocal(value);
  else
    AMatrix::setDiagonalToConstant(value);
}

void AMatrixDense::addScalar(double v)
{
  if (isFlagEigen())
    _addScalarLocal(v);
  else
    AMatrix::addScalar(v);
}

void AMatrixDense::addScalarDiag(double v)
{
  if (isFlagEigen())
    _addScalarDiagLocal(v);
  else
    AMatrix::addScalarDiag(v);
}

void AMatrixDense::prodScalar(double v)
{
  if (isFlagEigen())
    _prodScalarLocal(v);
  else
    AMatrix::prodScalar(v);
}

void AMatrixDense::addMatInPlace(const AMatrixDense& y, double cx, double cy)
{
  if (isFlagEigen() && y.isFlagEigen())
    _addMatInPlaceLocal(y, cx, cy);
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
    _prodMatMatInPlaceLocal(xm, ym, transposeX, transposeY);
  }
  else
    AMatrix::prodMatMatInPlace(x, y, transposeX, transposeY);
}

void AMatrixDense::prodNormMatMatInPlace(const AMatrixDense &a,
                                         const AMatrixDense &m,
                                         bool transpose)
{
  if (isFlagEigen() && a.isFlagEigen() && m.isFlagEigen())
    _prodNormMatMatInPlaceLocal(a, m, transpose);
  else
    AMatrix::prodNormMatMatInPlace(a, m, transpose);
}

/*! Product 't(A)' %*% ['vec'] %*% 'A' or 'A' %*% ['vec'] %*% 't(A)' stored in 'this'*/
void AMatrixDense::prodNormMatInPlace(const AMatrixDense &a, const VectorDouble& vec, bool transpose)
{
  if (isFlagEigen() && a.isFlagEigen())
    _prodNormMatInPlaceLocal(a, vec, transpose);
  else
    AMatrix::prodNormMatInPlace(a, vec, transpose);
}

void AMatrixDense::fill(double value)
{
  if (isFlagEigen())
    _fillLocal(value);
  else
    AMatrix::fill(value);
}

/*! Multiply a Matrix row-wise */
void AMatrixDense::multiplyRow(const VectorDouble& vec)
{
  if (isFlagEigen())
  {
    _multiplyRowLocal(vec);
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
    _multiplyColumnLocal(vec);
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
    _divideRowLocal(vec);
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
    _divideColumnLocal(vec);
  }
  else
  {
    AMatrix::divideColumn(vec);
  }
}

void AMatrixDense::_prodMatVecInPlacePtrLocal(const double *x, double *y, bool transpose) const
{
  Eigen::Map<const Eigen::VectorXd> xm(x, getNCols());
  Eigen::Map<Eigen::VectorXd> ym(y, getNRows());
  if (transpose)
    ym.noalias() = _eigenMatrix.transpose() * xm;
  else
    ym.noalias() = _eigenMatrix * xm;
}

void AMatrixDense::_prodVecMatInPlacePtrLocal(const double *x, double *y, bool transpose) const
{
  Eigen::Map<const Eigen::VectorXd> xm(x, getNRows());
  Eigen::Map<Eigen::VectorXd> ym(y, getNCols());
  if (transpose)
    ym.noalias() = xm.transpose() * _eigenMatrix.transpose();
  else
    ym.noalias() = xm.transpose() * _eigenMatrix;
}

/*! Perform 'vec' * 'this' */
VectorDouble AMatrixDense::prodVecMat(const VectorDouble& x, bool transpose) const
{
  if (isFlagEigen())
    return _prodVecMatLocal(x, transpose);
  else
    return AMatrix::prodVecMat(x, transpose);
}

/*! Perform 'this' * 'vec' */
VectorDouble AMatrixDense::prodMatVec(const VectorDouble& x, bool transpose) const
{
  if (isFlagEigen())
    return _prodMatVecLocal(x, transpose);
  else
    return AMatrix::prodMatVec(x, transpose);
}

/*! Extract a Row */
VectorDouble AMatrixDense::getRow(int irow) const
{
  if (isFlagEigen())
  {
    return _getRowLocal(irow);
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
    return _getColumnLocal(icol);
  }
  else
  {
    return AMatrix::getColumn(icol);
  }
}

int AMatrixDense::_computeEigen(bool optionPositive)
{
  if (!isSquare())
  {
    messerr("The current Matrix does not seems to be square");
    return 1;
  }

  if (isFlagEigen())
    return _computeEigenLocal(optionPositive);
  else
    my_throw("'_computeEigen' should never be called here");
  return ITEST;
}

int AMatrixDense::_computeGeneralizedEigen(const MatrixSquareSymmetric& b, bool optionPositive)
{
  if (!isSquare())
  {
    messerr("The current Matrix does not seems to be square");
    return 1;
  }

  if (isFlagEigen())
    return _computeGeneralizedEigenLocal(b, optionPositive);
  else
    my_throw("'_computeGeneralizedEigen' should never be called here");
  return ITEST;
}

/// =========================================================================
/// The subsequent methods rely on the specific local storage ('eigenMatrix')
/// =========================================================================
int AMatrixDense::_solveLocal(const VectorDouble &b, VectorDouble &x) const
{
  if (!isSquare())
  {
    my_throw("Invert method is limited to Square Matrices");
    return 1;
  }
  Eigen::Map<const Eigen::VectorXd> bm(b.data(), getNCols());
  Eigen::Map<Eigen::VectorXd> xm(x.data(), getNRows());
  xm = _eigenMatrix.inverse() * bm;
  return 0;
}

int AMatrixDense::_invertLocal()
{
  if (!isSquare())
  {
    my_throw("Invert method is limited to Square matrices");
    return 1;
  }
  _eigenMatrix = _eigenMatrix.inverse();
  return 0;
}

void AMatrixDense::_transposeInPlaceLocal()
{
  _eigenMatrix.transposeInPlace();
}

// Default storage in Eigen is column-major (see https://eigen.tuxfamily.org/dox/group__TopicStorageOrders.html)
int AMatrixDense::_getIndexToRankLocal(int irow, int icol) const
{
  return (icol * getNRows() + irow);
}

int AMatrixDense::_getMatrixPhysicalSizeLocal() const
{
  return _eigenMatrix.size();
}

double& AMatrixDense::_getValueRefLocal(int irow, int icol)
{
  return *(_eigenMatrix.data() + _getIndexToRank(irow, icol));
}

void AMatrixDense::_setValueLocal(int irow, int icol, double value)
{
  _eigenMatrix(irow, icol) = value;
}

void AMatrixDense::_updValueLocal(int irow, int icol, const EOperator& oper, double value)
{
  _eigenMatrix(irow, icol) = modifyOperator(oper, _eigenMatrix(irow, icol), value);
}

void AMatrixDense::_setValueLocal(int irank, double value)
{
  *(_eigenMatrix.data() + irank) = value;
}

double AMatrixDense::_getValueLocal(int irank) const
{
  return *(_eigenMatrix.data() + irank);
}

double AMatrixDense::_getValueLocal(int irow, int icol) const
{
  return _eigenMatrix(irow, icol);
}

void AMatrixDense::_allocateLocal()
{
  if (isMultiThread()) omp_set_num_threads(getMultiThread()); // TODO Move to multithread handling class
  _eigenMatrix = Eigen::MatrixXd::Constant(getNRows(),getNCols(),0.);
}

void AMatrixDense::_recopyLocal(const AMatrixDense &r)
{
  _eigenMatrix = r._eigenMatrix;
  _flagEigenDecompose = r._flagEigenDecompose;
  _eigenValues = r._eigenValues;
  if (_eigenVectors != nullptr) _eigenVectors = r._eigenVectors->clone();
}

void AMatrixDense::_setColumnLocal(int icol, const VectorDouble& tab)
{
  Eigen::Map<const Eigen::VectorXd> tabm(tab.data(), getNRows());
  _eigenMatrix.col(icol) = tabm;
}

void AMatrixDense::_setRowLocal(int irow, const VectorDouble& tab)
{
  Eigen::Map<const Eigen::VectorXd> tabm(tab.data(), getNCols());
  _eigenMatrix.row(irow) = tabm;
}

void AMatrixDense::_setDiagonalLocal(const VectorDouble& tab)
{
  _eigenMatrix.setZero();
  Eigen::Map<const Eigen::VectorXd> tabm(tab.data(), getNRows());
  _eigenMatrix.diagonal() = tabm;
}

void AMatrixDense::_setDiagonalToConstantLocal(double value)
{
  _eigenMatrix.setZero();
  _eigenMatrix.diagonal() = Eigen::VectorXd::Constant(getNRows(),value);
}

void AMatrixDense::_addScalarLocal(double v)
{
  _eigenMatrix.array() += v;
}

void AMatrixDense::_addScalarDiagLocal(double v)
{
  _eigenMatrix.diagonal() += Eigen::VectorXd::Constant(getNRows(),v);
}

void AMatrixDense::_prodScalarLocal(double v)
{
  _eigenMatrix.array() *= v;
}

void AMatrixDense::_addMatInPlaceLocal(const AMatrixDense& y, double cx, double cy)
{
  _eigenMatrix.noalias() = cx * _eigenMatrix + cy * y._eigenMatrix;
}

void AMatrixDense::_prodMatMatInPlaceLocal(const AMatrixDense *x,
                                           const AMatrixDense *y,
                                           bool transposeX,
                                           bool transposeY)
{
  if (transposeX)
  {
    if (transposeY)
    {
      _eigenMatrix.noalias() = x->_eigenMatrix.transpose() * y->_eigenMatrix.transpose();
    }
    else
    {
      _eigenMatrix.noalias() = x->_eigenMatrix.transpose() * y->_eigenMatrix;
    }
  }
  else
  {
    if (transposeY)
    {
      _eigenMatrix.noalias() = x->_eigenMatrix * y->_eigenMatrix.transpose();
    }
    else
    {
      _eigenMatrix.noalias() = x->_eigenMatrix * y->_eigenMatrix;
    }
  }
}

void AMatrixDense::_prodNormMatMatInPlaceLocal(const AMatrixDense& a, const AMatrixDense& m, bool transpose)
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

void AMatrixDense::_prodNormMatInPlaceLocal(const AMatrixDense& a, const VectorDouble& vec, bool transpose)
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

void AMatrixDense::_fillLocal(double value)
{
  _eigenMatrix.setConstant(value);
}

void AMatrixDense::_multiplyRowLocal(const VectorDouble& vec)
{
  Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), getNCols());
  _eigenMatrix = vecm.asDiagonal() * _eigenMatrix;
}

void AMatrixDense::_multiplyColumnLocal(const VectorDouble& vec)
{
  Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), getNRows());
  _eigenMatrix = _eigenMatrix * vecm.asDiagonal();
}

void AMatrixDense::_divideRowLocal(const VectorDouble& vec)
{
  VectorDouble temp = VH::inverse(vec);
  Eigen::Map<const Eigen::VectorXd> vecm(temp.data(), getNCols());
  _eigenMatrix = vecm.asDiagonal() * _eigenMatrix;
}

void AMatrixDense::_divideColumnLocal(const VectorDouble& vec)
{
  VectorDouble temp = VH::inverse(vec);
  Eigen::Map<const Eigen::VectorXd> vecm(temp.data(), getNRows());
  _eigenMatrix = _eigenMatrix * vecm.asDiagonal();
}

VectorDouble AMatrixDense::_prodMatVecLocal(const VectorDouble& x, bool transpose) const
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

VectorDouble AMatrixDense::_prodVecMatLocal(const VectorDouble& x, bool transpose) const
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

/*! Extract a Row */
VectorDouble AMatrixDense::_getRowLocal(int irow) const
{
  Eigen::VectorXd resm = _eigenMatrix.row(irow);
  VectorDouble res(resm.data(), resm.data() + resm.size());
  return res;
}

/*! Extract a Column */
VectorDouble AMatrixDense::_getColumnLocal(int icol) const
{
  Eigen::VectorXd resm = _eigenMatrix.col(icol);
  VectorDouble res(resm.data(), resm.data() + resm.size());
  return res;
}

int AMatrixDense::_computeEigenLocal(bool optionPositive)
{
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(_eigenMatrix);

  _flagEigenDecompose = true;

  Eigen::VectorXd eigenValues  = solver.eigenvalues().real();
  Eigen::MatrixXd eigenVectors = solver.eigenvectors().real();

  int nrows = getNRows();
  int ncols = getNCols();

  _eigenValues = VectorDouble(nrows);
  Eigen::Map<Eigen::VectorXd>(&_eigenValues[0], eigenValues.size()) = eigenValues;
  std::reverse(_eigenValues.begin(), _eigenValues.end());

  VectorDouble vec(nrows * ncols);
  Eigen::Map<Eigen::MatrixXd>(&vec[0], nrows, ncols) = eigenVectors;
  _eigenVectors = MatrixSquareGeneral::createFromVD(vec, nrows, false, 1, true);

  if (optionPositive) _eigenVectors->makePositiveColumn();

  return 0;
}

int AMatrixDense::_computeGeneralizedEigenLocal(const MatrixSquareSymmetric &b, bool optionPositive)
{
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(_eigenMatrix, b._eigenMatrix);

  _flagEigenDecompose = true;

  Eigen::VectorXd eigenValues  = solver.eigenvalues().real();
  Eigen::MatrixXd eigenVectors = solver.eigenvectors().real();

  int nrows = getNRows();
  int ncols = getNCols();

  _eigenValues = VectorDouble(nrows);
  Eigen::Map<Eigen::VectorXd>(&_eigenValues[0], eigenValues.size()) = eigenValues;
  std::reverse(_eigenValues.begin(), _eigenValues.end());

  VectorDouble vec(nrows * ncols);
  Eigen::Map<Eigen::MatrixXd>(&vec[0], nrows, ncols) = eigenVectors;
  _eigenVectors = MatrixSquareGeneral::createFromVD(vec, nrows, false, 1, true);

  if (optionPositive) _eigenVectors->makePositiveColumn();

  return 0;
}

