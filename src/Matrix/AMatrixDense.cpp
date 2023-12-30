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
#include "Matrix/MatrixFactory.hpp"
#include "Matrix/AMatrix.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/AException.hpp"
#include <math.h>

AMatrixDense::AMatrixDense(int nrow, int ncol)
  : AMatrix(nrow, ncol),
    _eigenMatrix()
{
  _allocate();
}

AMatrixDense::AMatrixDense(const AMatrixDense &r)
  : AMatrix(r),
    _eigenMatrix()
{
  // operator "=" is faster than the copy constructor
  // (see https://stackoverflow.com/questions/47644021/eigen-copy-constructor-vs-operator-performance)
  _allocate();
  _recopyLocal(r);
}

AMatrixDense::AMatrixDense(const AMatrix &m)
    : AMatrix(m),
      _eigenMatrix()
{
  if (m.isEmpty())
  {
    messerr("The input matrix should be non-empty");
    _clear();
    return;
  }
  _allocate();
  if (isFlagEigen()) copyElements(m);
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

void AMatrixDense::_prodVectorInPlace(const double *inv, double *outv) const
{
  if (isFlagEigen())
    _prodVectorLocal(inv, outv);
  else
    my_throw("_prodVector should never be called here");
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

void AMatrixDense::addMatrix(const AMatrixDense& y, double value)
{
  if (isFlagEigen())
    _addMatrixLocal(y);
  else
    AMatrix::addMatrix(y, value);
}

void AMatrixDense::prodMatrix(const AMatrixDense& x, const AMatrixDense& y)
{
  if (isFlagEigen())
    _prodMatrixLocal(x, y);
  else
    AMatrix::prodMatrix(x, y);
}

void AMatrixDense::prodTMatrix(const AMatrixDense& x, const AMatrixDense& y)
{
  if (isFlagEigen())
    _prodTMatrixLocal(x, y);
  else
    AMatrix::prodTMatrix(x, y);
}

void AMatrixDense::linearCombination(double cx, double cy, const AMatrixDense& y)
{
  if (isFlagEigen())
    _linearCombinationLocal(cx, cy, y);
  else
    AMatrix::linearCombination(cx, cy, y);
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

/*! Perform M * 'vec' */
VectorDouble AMatrixDense::prodVector(const VectorDouble& vec) const
{
  if (isFlagEigen())
  {
    return _prodVectorLocal(vec);
  }
  else
  {
    return AMatrix::prodVector(vec);
  }
}

/*! Perform 'vec'^T * M */
VectorDouble AMatrixDense::prodTVector(const VectorDouble& vec) const
{
  if (isFlagEigen())
  {
    return _prodTVectorLocal(vec);
  }
  else
  {
    return AMatrix::prodTVector(vec);
  }
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

/// =========================================================================
/// The subsequent methods rely on the specific local storage ('eigenMatrix')
/// =========================================================================
int AMatrixDense::_solveLocal(const VectorDouble &b, VectorDouble &x) const
{
  if (!isSquare() || !isSymmetric())
  {
    my_throw("Invert method is limited to Square Symmetrical Matrices");
    return 1;
  }
  Eigen::Map<const Eigen::VectorXd> bm(b.data(), getNCols());
  Eigen::Map <Eigen::VectorXd> xm(x.data(), getNRows());
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

void AMatrixDense::_prodVectorLocal(const double *inv, double *outv) const
{
  Eigen::Map<const Eigen::VectorXd> inm(inv, getNCols());
  Eigen::Map<Eigen::VectorXd> outm(outv, getNRows());
  outm.noalias() = _eigenMatrix * inm;
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

void AMatrixDense::_addMatrixLocal(const AMatrixDense& y, double value)
{
  _eigenMatrix.noalias() += y._eigenMatrix * value;
}

void AMatrixDense::_prodMatrixLocal(const AMatrixDense& x, const AMatrixDense& y)
{
  _eigenMatrix.noalias() = x._eigenMatrix * y._eigenMatrix;
}

void AMatrixDense::_prodTMatrixLocal(const AMatrixDense& x, const AMatrixDense& y)
{
  _eigenMatrix.noalias() = x._eigenMatrix.transpose() * y._eigenMatrix;
}

/**
 * Perform: 'this' = cx * 'this' + cy * y (cannot be noalias'ed)
 * @param cx
 * @param cy
 * @param y
 */
void AMatrixDense::_linearCombinationLocal(double cx, double cy,const AMatrixDense &y)
{
  _eigenMatrix = cx * _eigenMatrix + cy * y._eigenMatrix;
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

VectorDouble AMatrixDense::_prodVectorLocal(const VectorDouble& vec) const
{
  Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), getNCols());
  Eigen::VectorXd resm = _eigenMatrix * vecm;
  VectorDouble res(resm.data(), resm.data() + resm.size());
  return res;
}

VectorDouble AMatrixDense::_prodTVectorLocal(const VectorDouble& vec) const
{
  Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), getNRows());
  Eigen::VectorXd resm = vecm.transpose() * _eigenMatrix;
  VectorDouble res(resm.data(), resm.data() + resm.size());
  return res;
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
