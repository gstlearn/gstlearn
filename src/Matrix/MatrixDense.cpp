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
#include "Matrix/MatrixDense.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Matrix/AMatrix.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Matrix/MatrixSquare.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "geoslib_define.h"

#include <cmath>
#include <omp.h>

namespace gstlrn
{
MatrixDense::MatrixDense(Id nrow, Id ncol)
  : AMatrix(nrow, ncol)
  , _eigenValues()
  , _eigenVectors(nullptr)
  , _eigenMatrix()
  , _maxSize(nrow * ncol)
{
  _allocate();
}

MatrixDense::MatrixDense(const MatrixDense& r)
  : AMatrix(r)
  , _eigenValues()
  , _eigenVectors(nullptr)
  , _eigenMatrix()
  , _maxSize(r._maxSize)

{
  _allocate();
  _recopy(r);
}

MatrixDense::MatrixDense(const AMatrix& r)
  : AMatrix(r)
  , _eigenValues()
  , _eigenVectors(nullptr)
  , _eigenMatrix()
  , _maxSize(r.getNRows() * r.getNCols())
{
  _allocate();
  copyElements(r);
}

MatrixDense& MatrixDense::operator=(const MatrixDense& r)
{
  if (this != &r)
  {
    AMatrix::operator=(r);
    _maxSize = r._maxSize;
    _allocate();
    _recopy(r);
  }
  return *this;
}

MatrixDense::~MatrixDense()
{
  _deallocate();
}

void MatrixDense::_allocate()
{
  _maxSize = getNRows() * getNCols();
  _deallocate();
  _eigenMatrix.resize(getNRows() * getNCols());
}

void MatrixDense::_deallocate()
{
  if (_eigenVectors != nullptr)
  {
    delete _eigenVectors;
    _eigenVectors = nullptr;
  }
  _eigenValues.clear();
  _eigenMatrix.clear();
}

double MatrixDense::getValue(Id irow, Id icol) const
{
  if (getFlagMatrixCheck() && !_isIndexValid(irow, icol)) return TEST;
  return eigenMat()(irow, icol);
}

double MatrixDense::_getValueByRank(Id irank) const
{
  return *(eigenMat().data() + irank);
}

constvect MatrixDense::getColumnPtr(Id icol) const
{
  const auto a = eigenMat().col(icol);
  size_t n     = getNRows();
  return {a.data(), n};
}
void MatrixDense::_setValueByRank(Id irank, double value)
{
  *(eigenMat().data() + irank) = value;
}

void MatrixDense::setValue(Id irow, Id icol, double value)
{
  if (getFlagMatrixCheck() && !_isIndexValid(irow, icol)) return;
  eigenMat().coeffRef(irow, icol) = value;
  if (mustBeSymmetric() && irow != icol)
    eigenMat().coeffRef(icol, irow) = value;
}

double MatrixDense::traceProd(const MatrixDense& a, MatrixDense& b)
{
  if (a.getNRows() != b.getNRows() || a.getNCols() != b.getNCols())
  {
    messerr("MatrixDense::traceProd: incompatible matrix sizes");
    return TEST;
  }

  b.eigenMat().array() *= a.eigenMat().transpose().array();
  return b.eigenMat().sum();
}
void MatrixDense::updValue(Id irow,
                           Id icol,
                           const EOperator& oper,
                           double value)
{
  if (getFlagMatrixCheck() && !_isIndexValid(irow, icol)) return;
  double result                   = modifyOperator(oper, eigenMat()(irow, icol), value);
  eigenMat().coeffRef(irow, icol) = result;
  eigenMat()(irow, icol)          = result;
  if (mustBeSymmetric() && irow != icol)
    eigenMat()(icol, irow) = result;
}

double& MatrixDense::_getValueRef(Id irow, Id icol)
{
  return *(eigenMat().data() + _getIndexToRank(irow, icol));
}

Id MatrixDense::_getMatrixPhysicalSize() const
{
  return eigenMat().size();
}

// Default storage in Eigen is column-major (see https://eigen.tuxfamily.org/dox/group__TopicStorageOrders.html)
Id MatrixDense::_getIndexToRank(Id irow, Id icol) const
{
  return (icol * getNRows() + irow);
}

void MatrixDense::_transposeInPlace()
{
  auto nrows = getNRows();
  auto ncols = getNCols();
  eigenMat().transposeInPlace();
  _setNCols(nrows);
  _setNRows(ncols);
}

void MatrixDense::_addProdMatVecInPlacePtr(constvect x, vect y, bool transpose) const
{
  if (transpose)
  {
    Eigen::Map<const Eigen::VectorXd> xm(x.data(), getNRows());
    Eigen::Map<Eigen::VectorXd> ym(y.data(), getNCols());
    ym.noalias() += eigenMat().transpose() * xm;
  }
  else
  {
    Eigen::Map<const Eigen::VectorXd> xm(x.data(), getNCols());
    Eigen::Map<Eigen::VectorXd> ym(y.data(), getNRows());
    ym.noalias() += eigenMat() * xm;
  }
}

void MatrixDense::_addProdVecMatInPlacePtr(constvect x, vect y, bool transpose) const
{
  if (transpose)
  {
    Eigen::Map<const Eigen::VectorXd> xm(x.data(), getNCols());
    Eigen::Map<Eigen::VectorXd> ym(y.data(), getNRows());
    ym.noalias() += xm.transpose() * eigenMat().transpose();
  }
  else
  {
    Eigen::Map<const Eigen::VectorXd> xm(x.data(), getNRows());
    Eigen::Map<Eigen::VectorXd> ym(y.data(), getNCols());
    ym.noalias() += xm.transpose() * eigenMat();
  }
}

Id MatrixDense::_invert()
{
  /// TODO : check beforehand if matrix is invertible ?
  eigenMat() = eigenMat().inverse();
  return 0;
}

Id MatrixDense::invert2(MatrixDense& res) const
{
  /// TODO : check beforehand if matrix is invertible ?
  res.eigenMat() = eigenMat().inverse();
  return 0;
}
Id MatrixDense::_solve(const VectorDouble& b, VectorDouble& x) const
{
  /// TODO : check beforehand if matrix is invertible ?
  Eigen::Map<const Eigen::VectorXd> bm(b.data(), getNCols());
  Eigen::Map<Eigen::VectorXd> xm(x.data(), getNRows());
  xm = eigenMat().inverse() * bm;
  return 0;
}

void MatrixDense::setColumn(Id icol, const VectorDouble& tab)
{
  if (getFlagMatrixCheck())
  {
    if (!_isColumnValid(icol)) return;
    if (!_isColumnSizeConsistent(tab)) return;
  }
  Eigen::Map<const Eigen::VectorXd> tabm(tab.data(), getNRows());
  eigenMat().col(icol) = tabm;
}

void MatrixDense::setColumnToConstant(Id icol, double value)
{
  if (getFlagMatrixCheck())
  {
    if (!_isColumnValid(icol)) return;
  }
  eigenMat().col(icol) = Eigen::VectorXd::Constant(getNRows(), value);
}

void MatrixDense::setRow(Id irow, const VectorDouble& tab)
{
  if (getFlagMatrixCheck())
  {
    if (!_isRowValid(irow)) return;
    if (!_isRowSizeConsistent(tab)) return;
  }
  Eigen::Map<const Eigen::VectorXd> tabm(tab.data(), getNCols());
  eigenMat().row(irow) = tabm;
}

void MatrixDense::setRowToConstant(Id irow, double value)
{
  if (getFlagMatrixCheck())
  {
    if (!_isRowValid(irow)) return;
  }
  eigenMat().row(irow) = Eigen::VectorXd::Constant(getNCols(), value);
}

void MatrixDense::setDiagonal(const VectorDouble& tab)
{
  if (getFlagMatrixCheck())
  {
    if (!_isRowSizeConsistent(tab)) return;
  }
  eigenMat().setZero();
  Eigen::Map<const Eigen::VectorXd> tabm(tab.data(), getNRows());
  eigenMat().diagonal() = tabm;
}

void MatrixDense::setDiagonalToConstant(double value)
{
  eigenMat().setZero();
  eigenMat().diagonal() = Eigen::VectorXd::Constant(getNRows(), value);
}

void MatrixDense::addScalar(double v)
{
  eigenMat().array() += v;
}

void MatrixDense::addScalarDiag(double v)
{
  if (isZero(v)) return;
  eigenMat().diagonal() += Eigen::VectorXd::Constant(getNRows(), v);
}

void MatrixDense::prodScalar(double v)
{
  if (isOne(v)) return;
  eigenMat().array() *= v;
}

void MatrixDense::addMat(const AMatrix& y, double cx, double cy)
{
  const auto* ym = dynamic_cast<const MatrixDense*>(&y);
  if (ym == nullptr || ym == this)
  {
    AMatrix::addMat(y, cx, cy);
  }
  else
  {
    eigenMat().noalias() = cx * eigenMat();
    if (cy == 0. || (getFlagMatrixCheck() && !isSameSize(y))) return;
    eigenMat().noalias() += cy * ym->eigenMat();
  }
}

void MatrixDense::prodMatMatInPlace(const AMatrix* x,
                                    const AMatrix* y,
                                    bool transposeX,
                                    bool transposeY)
{
  if (getFlagMatrixCheck() &&
      !_isMatrixCompatible("MatrixDense::prodMatMatInPlace",
                           x, 0, transposeX,
                           y, 0, transposeY)) return;

  const auto* xm = dynamic_cast<const MatrixDense*>(x);
  const auto* ym = dynamic_cast<const MatrixDense*>(y);
  if (xm == nullptr || ym == nullptr)
  {
    AMatrix::prodMatMatInPlace(x, y, transposeX, transposeY);
  }
  else
  {
    if (transposeX)
    {
      if (transposeY)
      {
        eigenMat().noalias() = xm->eigenMat().transpose() * ym->eigenMat().transpose();
      }
      else
      {
        eigenMat().noalias() = xm->eigenMat().transpose() * ym->eigenMat();
      }
    }
    else
    {
      if (transposeY)
      {
        eigenMat().noalias() = xm->eigenMat() * ym->eigenMat().transpose();
      }
      else
      {
        eigenMat().noalias() = xm->eigenMat() * ym->eigenMat();
      }
    }
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
void MatrixDense::prodNormMatMatInPlace(const AMatrix* a,
                                        const AMatrix* m,
                                        bool transpose)
{
  if (getFlagMatrixCheck() &&
      !_isMatrixCompatible("MatrixSparse::prodNormMatMatInPlace",
                           a, 0, transpose,
                           m, 0, false,
                           a, 0, !transpose)) return;

  const auto* am = dynamic_cast<const MatrixDense*>(a);
  const auto* mm = dynamic_cast<const MatrixDense*>(m);
  if (am == nullptr || mm == nullptr)
  {
    AMatrix::prodNormMatMatInPlace(a, m, transpose);
  }
  else
  {
    if (transpose)
    {
      eigenMat().noalias() = am->eigenMat().transpose() * mm->eigenMat() * am->eigenMat();
    }
    else
    {
      eigenMat().noalias() = am->eigenMat() * mm->eigenMat() * am->eigenMat().transpose();
    }
  }
}

/**
 * Product 't(A)' %*% ['vec'] %*% 'A' or 'A' %*% ['vec'] %*% 't(A)' stored in 'this'
 *
 * @param a Input matrix
 * @param vec Input vector
 * @param transpose When True, the input Matrix is transposed
 */
void MatrixDense::prodNormMatVecInPlace(const AMatrix* a, const VectorDouble& vec, bool transpose)
{
  if (getFlagMatrixCheck() &&
      !_isMatrixCompatible("MatrixDense::prodNormMatVecInPlace",
                           a, 0, transpose,
                           a, 0, !transpose)) return;

  const auto* am = dynamic_cast<const MatrixDense*>(a);
  if (am == nullptr)
  {
    AMatrix::prodNormMatVecInPlace(a, vec, transpose);
  }
  else
  {
    if (transpose)
    {
      Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), vec.size());
      eigenMat().noalias() = am->eigenMat().transpose() * vecm * am->eigenMat();
    }
    else
    {
      Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), vec.size());
      eigenMat().noalias() = am->eigenMat() * vecm * am->eigenMat().transpose();
    }
  }
}

void MatrixDense::prodNormMatInPlace(const AMatrix* a, bool transpose)
{
  if (getFlagMatrixCheck() &&
      !_isMatrixCompatible("MatrixDense::prodNormMatInPlace",
                           a, 0, transpose,
                           a, 0, !transpose)) return;

  const auto* am = dynamic_cast<const MatrixDense*>(a);
  if (am == nullptr)
  {
    AMatrix::prodNormMatInPlace(a, transpose);
  }
  else
  {
    if (transpose)
    {
      eigenMat() = am->eigenMat().transpose() * am->eigenMat();
    }
    else
    {
      eigenMat() = am->eigenMat() * am->eigenMat().transpose();
    }
  }
}

void MatrixDense::linearCombination(double val1,
                                    const AMatrix* mat1,
                                    double val2,
                                    const AMatrix* mat2,
                                    double val3,
                                    const AMatrix* mat3)
{
  const auto* mmat1 = dynamic_cast<const MatrixDense*>(mat1);
  const auto* mmat2 = dynamic_cast<const MatrixDense*>(mat2);
  const auto* mmat3 = dynamic_cast<const MatrixDense*>(mat3);

  if ((mat1 != nullptr && mmat1 == nullptr) ||
      (mat2 != nullptr && mmat2 == nullptr) || (mat2 == this) ||
      (mat3 != nullptr && mmat3 == nullptr) || (mat3 == this))
  {
    AMatrix::linearCombination(val1, mat1, val2, mat2, val3, mat3);
  }
  else
  {
    if (mat1 != nullptr && val1 != 0.)
      eigenMat() = val1 * mmat1->eigenMat();
    if (mat2 != nullptr && val2 != 0.)
      eigenMat() += val2 * mmat2->eigenMat();
    if (mat3 != nullptr && val3 != 0.)
      eigenMat() += val3 * mmat3->eigenMat();
  }
}

void MatrixDense::fill(double value)
{
  eigenMat().setConstant(value);
}

/*! Multiply a Matrix row-wise */
void MatrixDense::multiplyRow(const VectorDouble& vec)
{
  if (getFlagMatrixCheck() && getNRows() != static_cast<Id>(vec.size()))
  {
    messerr("The size of 'vec' must match the number of rows. Nothing is done");
    return;
  }
  Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), getNCols());
  eigenMat() = vecm.asDiagonal() * eigenMat();
}

/*! Multiply a Matrix column-wise */
void MatrixDense::multiplyColumn(const VectorDouble& vec)
{
  if (getFlagMatrixCheck() && getNCols() != static_cast<Id>(vec.size()))
  {
    messerr("The size of 'vec' must match the number of columns. Nothing is done");
    return;
  }
  Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), getNRows());
  eigenMat() = eigenMat() * vecm.asDiagonal();
}

/*! Divide a Matrix row-wise */
void MatrixDense::divideRow(const VectorDouble& vec)
{
  if (getFlagMatrixCheck() && getNRows() != static_cast<Id>(vec.size()))
  {
    messerr("The size of 'vec' must match the number of rows. Nothing is done");
    return;
  }
  VectorDouble temp = VH::inverse(vec);
  Eigen::Map<const Eigen::VectorXd> vecm(temp.data(), getNCols());
  eigenMat() = vecm.asDiagonal() * eigenMat();
}

/*! Divide a Matrix column-wise */
void MatrixDense::divideColumn(const VectorDouble& vec)
{
  if (getFlagMatrixCheck() && getNCols() != static_cast<Id>(vec.size()))
  {
    messerr("The size of 'vec' must match the number of columns. Nothing is done");
    return;
  }
  VectorDouble temp = VH::inverse(vec);
  Eigen::Map<const Eigen::VectorXd> vecm(temp.data(), getNRows());
  eigenMat() = eigenMat() * vecm.asDiagonal();
}

/*! Extract a Row */
VectorDouble MatrixDense::getRow(Id irow) const
{
  VectorDouble res(getNCols());
  for (size_t i = 0; i < res.size(); ++i)
  {
    res[i] = eigenMat().row(irow)[i];
  }
  return res;
}

/*! Extract a Column */
VectorDouble MatrixDense::getColumn(Id icol) const
{
  VectorDouble res(getNRows());
  for (size_t i = 0; i < res.size(); ++i)
  {
    res[i] = eigenMat().col(icol)[i];
  }
  return res;
}

constvect MatrixDense::getViewOnColumn(Id icol) const
{
  constvect res(eigenMat().col(icol).data(), getNRows());
  return res;
}

vect MatrixDense::getViewOnColumnModify(Id icol)
{
  vect res(eigenMat().col(icol).data(), getNRows());
  return res;
}
Id MatrixDense::_terminateEigen(const Eigen::VectorXd& eigenValues,
                                 const Eigen::MatrixXd& eigenVectors,
                                 bool optionPositive,
                                 bool changeOrder)
{
  auto nrows = getNRows();
  auto ncols = getNCols();

  _eigenValues                                                         = VectorDouble(nrows);
  Eigen::Map<Eigen::VectorXd>(_eigenValues.data(), eigenValues.size()) = eigenValues;

  if (changeOrder)
    std::reverse(_eigenValues.begin(), _eigenValues.end());

  delete _eigenVectors;

  VectorDouble vec(nrows * ncols);
  Eigen::Map<Eigen::MatrixXd>(vec.data(), nrows, ncols) = eigenVectors;

  _eigenVectors = MatrixSquare::createFromVD(vec, nrows, false, changeOrder);

  if (optionPositive)
    _eigenVectors->makePositiveColumn();

  return 0;
}

Id MatrixDense::_computeGeneralizedEigen(const MatrixSymmetric& b, bool optionPositive)
{
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(eigenMat(), b.eigenMat());
  Eigen::VectorXd eigenValues  = solver.eigenvalues().real();
  Eigen::MatrixXd eigenVectors = solver.eigenvectors().real();

  return _terminateEigen(eigenValues, eigenVectors, optionPositive, true);
}

Id MatrixDense::_computeEigen(bool optionPositive)
{
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(eigenMat());
  Eigen::VectorXd eigenValues  = solver.eigenvalues().real();
  Eigen::MatrixXd eigenVectors = solver.eigenvectors().real();

  return _terminateEigen(eigenValues, eigenVectors, optionPositive, true);
}

bool MatrixDense::_needToReset(Id nrows, Id ncols)
{
  Id newsize = nrows * ncols;

  if (newsize > _maxSize)
  {
    _maxSize = newsize;
    return true;
  }
  return false;
}

void MatrixDense::_recopy(const MatrixDense& r)
{
  _eigenMatrix = r._eigenMatrix;
  _eigenValues = r._eigenValues;
  delete _eigenVectors;
  _eigenVectors = nullptr;
  if (r._eigenVectors != nullptr)
  {
    _eigenVectors = r._eigenVectors->clone();
  }
}

MatrixDense* MatrixDense::create(const MatrixDense* mat)
{
  return new MatrixDense(*mat);
}
MatrixDense* MatrixDense::create(Id nrow, Id ncol)
{
  return new MatrixDense(nrow, ncol);
}

/**
 * Converts a VectorVectorDouble into a Matrix
 * Note: the input argument is stored by row (if coming from [] specification)
 * @param  X Input VectorVectorDouble argument
 * @return The returned rectangular matrix
 *
 * @remark: the matrix is transposed implicitly while reading
 */
MatrixDense* MatrixDense::createFromVVD(const VectorVectorDouble& X)
{
  Id nrow = static_cast<Id>(X.size());
  Id ncol = static_cast<Id>(X[0].size());

  auto* mat = new MatrixDense(nrow, ncol);
  mat->_fillFromVVD(X);
  return mat;
}

MatrixDense* MatrixDense::createFromVD(const VectorDouble& X,
                                       Id nrow,
                                       Id ncol,
                                       bool byCol,
                                       bool invertColumnOrder)
{
  if (nrow * ncol != static_cast<Id>(X.size()))
  {
    messerr("Inconsistency between arguments 'nrow'(%d) and 'ncol'(%d)", nrow, ncol);
    messerr("and the dimension of the input Vector (%d)", static_cast<Id>(X.size()));
  }
  auto* mat = new MatrixDense(nrow, ncol);

  Id lec = 0;
  if (byCol)
  {
    for (Id irow = 0; irow < nrow; irow++)
      for (Id icol = 0; icol < ncol; icol++)
      {
        Id jcol = (invertColumnOrder) ? ncol - icol - 1 : icol;
        mat->setValue(irow, jcol, X[lec++]);
      }
  }
  else
  {
    for (Id icol = 0; icol < ncol; icol++)
      for (Id irow = 0; irow < nrow; irow++)
      {
        Id jcol = (invertColumnOrder) ? ncol - icol - 1 : icol;
        mat->setValue(irow, jcol, X[lec++]);
      }
  }
  return mat;
}

MatrixDense* MatrixDense::glue(const AMatrix* A1,
                               const AMatrix* A2,
                               bool flagShiftRow,
                               bool flagShiftCol)
{
  // Create the new matrix
  Id shiftRow = (flagShiftRow) ? A1->getNRows() : 0;
  Id shiftCol = (flagShiftCol) ? A1->getNCols() : 0;

  Id nrows = (flagShiftRow) ? A1->getNRows() + A2->getNRows() : MAX(A1->getNRows(), A2->getNRows());
  Id ncols = (flagShiftCol) ? A1->getNCols() + A2->getNCols() : MAX(A1->getNCols(), A2->getNCols());

  auto* mat = new MatrixDense(nrows, ncols);
  mat->fill(0.);

  // Copy the first input matrix

  for (Id irow = 0; irow < A1->getNRows(); irow++)
    for (Id icol = 0; icol < A1->getNCols(); icol++)
      mat->setValue(irow, icol, A1->getValue(irow, icol));

  // Copy the second input matrix

  for (Id irow = 0; irow < A2->getNRows(); irow++)
    for (Id icol = 0; icol < A2->getNCols(); icol++)
      mat->setValue(irow + shiftRow, icol + shiftCol, A2->getValue(irow, icol));

  return mat;
}

void MatrixDense::addRow(Id nrow_added)
{
  auto nrows = getNRows();
  auto ncols = getNCols();

  AMatrix* statsSave = this->clone();
  reset(nrows + nrow_added, ncols);
  for (Id irow = 0; irow < nrows; irow++)
    for (Id icol = 0; icol < ncols; icol++)
      setValue(irow, icol, statsSave->getValue(irow, icol));
  delete statsSave;
}

void MatrixDense::addColumn(Id ncolumn_added)
{
  auto nrows = getNRows();
  auto ncols = getNCols();

  AMatrix* statsSave = this->clone();
  reset(nrows, ncols + ncolumn_added);
  for (Id irow = 0; irow < nrows; irow++)
    for (Id icol = 0; icol < ncols; icol++)
      setValue(irow, icol, statsSave->getValue(irow, icol));
  delete statsSave;
}

/**
 * @brief Create an output Rectangular Matrix by selecting some rows and columns
 *        of the Input matrix 'A'
 *
 * @param res      Output Rectangular Matrix
 * @param A        Input Rectangular Matrix
 * @param rowKeep  Set of Rows to be kept (all if not defined)
 * @param colKeep  Set of Columns to be kept (all if not defined)
 * @param flagInvertRow when True, transform 'rowKeep' into 'rowDrop'
 * @param flagInvertCol when True, transform 'colKeep' into 'colDrop'
 */
bool MatrixDense::sample(MatrixDense& res,
                         const AMatrix& A,
                         const VectorInt& rowKeep,
                         const VectorInt& colKeep,
                         bool flagInvertRow,
                         bool flagInvertCol)
{
  auto nrowtotal = A.getNRows();
  auto ncoltotal = A.getNCols();
  VectorInt rows = rowKeep;
  if (rows.empty()) rows = VH::sequence(nrowtotal);
  if (flagInvertRow) rows = VH::complement(VH::sequence(nrowtotal), rows);
  VectorInt cols = colKeep;
  if (cols.empty()) cols = VH::sequence(ncoltotal);
  if (flagInvertCol) cols = VH::complement(VH::sequence(ncoltotal), cols);

  Id nrows = static_cast<Id>(rows.size());
  Id ncols = static_cast<Id>(cols.size());
  if (nrows <= 0 || ncols <= 0) return false;

  for (Id irow = 0; irow < nrows; irow++)
  {
    if (!checkArg("Selected Row index", rows[irow], nrowtotal)) return false;
  }
  for (Id icol = 0; icol < ncols; icol++)
  {
    if (!checkArg("Selected Column index", cols[icol], ncoltotal)) return false;
  }

  res.resize(nrows, ncols);
  for (Id irow = 0; irow < nrows; irow++)
    for (Id icol = 0; icol < ncols; icol++)
      res.setValue(irow, icol, A.getValue(rows[irow], cols[icol]));
  return true;
}

/**
 * @brief Set the values contained in 'A' into the current matrix
 *
 * @param A Input Matrix
 * @param rowFetch Set of row indices of 'this' where values of 'A' should be stored
 * @param colFetch Set of column indices of'this' where values of 'A' should be stored
 * @param flagInvertRow when True, transform 'rowFetch' into 'rowAvoid'
 * @param flagInvertCol when True, transform 'colFetch' into 'colAvoid'
 */
void MatrixDense::unsample(const AMatrix* A,
                           const VectorInt& rowFetch,
                           const VectorInt& colFetch,
                           bool flagInvertRow,
                           bool flagInvertCol)
{
  auto nrowtotal = getNRows();
  auto ncoltotal = getNCols();
  VectorInt rows = rowFetch;
  if (rows.empty()) rows = VH::sequence(nrowtotal);
  if (flagInvertRow) rows = VH::complement(VH::sequence(nrowtotal), rows);
  VectorInt cols = colFetch;
  if (cols.empty()) cols = VH::sequence(ncoltotal);
  if (flagInvertCol) cols = VH::complement(VH::sequence(ncoltotal), cols);

  Id nrows = static_cast<Id>(rows.size());
  Id ncols = static_cast<Id>(cols.size());
  if (nrows <= 0 || ncols <= 0) return;

  for (Id irow = 0; irow < nrows; irow++)
  {
    if (!checkArg("Selected Row index", rows[irow], getNRows()))
      return;
  }
  for (Id icol = 0; icol < ncols; icol++)
  {
    if (!checkArg("Selected Column index", cols[icol], getNCols()))
      return;
  }

  for (Id irow = 0; irow < nrows; irow++)
    for (Id icol = 0; icol < ncols; icol++)
      setValue(rows[irow], cols[icol], A->getValue(irow, icol));
}

/**
 * @brief Perform the compressing product, according to 'transpose'
 * - False: 'this'[nrows,ncols] %*% t('matLC')[ncolsCL,nrowsCL] -> mat[nrows,ncolsCL]
 * - True:  t('matLC')[ncolsCL,nrowsCL] %*% 'this'[nrows,ncols] -> mat[ncolsCL,ncols]
 *
 * @param matLC
 * @param transpose
 * @return MatrixDense
 */
MatrixDense MatrixDense::compressMatLC(const MatrixDense& matLC,
                                       bool transpose)
{
  auto nrows = getNRows();
  auto ncols = getNCols();
  auto nrowCL = matLC.getNRows();
  auto ncolCL = matLC.getNCols();
  MatrixDense mat;

  if (!transpose)
  {
    if (ncols != ncolCL)
    {
      messerr("Number of Columns (%d) should match number of columns of 'matLC' (%d)",
              ncols, ncolCL);
      return mat;
    }
    mat.resize(nrows, nrowCL);
    for (Id irow = 0; irow < nrows; irow++)
      for (Id irowCL = 0; irowCL < nrowCL; irowCL++)
      {
        double value = 0.;
        for (Id k = 0; k < ncols; k++)
          value += getValue(irow, k) * matLC.getValue(irowCL, k);
        mat.setValue(irow, irowCL, value);
      }
  }
  else
  {
    if (ncolCL != nrows)
    {
      messerr("Number of Rows (%d) should match number of Columns of 'matLC' (%d)",
              nrows, ncolCL);
      return mat;
    }
    mat.resize(nrowCL, ncols);
    for (Id irowCL = 0; irowCL < nrowCL; irowCL++)
      for (Id icol = 0; icol < ncols; icol++)
      {
        double value = 0.;
        for (Id k = 0; k < nrows; k++)
          value += matLC.getValue(irowCL, k) * getValue(k, icol);
        mat.setValue(irowCL, icol, value);
      }
  }

  return mat;
}

void MatrixDense::sum(const MatrixDense* mat1,
                      const MatrixDense* mat2,
                      MatrixDense* mat3)
{
  mat3->eigenMat().noalias() = mat1->eigenMat() + mat2->eigenMat();
}

} // namespace gstlrn
