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
#include "Matrix/MatrixSquare.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Matrix/AMatrix.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Utilities.hpp"
#include "geoslib_define.h"

#include <math.h>
#include <omp.h>

MatrixDense::MatrixDense(int nrow, int ncol)
  : AMatrix(nrow, ncol)
  , _flagEigenDecompose(false)
  , _eigenValues()
  , _eigenVectors(nullptr)
  , _eigenMatrix()
{
  _allocate();
}

MatrixDense::MatrixDense(const MatrixDense& r)
  : AMatrix(r)
  , _flagEigenDecompose(false)
  , _eigenValues()
  , _eigenVectors(nullptr)
  , _eigenMatrix()
{
  _allocate();
  _recopy(r);
}

MatrixDense::MatrixDense(const AMatrix& r)
  : AMatrix(r)
  , _flagEigenDecompose(false)
  , _eigenValues()
  , _eigenVectors(nullptr)
  , _eigenMatrix()
{
  _allocate();
  copyElements(r);
}

MatrixDense& MatrixDense::operator=(const MatrixDense& r)
{
  if (this != &r)
  {
    AMatrix::operator=(r);
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

double MatrixDense::getValue(int irow, int icol, bool flagCheck) const
{
  if (flagCheck && !_isIndexValid(irow, icol)) return TEST;
  return getEigenMat()(irow, icol);
}

double MatrixDense::_getValueByRank(int irank) const
{
  return *(getEigenMat().data() + irank);
}

constvect MatrixDense::getColumnPtr(int icol) const
{
  const auto a = getEigenMat().col(icol);
  size_t n     = getNRows();
  return {a.data(), n};
}
void MatrixDense::_setValueByRank(int irank, double value)
{
  *(getEigenMat().data() + irank) = value;
}

void MatrixDense::setValue(int irow, int icol, double value, bool flagCheck)
{
  if (flagCheck && !_isIndexValid(irow, icol)) return;
  getEigenMat()(irow, icol) = value;
  if (mustBeSymmetric() && irow != icol) getEigenMat()(icol, irow) = value;
}

void MatrixDense::updValue(int irow,
                           int icol,
                           const EOperator& oper,
                           double value,
                           bool flagCheck)
{
  if (flagCheck && !_isIndexValid(irow, icol)) return;
  double result             = modifyOperator(oper, getEigenMat()(irow, icol), value);
  getEigenMat()(irow, icol) = result;
  if (mustBeSymmetric() && irow != icol)
    getEigenMat()(icol, irow) = result;
}

double& MatrixDense::_getValueRef(int irow, int icol)
{
  return *(getEigenMat().data() + _getIndexToRank(irow, icol));
}

int MatrixDense::_getMatrixPhysicalSize() const
{
  return getEigenMat().size();
}

// Default storage in Eigen is column-major (see https://eigen.tuxfamily.org/dox/group__TopicStorageOrders.html)
int MatrixDense::_getIndexToRank(int irow, int icol) const
{
  return (icol * getNRows() + irow);
}

void MatrixDense::_transposeInPlace()
{
  int nrows = getNRows();
  int ncols = getNCols();
  getEigenMat().transposeInPlace();
  _setNCols(nrows);
  _setNRows(ncols);
}

void MatrixDense::_addProdMatVecInPlaceToDestPtr(const double* x, double* y, bool transpose) const
{
  Eigen::Map<const Eigen::VectorXd> xm(x, getNCols());
  Eigen::Map<Eigen::VectorXd> ym(y, getNRows());
  if (transpose)
    ym.noalias() += getEigenMat().transpose() * xm;
  else
    ym.noalias() += getEigenMat() * xm;
}

// TODO supress this method and implement it in the virtual class AMatrix
void MatrixDense::_prodMatVecInPlacePtr(const double* x, double* y, bool transpose) const
{
  int nc = transpose ? getNRows() : getNCols();
  int nr = transpose ? getNCols() : getNRows();
  Eigen::Map<const Eigen::VectorXd> xm(x, nc);
  Eigen::Map<Eigen::VectorXd> ym(y, nr);
  if (transpose)
    ym.noalias() = getEigenMat().transpose() * xm;
  else
    ym.noalias() = getEigenMat() * xm;
}

void MatrixDense::_prodVecMatInPlacePtr(const double* x, double* y, bool transpose) const
{
  Eigen::Map<const Eigen::VectorXd> xm(x, getNRows());
  Eigen::Map<Eigen::VectorXd> ym(y, getNCols());
  if (transpose)
    ym.noalias() = xm.transpose() * getEigenMat().transpose();
  else
    ym.noalias() = xm.transpose() * getEigenMat();
}

int MatrixDense::_invert()
{
  /// TODO : check beforehand if matrix is invertible ?
  getEigenMat() = getEigenMat().inverse();
  return 0;
}

int MatrixDense::invert2(MatrixDense& res) const
{
  /// TODO : check beforehand if matrix is invertible ?
  res.getEigenMat() = getEigenMat().inverse();
  return 0;
}
int MatrixDense::_solve(const VectorDouble& b, VectorDouble& x) const
{
  /// TODO : check beforehand if matrix is invertible ?
  Eigen::Map<const Eigen::VectorXd> bm(b.data(), getNCols());
  Eigen::Map<Eigen::VectorXd> xm(x.data(), getNRows());
  xm = getEigenMat().inverse() * bm;
  return 0;
}

void MatrixDense::setColumn(int icol, const VectorDouble& tab, bool flagCheck)
{
  if (flagCheck)
  {
    if (!_isColumnValid(icol)) return;
    if (!_isColumnSizeConsistent(tab)) return;
  }
  Eigen::Map<const Eigen::VectorXd> tabm(tab.data(), getNRows());
  getEigenMat().col(icol) = tabm;
}

void MatrixDense::setRow(int irow, const VectorDouble& tab, bool flagCheck)
{
  if (flagCheck)
  {
    if (!_isRowValid(irow)) return;
    if (!_isRowSizeConsistent(tab)) return;
  }
  Eigen::Map<const Eigen::VectorXd> tabm(tab.data(), getNCols());
  getEigenMat().row(irow) = tabm;
}

void MatrixDense::setDiagonal(const VectorDouble& tab, bool flagCheck)
{
  if (flagCheck)
  {
    if (!_isRowSizeConsistent(tab)) return;
  }
  getEigenMat().setZero();
  Eigen::Map<const Eigen::VectorXd> tabm(tab.data(), getNRows());
  getEigenMat().diagonal() = tabm;
}

void MatrixDense::setDiagonalToConstant(double value)
{
  getEigenMat().setZero();
  getEigenMat().diagonal() = Eigen::VectorXd::Constant(getNRows(), value);
}

void MatrixDense::addScalar(double v)
{
  getEigenMat().array() += v;
}

void MatrixDense::addScalarDiag(double v)
{
  getEigenMat().diagonal() += Eigen::VectorXd::Constant(getNRows(), v);
}

void MatrixDense::prodScalar(double v)
{
  getEigenMat().array() *= v;
}

void MatrixDense::addMatInPlace(const MatrixDense& y, double cx, double cy)
{
  getEigenMat().noalias() = cx * getEigenMat() + cy * y.getEigenMat();
}

void MatrixDense::prodMatMatInPlace(const AMatrix* x,
                                    const AMatrix* y,
                                    bool transposeX,
                                    bool transposeY)
{
  const MatrixDense* xm = dynamic_cast<const MatrixDense*>(x);
  const MatrixDense* ym = dynamic_cast<const MatrixDense*>(y);
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

void MatrixDense::prodMatMatInPlaceOptim(const MatrixDense* x,
                                         const MatrixDense* y,
                                         bool transposeX,
                                         bool transposeY)
{

  if (transposeX)
  {
    if (transposeY)
    {
      getEigenMat().noalias() = x->getEigenMat().transpose() * y->getEigenMat().transpose();
    }
    else
    {
      getEigenMat().noalias() = x->getEigenMat().transpose() * y->getEigenMat();
    }
  }
  else
  {
    if (transposeY)
    {
      getEigenMat().noalias() = x->getEigenMat() * y->getEigenMat().transpose();
    }
    else
    {
      getEigenMat().noalias() = x->getEigenMat() * y->getEigenMat();
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
void MatrixDense::prodNormMatMatInPlace(const MatrixDense* a,
                                        const MatrixDense* m,
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
void MatrixDense::prodNormMatVecInPlace(const MatrixDense& a, const VectorDouble& vec, bool transpose)
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

void MatrixDense::fill(double value)
{
  getEigenMat().setConstant(value);
}

/*! Multiply a Matrix row-wise */
void MatrixDense::multiplyRow(const VectorDouble& vec)
{
  Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), getNCols());
  getEigenMat() = vecm.asDiagonal() * getEigenMat();
}

/*! Multiply a Matrix column-wise */
void MatrixDense::multiplyColumn(const VectorDouble& vec)
{
  Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), getNRows());
  getEigenMat() = getEigenMat() * vecm.asDiagonal();
}

/*! Divide a Matrix row-wise */
void MatrixDense::divideRow(const VectorDouble& vec)
{
  VectorDouble temp = VH::inverse(vec);
  Eigen::Map<const Eigen::VectorXd> vecm(temp.data(), getNCols());
  getEigenMat() = vecm.asDiagonal() * getEigenMat();
}

/*! Divide a Matrix column-wise */
void MatrixDense::divideColumn(const VectorDouble& vec)
{
  VectorDouble temp = VH::inverse(vec);
  Eigen::Map<const Eigen::VectorXd> vecm(temp.data(), getNRows());
  getEigenMat() = getEigenMat() * vecm.asDiagonal();
}

/*! Perform 'vec' * 'this' */
VectorDouble MatrixDense::prodVecMat(const VectorDouble& x, bool transpose) const
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
VectorDouble MatrixDense::prodMatVec(const VectorDouble& x, bool transpose) const
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
VectorDouble MatrixDense::getRow(int irow) const
{
  VectorDouble res(getNCols());
  for (size_t i = 0; i < res.size(); ++i)
  {
    res[i] = getEigenMat().row(irow)[i];
  }
  return res;
}

/*! Extract a Column */
VectorDouble MatrixDense::getColumn(int icol) const
{
  VectorDouble res(getNRows());
  for (size_t i = 0; i < res.size(); ++i)
  {
    res[i] = getEigenMat().col(icol)[i];
  }
  return res;
}

constvect MatrixDense::getViewOnColumn(int icol) const
{
  constvect res(getEigenMat().col(icol).data(), getNRows());
  return res;
}

vect MatrixDense::getViewOnColumnModify(int icol)
{
  vect res(getEigenMat().col(icol).data(), getNRows());
  return res;
}
int MatrixDense::_terminateEigen(const Eigen::VectorXd& eigenValues,
                                 const Eigen::MatrixXd& eigenVectors,
                                 bool optionPositive,
                                 bool changeOrder)
{
  int nrows = getNRows();
  int ncols = getNCols();

  _eigenValues                                                         = VectorDouble(nrows);
  Eigen::Map<Eigen::VectorXd>(_eigenValues.data(), eigenValues.size()) = eigenValues;

  if (changeOrder)
    std::reverse(_eigenValues.begin(), _eigenValues.end());

  delete _eigenVectors;

  VectorDouble vec(nrows * ncols);
  Eigen::Map<Eigen::MatrixXd>(vec.data(), nrows, ncols) = eigenVectors;

  _eigenVectors = MatrixSquare::createFromVD(vec, nrows, false, changeOrder);

  if (optionPositive) _eigenVectors->makePositiveColumn();

  _flagEigenDecompose = true;

  return 0;
}

int MatrixDense::_computeGeneralizedEigen(const MatrixSymmetric& b, bool optionPositive)
{
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(getEigenMat(), b.getEigenMat());
  Eigen::VectorXd eigenValues  = solver.eigenvalues().real();
  Eigen::MatrixXd eigenVectors = solver.eigenvectors().real();

  return _terminateEigen(eigenValues, eigenVectors, optionPositive, true);
}

int MatrixDense::_computeEigen(bool optionPositive)
{
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(getEigenMat());
  Eigen::VectorXd eigenValues  = solver.eigenvalues().real();
  Eigen::MatrixXd eigenVectors = solver.eigenvectors().real();

  return _terminateEigen(eigenValues, eigenVectors, optionPositive, true);
}

bool MatrixDense::_needToReset(int nrows, int ncols)
{
  int newsize = nrows * ncols;

  if (newsize > _maxSize)
  {
    _maxSize = newsize;
    return true;
  }
  return false;
}

void MatrixDense::_recopy(const MatrixDense& r)
{
  _eigenMatrix        = r._eigenMatrix;
  _flagEigenDecompose = r._flagEigenDecompose;
  _eigenValues        = r._eigenValues;
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
MatrixDense* MatrixDense::create(int nrow, int ncol)
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
  int nrow = (int)X.size();
  int ncol = (int)X[0].size();

  MatrixDense* mat = new MatrixDense(nrow, ncol);
  mat->_fillFromVVD(X);
  return mat;
}

MatrixDense* MatrixDense::createFromVD(const VectorDouble& X,
                                       int nrow,
                                       int ncol,
                                       bool byCol,
                                       bool invertColumnOrder)
{
  if (nrow * ncol != (int)X.size())
  {
    messerr("Inconsistency between arguments 'nrow'(%d) and 'ncol'(%d)", nrow, ncol);
    messerr("and the dimension of the input Vector (%d)", (int)X.size());
  }
  MatrixDense* mat = new MatrixDense(nrow, ncol);

  int lec = 0;
  if (byCol)
  {
    for (int irow = 0; irow < nrow; irow++)
      for (int icol = 0; icol < ncol; icol++)
      {
        int jcol = (invertColumnOrder) ? ncol - icol - 1 : icol;
        mat->setValue(irow, jcol, X[lec++]);
      }
  }
  else
  {
    for (int icol = 0; icol < ncol; icol++)
      for (int irow = 0; irow < nrow; irow++)
      {
        int jcol = (invertColumnOrder) ? ncol - icol - 1 : icol;
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
  int shiftRow = (flagShiftRow) ? A1->getNRows() : 0;
  int shiftCol = (flagShiftCol) ? A1->getNCols() : 0;

  int nrows = (flagShiftRow) ? A1->getNRows() + A2->getNRows() : MAX(A1->getNRows(), A2->getNRows());
  int ncols = (flagShiftCol) ? A1->getNCols() + A2->getNCols() : MAX(A1->getNCols(), A2->getNCols());

  MatrixDense* mat = new MatrixDense(nrows, ncols);
  mat->fill(0.);

  // Copy the first input matrix

  for (int irow = 0; irow < A1->getNRows(); irow++)
    for (int icol = 0; icol < A1->getNCols(); icol++)
      mat->setValue(irow, icol, A1->getValue(irow, icol));

  // Copy the second input matrix

  for (int irow = 0; irow < A2->getNRows(); irow++)
    for (int icol = 0; icol < A2->getNCols(); icol++)
      mat->setValue(irow + shiftRow, icol + shiftCol, A2->getValue(irow, icol));

  return mat;
}

void MatrixDense::addRow(int nrow_added)
{
  int nrows = getNRows();
  int ncols = getNCols();

  AMatrix* statsSave = this->clone();
  reset(nrows + nrow_added, ncols);
  for (int irow = 0; irow < nrows; irow++)
    for (int icol = 0; icol < ncols; icol++)
      setValue(irow, icol, statsSave->getValue(irow, icol));
  delete statsSave;
}

void MatrixDense::addColumn(int ncolumn_added)
{
  int nrows = getNRows();
  int ncols = getNCols();

  AMatrix* statsSave = this->clone();
  reset(nrows, ncols + ncolumn_added);
  for (int irow = 0; irow < nrows; irow++)
    for (int icol = 0; icol < ncols; icol++)
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
  int nrowtotal  = A.getNRows();
  int ncoltotal  = A.getNCols();
  VectorInt rows = rowKeep;
  if (rows.empty()) rows = VH::sequence(nrowtotal);
  if (flagInvertRow) rows = VH::complement(VH::sequence(nrowtotal), rows);
  VectorInt cols = colKeep;
  if (cols.empty()) cols = VH::sequence(ncoltotal);
  if (flagInvertCol) cols = VH::complement(VH::sequence(ncoltotal), cols);

  int nrows = (int)rows.size();
  int ncols = (int)cols.size();
  if (nrows <= 0 || ncols <= 0) return false;

  for (int irow = 0; irow < nrows; irow++)
  {
    if (!checkArg("Selected Row index", rows[irow], nrowtotal)) return false;
  }
  for (int icol = 0; icol < ncols; icol++)
  {
    if (!checkArg("Selected Column index", cols[icol], ncoltotal)) return false;
  }

  res.resize(nrows, ncols);
  for (int irow = 0; irow < nrows; irow++)
    for (int icol = 0; icol < ncols; icol++)
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
  int nrowtotal  = getNRows();
  int ncoltotal  = getNCols();
  VectorInt rows = rowFetch;
  if (rows.empty()) rows = VH::sequence(nrowtotal);
  if (flagInvertRow) rows = VH::complement(VH::sequence(nrowtotal), rows);
  VectorInt cols = colFetch;
  if (cols.empty()) cols = VH::sequence(ncoltotal);
  if (flagInvertCol) cols = VH::complement(VH::sequence(ncoltotal), cols);

  int nrows = (int)rows.size();
  int ncols = (int)cols.size();
  if (nrows <= 0 || ncols <= 0) return;

  for (int irow = 0; irow < nrows; irow++)
  {
    if (!checkArg("Selected Row index", rows[irow], getNRows()))
      return;
  }
  for (int icol = 0; icol < ncols; icol++)
  {
    if (!checkArg("Selected Column index", cols[icol], getNCols()))
      return;
  }

  for (int irow = 0; irow < nrows; irow++)
    for (int icol = 0; icol < ncols; icol++)
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
  int nrows  = getNRows();
  int ncols  = getNCols();
  int nrowCL = matLC.getNRows();
  int ncolCL = matLC.getNCols();
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
    for (int irow = 0; irow < nrows; irow++)
      for (int irowCL = 0; irowCL < nrowCL; irowCL++)
      {
        double value = 0.;
        for (int k = 0; k < ncols; k++)
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
    for (int irowCL = 0; irowCL < nrowCL; irowCL++)
      for (int icol = 0; icol < ncols; icol++)
      {
        double value = 0.;
        for (int k = 0; k < nrows; k++)
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
  mat3->getEigenMat().noalias() = mat1->getEigenMat() + mat2->getEigenMat();
}
