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
#include "geoslib_old_f.h"

#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/AMatrix.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Utilities.hpp"

MatrixRectangular::MatrixRectangular(int nrows, int ncols, int opt_eigen)
    : AMatrixDense(nrows, ncols, opt_eigen),
      _rectMatrix()
{
  _allocate();
}

MatrixRectangular::MatrixRectangular(const MatrixRectangular &r)
  : AMatrixDense(r),
    _rectMatrix()
{
  _recopyLocal(r);
}

MatrixRectangular::MatrixRectangular(const AMatrix &m)
    : AMatrixDense(m),
      _rectMatrix()
{
  const MatrixRectangular* matrixLoc = dynamic_cast<const MatrixRectangular*>(&m);
  if (matrixLoc != nullptr)
    _recopyLocal(*matrixLoc);
  else
  {
    _allocate();
    AMatrix::copyElements(m);
  }
}

MatrixRectangular& MatrixRectangular::operator= (const MatrixRectangular &r)
{
  if (this != &r)
  {
    AMatrixDense::operator=(r);
    _recopyLocal(r);
  }
  return *this;
}

MatrixRectangular::~MatrixRectangular()
{
  _deallocate();
}

/**
 * Converts a VectorVectorDouble into a Matrix
 * Note: the input argument is stored by row (if coming from [] specification)
 * @param  X Input VectorVectorDouble argument
 * @param opt_eigen Option for use of Eigen Library
 * @return The returned rectangular matrix
 *
 * @remark: the matrix is transposed implicitly while reading
 */
MatrixRectangular* MatrixRectangular::createFromVVD(const VectorVectorDouble& X, int opt_eigen)
{
  int nrow = (int) X.size();
  int ncol = (int) X[0].size();

  MatrixRectangular* mat = new MatrixRectangular(nrow, ncol, opt_eigen);
  mat->_fillFromVVD(X);
  return mat;
}

MatrixRectangular* MatrixRectangular::createFromVD(const VectorDouble &X,
                                                   int nrow,
                                                   int ncol,
                                                   bool byCol,
                                                   int opt_eigen,
                                                   bool invertColumnOrder)
{
  if (nrow * ncol != (int) X.size())
  {
    messerr("Inconsistency between arguments 'nrow'(%d) and 'ncol'(%d)", nrow, ncol);
    messerr("and the dimension of the input Vector (%d)", (int) X.size());
  }
  MatrixRectangular* mat = new MatrixRectangular(nrow, ncol, opt_eigen);

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

double MatrixRectangular::_getValueByRank(int irank) const
{
  if (isFlagEigen())
    return AMatrixDense::_getValueByRank(irank);
  else
    return _getValueByRankLocal(irank);
}

void MatrixRectangular::_setValueByRank(int irank, double value)
{
  if (isFlagEigen())
    AMatrixDense::_setValueByRank(irank, value);
  else
    _setValueByRankLocal(irank, value);
}

void MatrixRectangular::_prodMatVecInPlacePtr(const double *x, double *y, bool transpose) const
{
  if (isFlagEigen())
    AMatrixDense::_prodMatVecInPlacePtr(x, y, transpose);
  else
    _prodMatVecInPlacePtrLocal(x, y, transpose);
}

void MatrixRectangular::_prodVecMatInPlacePtr(const double *x, double *y, bool transpose) const
{
  if (isFlagEigen())
    AMatrixDense::_prodVecMatInPlacePtr(x, y, transpose);
  else
    _prodVecMatInPlacePtrLocal(x, y, transpose);
}

void MatrixRectangular::_transposeInPlace()
{
  if (isFlagEigen())
    AMatrixDense::_transposeInPlace();
  else
    _transposeInPlaceLocal();
}

void MatrixRectangular::_allocate()
{
  if (isFlagEigen())
    AMatrixDense::_allocate();
  else
    _allocateLocal();
  fill(0.);
}

int MatrixRectangular::_getIndexToRank(int irow, int icol) const
{
  if (isFlagEigen())
    return AMatrixDense::_getIndexToRank(irow, icol);
  else
    return _getIndexToRankLocal(irow, icol);
}

int MatrixRectangular::_invert()
{
  my_throw("Invert method is limited to Square matrices");
  return 0;
}

int MatrixRectangular::_solve(const VectorDouble& /*b*/, VectorDouble& /*x*/) const
{
  my_throw("Invert method is limited to Square Symmetrical Matrices");
  return 0;
}

void MatrixRectangular::addRow(int nrow_added)
{
  int nrows = getNRows();
  int ncols = getNCols();

  AMatrix* statsSave = this->clone();
  reset(nrows+nrow_added, ncols);
  for (int irow=0; irow< nrows; irow++)
    for (int icol=0; icol<ncols; icol++)
      setValue(irow, icol, statsSave->getValue(irow, icol));
}

void MatrixRectangular::addColumn(int ncolumn_added)
{
  int nrows = getNRows();
  int ncols = getNCols();

  AMatrix* statsSave = this->clone();
  reset(nrows, ncols+ncolumn_added);
  for (int irow=0; irow< nrows; irow++)
    for (int icol=0; icol<ncols; icol++)
      setValue(irow, icol, statsSave->getValue(irow, icol));
}

int MatrixRectangular::_getMatrixPhysicalSize() const
{
  if (isFlagEigen())
    return AMatrixDense::_getMatrixPhysicalSize();
  else
    return AMatrix::_getMatrixPhysicalSize();
}

double& MatrixRectangular::_getValueRef(int irow, int icol)
{
  if (isFlagEigen())
    return AMatrixDense::_getValueRef(irow, icol);
  else
    return _getValueRefLocal(irow, icol);
}

/// ========================================================================
/// The subsequent methods rely on the specific local storage ('rectMatrix')
/// ========================================================================

void MatrixRectangular::_allocateLocal()
{
  _rectMatrix.resize(_getMatrixPhysicalSize(),0.);
}

void MatrixRectangular::_recopyLocal(const MatrixRectangular& r)
{
  _rectMatrix = r._rectMatrix;
}

double MatrixRectangular::_getValue(int irow, int icol) const
{
  if (! _isIndexValid(irow,icol)) return TEST;
  int rank = _getIndexToRank(irow,icol);
  return _rectMatrix[rank];
}

double MatrixRectangular::_getValueByRankLocal(int irank) const
{
  if (! _isRankValid(irank)) return TEST;
  return _rectMatrix[irank];
}

/*! Gets a reference to the value at row 'irow' and column 'icol' */
double& MatrixRectangular::_getValueRefLocal(int irow, int icol)
{
  int rank = _getIndexToRank(irow,icol);
  return _rectMatrix[rank];
}

void MatrixRectangular::_setValueByRankLocal(int irank, double value)
{
  if (! _isRankValid(irank)) return;
  _rectMatrix[irank] = value;
}

void MatrixRectangular::_setValue(int irow, int icol, double value)
{
  if (! _isIndexValid(irow, icol)) return;
  int rank = _getIndexToRank(irow, icol);
  _rectMatrix[rank] = value;
}

void MatrixRectangular::_updValue(int irow, int icol, const EOperator& oper, double value)
{
  if (! _isIndexValid(irow, icol)) return;
  int rank = _getIndexToRank(irow, icol);
  _rectMatrix[rank] = modifyOperator(oper, _rectMatrix[rank], value);
}

void MatrixRectangular::_prodMatVecInPlacePtrLocal(const double *x, double *y, bool transpose) const
{
  if (! transpose)
    matrix_product_safe(getNRows(), getNCols(), 1, _rectMatrix.data(), x, y);
  else
    matrix_product_safe(1, getNRows(), getNCols(), x, _rectMatrix.data(), y);
}

void MatrixRectangular::_prodVecMatInPlacePtrLocal(const double *x, double *y, bool transpose) const
{
  if (! transpose)
    matrix_product_safe(1, getNRows(), getNCols(), x, _rectMatrix.data(), y);
  else
    matrix_product_safe(getNRows(), getNCols(), 1, _rectMatrix.data(), x, y);
}

void MatrixRectangular::_transposeInPlaceLocal()
{
  VectorDouble old;
  old.resize(getNRows() * getNCols());
  matrix_transpose(getNRows(), getNCols(), _rectMatrix, old);
  _rectMatrix = old;
  int temp = getNCols();
  _setNCols(getNRows());
  _setNRows(temp);
}

int MatrixRectangular::_getIndexToRankLocal(int irow, int icol) const
{
  return (icol * getNRows() + irow);
}

MatrixRectangular* MatrixRectangular::glue(const AMatrix *A1,
                                           const AMatrix *A2,
                                           bool flagShiftRow,
                                           bool flagShiftCol)
{
  // Create the new matrix
  int shiftRow = (flagShiftRow) ? A1->getNRows() : 0;
  int shiftCol = (flagShiftCol) ? A1->getNCols() : 0;

  int nrows = (flagShiftRow) ? A1->getNRows() + A2->getNRows() : MAX(A1->getNRows(), A2->getNRows());
  int ncols = (flagShiftCol) ? A1->getNCols() + A2->getNCols() : MAX(A1->getNCols(), A2->getNCols());

  MatrixRectangular* mat = new MatrixRectangular(nrows, ncols);
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
