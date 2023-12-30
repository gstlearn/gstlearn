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

MatrixRectangular::MatrixRectangular(int nrows, int ncols)
    : AMatrixDense(nrows, ncols),
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
 * @return The returned rectangular matrix
 *
 * @remark: the matrix is transposed implicitly while reading
 */
MatrixRectangular* MatrixRectangular::createFromVVD(const VectorVectorDouble& X)
{
  int nrow = (int) X.size();
  int ncol = (int) X[0].size();

  MatrixRectangular* mat = new MatrixRectangular(nrow, ncol);
  mat->_fillFromVVD(X);
  return mat;
}

MatrixRectangular* MatrixRectangular::createFromVD(const VectorDouble &X,
                                                   int nrow,
                                                   int ncol,
                                                   bool byCol)
{
  if (nrow * ncol != (int) X.size())
  {
    messerr("Inconsistency between arguments 'nrow'(%d) and 'ncol'(%d)", nrow, ncol);
    messerr("and the dimension of the input Vector (%d)", (int) X.size());
  }
  MatrixRectangular* mat = new MatrixRectangular(nrow, ncol);

  int lec = 0;
  if (byCol)
  {
    for (int irow = 0; irow < nrow; irow++)
      for (int icol = 0; icol < ncol; icol++)
        mat->setValue(irow, icol, X[lec++]);
  }
  else
  {
    for (int icol = 0; icol < ncol; icol++)
      for (int irow = 0; irow < nrow; irow++)
        mat->setValue(irow, icol, X[lec++]);
  }
  return mat;
}

double MatrixRectangular::_getValue(int irow, int icol) const
{
  if (isFlagEigen())
    return AMatrixDense::_getValue(irow, icol);
  else
    return _getValueLocal(irow, icol);
}

double MatrixRectangular::_getValueByRank(int irank) const
{
  if (isFlagEigen())
    return AMatrixDense::_getValueByRank(irank);
  else
    return _getValueLocal(irank);
}

void MatrixRectangular::_setValueByRank(int irank, double value)
{
  if (isFlagEigen())
    AMatrixDense::_setValueByRank(irank, value);
  else
    _setValueLocal(irank, value);
}

void MatrixRectangular::_setValue(int irow, int icol, double value)
{
  if (isFlagEigen())
    AMatrixDense::_setValue(irow, icol, value);
  else
    _setValueLocal(irow, icol, value);
}

void MatrixRectangular::_prodVectorInPlace(const double *inv, double *outv) const
{
  if (isFlagEigen())
    AMatrixDense::_prodVectorInPlace(inv, outv);
  else
    _prodVectorLocal(inv, outv);
}

void MatrixRectangular::_transposeInPlace()
{
  if (isFlagEigen())
    AMatrixDense::_transposeInPlace();
  else
    _transposeInPlaceLocal();
}

void MatrixRectangular::_deallocate()
{
  if (isFlagEigen())
    AMatrixDense::_deallocate();
  else
  {
    // Potential code for this class would be located here
  }
}

void MatrixRectangular::_allocate()
{
  if (isFlagEigen())
    AMatrixDense::_allocate();
  else
    _allocateLocal();
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

MatrixRectangular* MatrixRectangular::reduce(const VectorInt &validRows,
                                             const VectorInt &validCols) const
{
  // Order and shrink the input vectors
  VectorInt localValidRows = VH::filter(validRows, 0, getNRows());
  VectorInt localValidCols = VH::filter(validCols, 0, getNCols());
  int newNRows = (int) localValidRows.size();
  int newNCols = (int) localValidCols.size();
  if (newNRows <= 0)
  {
    messerr("The new Matrix has no Row left");
    return nullptr;
  }
  if (newNCols <= 0)
  {
    messerr("The new Matrix has no Column left");
    return nullptr;
  }

  MatrixRectangular* res = new MatrixRectangular(newNRows, newNCols);
  res->copyReduce(this, localValidRows, localValidCols);

  return res;
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

double MatrixRectangular::_getValueLocal(int irow, int icol) const
{
  if (! _isIndexValid(irow,icol)) return TEST;
  int rank = _getIndexToRank(irow,icol);
  return _rectMatrix[rank];
}

double MatrixRectangular::_getValueLocal(int irank) const
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

void MatrixRectangular::_setValueLocal(int irank, double value)
{
  if (! _isRankValid(irank)) return;
  _rectMatrix[irank] = value;
}

void MatrixRectangular::_setValueLocal(int irow, int icol, double value)
{
  if (! _isIndexValid(irow, icol)) return;
  int rank = _getIndexToRank(irow, icol);
  _rectMatrix[rank] = value;
}

void MatrixRectangular::_prodVectorLocal(const double *inv, double *outv) const
{
  matrix_product_safe(getNRows(), getNCols(), 1, _rectMatrix.data(), inv, outv);
}

void MatrixRectangular::_transposeInPlaceLocal()
{
  VectorDouble old;
  old.resize(getNRows() * getNCols());
  matrix_transpose(getNRows(), getNCols(), _rectMatrix.data(), old.data());
  _rectMatrix = old;
  int temp = getNCols();
  _setNCols(getNRows());
  _setNRows(temp);
}

int MatrixRectangular::_getIndexToRankLocal(int irow, int icol) const
{
  return (icol * getNRows() + irow);
}

