/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_old_f.h"

#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/AMatrix.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorHelper.hpp"

MatrixRectangular::MatrixRectangular(int nrows, int ncols)
    : AMatrix(nrows, ncols),
      _rectMatrix()
{
  _allocate();
}

MatrixRectangular::MatrixRectangular(const MatrixRectangular &r)
  : AMatrix(r),
    _rectMatrix(r._rectMatrix)
{
}

MatrixRectangular& MatrixRectangular::operator= (const MatrixRectangular &r)
{
  if (this != &r)
  {
    AMatrix::operator=(r);
    _rectMatrix = r._rectMatrix;
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
 * @param X Input VectorVectorDouble argument
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
  if (! _isIndexValid(irow,icol)) return TEST;
  int rank = _getIndexToRank(irow,icol);
  return _rectMatrix[rank];
}

double MatrixRectangular::_getValue(int irank) const
{
  if (! _isRankValid(irank)) return TEST;
  return _rectMatrix[irank];
}

void MatrixRectangular::_setValue(int irank, double value)
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

void MatrixRectangular::_prodVector(const double *inv, double *outv) const
{
  matrix_product_safe(getNRows(), getNCols(), 1, _rectMatrix.data(), inv, outv);
}

void MatrixRectangular::_transposeInPlace()
{
  VectorDouble old;
  old.resize(getNRows() * getNCols());
  matrix_transpose(getNRows(), getNCols(), _rectMatrix.data(), old.data());
  _rectMatrix = old;
  int temp = getNCols();
  _setNCols(getNRows());
  _setNRows(temp);
}

void MatrixRectangular::_setValues(const double *values, bool byCol)
{
  if (byCol)
  {
    int ecr = 0;
    for (int icol = 0; icol < getNCols(); icol++)
      for (int irow = 0; irow < getNRows(); irow++, ecr++)
      {
        setValue(irow, icol, values[ecr]);
      }
  }
  else
  {
    int ecr = 0;
    for (int irow = 0; irow < getNRows(); irow++)
      for (int icol = 0; icol < getNCols(); icol++, ecr++)
      {
        setValue(irow, icol, values[ecr]);
      }
  }
}

/*! Gets a reference to the value at row 'irow' and column 'icol' */
double& MatrixRectangular::_getValueRef(int irow, int icol)
{
  int rank = _getIndexToRank(irow,icol);
  return _rectMatrix[rank];
}

void MatrixRectangular::_deallocate()
{

}

void MatrixRectangular::_allocate()
{
  _rectMatrix.resize(_getMatrixSize(),0.);
}

int MatrixRectangular::_getIndexToRank(int irow, int icol) const
{
  int rank = icol * getNRows() + irow;
  return rank;
}

int MatrixRectangular::_getMatrixSize() const
{
  return (getNRows() * getNCols());
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
  int newNRows = localValidRows.size();
  int newNCols = localValidCols.size();
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
