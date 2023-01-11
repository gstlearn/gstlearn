/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_old_f.h"

#include "Matrix/MatrixInt.hpp"
#include "Matrix/AMatrix.hpp"
#include "Basic/AException.hpp"
#include "Basic/AStringable.hpp"
#include "csparse_d.h"

MatrixInt::MatrixInt(int nrows, int ncols)
    : AStringable(),
      _nRows(nrows),
      _nCols(ncols),
      _rectMatrix()
{
  _allocate();
}

MatrixInt::MatrixInt(const MatrixInt &r)
  : AStringable(r),
    _nRows(r._nRows),
    _nCols(r._nCols),
    _rectMatrix(r._rectMatrix)
{
}

MatrixInt& MatrixInt::operator= (const MatrixInt &r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _nRows = r._nRows;
    _nCols = r._nCols;
    _rectMatrix = r._rectMatrix;
  }
  return *this;
}

MatrixInt::~MatrixInt()
{
  _deallocate();
}

int MatrixInt::getValue(int irow, int icol) const
{
  _isIndexValid(irow,icol);
  int rank = _getIndexToRank(irow,icol);
  return _rectMatrix[rank];
}

int MatrixInt::getValue(int irank) const
{
  _isRankValid(irank);
  return _rectMatrix[irank];
}

void MatrixInt::setValue(int irank, int value)
{
  _isRankValid(irank);
  _rectMatrix[irank] = value;
}

void MatrixInt::setValue(int irow, int icol, int value)
{
  _isIndexValid(irow, icol);
  int rank = _getIndexToRank(irow, icol);
  _rectMatrix[rank] = value;
}

void MatrixInt::transposeInPlace()
{
  VectorInt old(getNRows() * getNCols());
  matrix_int_transpose(getNRows(), getNCols(), _rectMatrix.data(), old.data());
  _rectMatrix = old;
  int temp = getNCols();
  setNCols(getNRows());
  setNRows(temp);
}

void MatrixInt::fill(int value)
{
  int size = getMatrixSize();
  for (int i = 0; i < size; i++)
    _rectMatrix[i] = value;
}

void MatrixInt::_deallocate()
{

}

void MatrixInt::_allocate()
{
  _rectMatrix.resize(getMatrixSize(),0);
}

int MatrixInt::_getIndexToRank(int irow, int icol) const
{
  int rank = icol * getNRows() + irow;
  return rank;
}

int MatrixInt::getMatrixSize() const
{
  return (getNRows() * getNCols());
}

bool MatrixInt::_isIndexValid(int irow, int icol) const
{
  if (irow < 0 || irow >= getNRows())
  {
    mesArg("Row index invalid",irow,getNRows());
    return false;
  }
  if (icol < 0 || icol >= getNCols())
  {
    mesArg("Column index invalid",icol,getNCols());
    return false;
  }
  return true;
}

bool MatrixInt::_isRankValid(int rank) const
{
  return (rank >= 0 && rank < getMatrixSize());
}

VectorInt MatrixInt::getValues() const
{
  VectorInt vect;
  for (int icol = 0; icol < _nCols; icol++)
    for (int irow = 0; irow < _nRows; irow++)
    {
      int value = getValue(irow,icol);
      vect.push_back(value);
    }
  return vect;
}

VectorInt MatrixInt::getValuesPerRow(int irow) const
{
  VectorInt vect;
  for (int icol = 0; icol < _nCols; icol++)
  {
    int value = getValue(irow,icol);
    vect.push_back(value);
  }
  return vect;
}

VectorInt MatrixInt::getValuesPerColumn(int icol) const
{
  VectorInt vect;
  for (int irow = 0; irow < _nRows; irow++)
  {
    int value = getValue(irow,icol);
    vect.push_back(value);
  }
  return vect;
}

VectorVectorInt MatrixInt::getMatrix() const
{
  VectorVectorInt vect(_nRows);;
  for (int irow = 0; irow < _nRows; irow++)
  {
    vect[irow].resize(_nCols);
    for (int icol = 0; icol < _nCols; icol++)
    {
      int value = getValue(irow,icol);
      vect[irow][icol] = value;
    }
  }
  return vect;
}

/**
 * Filling the matrix with an array of values
 * Note that this array is ALWAYS dimensioned to the total number
 * of elements in the matrix.
 * Kept for compatibility with old code where matrix contents was stored as
 * a int* array
 * @param values Input array (Dimension: nrow * ncol)
 * @param byCol true for Column major; false for Row Major
 */
void MatrixInt::setValues(const VectorInt& values, bool byCol)
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

void MatrixInt::setValuesOldStyle(const int* values, bool byCol)
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

void MatrixInt::reset(int nrows, int ncols)
{
  _isNumbersValid(nrows, ncols);
  _deallocate();
  _nRows = nrows;
  _nCols = ncols;
  _allocate();
}

bool MatrixInt::_isNumbersValid(int nrows, int ncols) const
{
  if (nrows < 0)
  {
    messerr("Argument 'nrows' is not valid");
    return false;
  }
  if (ncols < 0)
  {
    messerr("Argument 'ncols' is not valid");
    return false;
  }
  return true;
}

String MatrixInt::toString(const AStringFormat* /* strfmt*/) const
{
  std::stringstream sstr;

  sstr << "- Number of rows    = " <<  _nRows << std::endl;
  sstr << "- Number of columns = " <<  _nCols << std::endl;

  sstr << toMatrix(String(), VectorString(), VectorString(), true, _nCols, _nRows,
                   getValues());
  return sstr.str();
}

/**
 * Converts a VectorVectorInt into a MatrixInt
 * Note: the input argument is stored by row (if coming from [] specification)
 * @param X Input VectorVectorInt argument
 * @return The returned rectangular matrix
 *
 * @remark: the matrix is transposed implicitly while reading
 */
MatrixInt* MatrixInt::createFromVVD(const VectorVectorInt& X)
{
  int nrow = (int) X.size();
  int ncol = (int) X[0].size();

  MatrixInt* mat = new MatrixInt(nrow, ncol);
  for (int irow = 0; irow < nrow; irow++)
    for (int icol = 0; icol < ncol; icol++)
      mat->setValue(irow, icol, X[irow][icol]);

  return mat;
}

