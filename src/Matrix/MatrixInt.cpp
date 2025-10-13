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
#include "Matrix/MatrixInt.hpp"
#include "Basic/AStringable.hpp"

namespace gstlrn
{
MatrixInt::MatrixInt(Id nrows, Id ncols)
  : AStringable()
  , _nRows(nrows)
  , _nCols(ncols)
  , _rectMatrix()
{
  _allocate();
}

MatrixInt::MatrixInt(const MatrixInt& r)
  : AStringable(r)
  , _nRows(r._nRows)
  , _nCols(r._nCols)
  , _rectMatrix(r._rectMatrix)
{
}

MatrixInt& MatrixInt::operator=(const MatrixInt& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _nRows      = r._nRows;
    _nCols      = r._nCols;
    _rectMatrix = r._rectMatrix;
  }
  return *this;
}

MatrixInt::~MatrixInt()
{
  _deallocate();
}

Id MatrixInt::getValue(Id irow, Id icol) const
{
  if (!_isIndexValid(irow, icol)) return ITEST;
  auto rank = _getIndexToRank(irow, icol);
  return _rectMatrix[rank];
}

Id MatrixInt::getValue(Id irank) const
{
  if (!_isRankValid(irank)) return ITEST;
  return _rectMatrix[irank];
}

Id& MatrixInt::_getValueRef(Id irow, Id icol)
{
  auto rank = _getIndexToRank(irow, icol);
  return _rectMatrix[rank];
}

void MatrixInt::setValueByRank(Id irank, Id value)
{
  _isRankValid(irank);
  _rectMatrix[irank] = value;
}

void MatrixInt::setValue(Id irow, Id icol, Id value)
{
  if (!_isIndexValid(irow, icol)) return;
  auto rank         = _getIndexToRank(irow, icol);
  _rectMatrix[rank] = value;
}

void MatrixInt::transposeInPlace()
{
  VectorInt old(getNRows() * getNCols());
  _transposeInPlace(getNRows(), getNCols(), _rectMatrix.data(), old.data());
  _rectMatrix = old;
  auto temp   = getNCols();
  setNCols(getNRows());
  setNRows(temp);
}

void MatrixInt::fill(Id value)
{
  auto size = getMatrixSize();
  for (Id i = 0; i < size; i++)
    _rectMatrix[i] = value;
}

void MatrixInt::_allocate()
{
  _rectMatrix.resize(getMatrixSize(), 0);
  fill(0);
}

void MatrixInt::_deallocate()
{
  _nRows = 0;
  _nCols = 0;
  _rectMatrix.clear();
}

Id MatrixInt::_getIndexToRank(Id irow, Id icol) const
{
  Id rank = icol * getNRows() + irow;
  return rank;
}

Id MatrixInt::getMatrixSize() const
{
  return (getNRows() * getNCols());
}

bool MatrixInt::_isIndexValid(Id irow, Id icol) const
{
  if (!checkArg("Row index invalid", irow, getNRows())) return false;
  if (!checkArg("Column index invalid", icol, getNCols())) return false;
  return true;
}

bool MatrixInt::_isRankValid(Id rank) const
{
  return (rank >= 0 && rank < getMatrixSize());
}

VectorInt MatrixInt::getValues() const
{
  VectorInt vect;
  for (Id icol = 0; icol < _nCols; icol++)
    for (Id irow = 0; irow < _nRows; irow++)
    {
      auto value = getValue(irow, icol);
      vect.push_back(value);
    }
  return vect;
}

VectorInt MatrixInt::getValuesPerRow(Id irow) const
{
  VectorInt vect;
  for (Id icol = 0; icol < _nCols; icol++)
  {
    auto value = getValue(irow, icol);
    vect.push_back(value);
  }
  return vect;
}

VectorInt MatrixInt::getValuesPerColumn(Id icol) const
{
  VectorInt vect;
  for (Id irow = 0; irow < _nRows; irow++)
  {
    auto value = getValue(irow, icol);
    vect.push_back(value);
  }
  return vect;
}

VectorVectorInt MatrixInt::getMatrix() const
{
  VectorVectorInt vect(_nRows);
  ;
  for (Id irow = 0; irow < _nRows; irow++)
  {
    vect[irow].resize(_nCols);
    for (Id icol = 0; icol < _nCols; icol++)
    {
      auto value       = getValue(irow, icol);
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
 * a Id* array
 * @param values Input array (Dimension: nrow * ncol)
 * @param byCol true for Column major; false for Row Major
 */
void MatrixInt::setValues(const VectorInt& values, bool byCol)
{
  if (byCol)
  {
    Id ecr = 0;
    for (Id icol = 0; icol < getNCols(); icol++)
      for (Id irow = 0; irow < getNRows(); irow++, ecr++)
      {
        setValue(irow, icol, values[ecr]);
      }
  }
  else
  {
    Id ecr = 0;
    for (Id irow = 0; irow < getNRows(); irow++)
      for (Id icol = 0; icol < getNCols(); icol++, ecr++)
      {
        setValue(irow, icol, values[ecr]);
      }
  }
}

void MatrixInt::setValuesOldStyle(const Id* values, bool byCol)
{
  if (byCol)
  {
    Id ecr = 0;
    for (Id icol = 0; icol < getNCols(); icol++)
      for (Id irow = 0; irow < getNRows(); irow++, ecr++)
      {
        setValue(irow, icol, values[ecr]);
      }
  }
  else
  {
    Id ecr = 0;
    for (Id irow = 0; irow < getNRows(); irow++)
      for (Id icol = 0; icol < getNCols(); icol++, ecr++)
      {
        setValue(irow, icol, values[ecr]);
      }
  }
}

void MatrixInt::reset(Id nrows, Id ncols)
{
  _isNumbersValid(nrows, ncols);
  _deallocate();
  _nRows = nrows;
  _nCols = ncols;
  _allocate();
}

void MatrixInt::resetFromArray(Id nrows, Id ncols, const Id* tab, bool byCol)
{
  reset(nrows, ncols);

  Id lec = 0;
  if (byCol)
  {
    for (Id icol = 0; icol < ncols; icol++)
      for (Id irow = 0; irow < nrows; irow++)
        setValue(irow, icol, tab[lec++]);
  }
  else
  {
    for (Id irow = 0; irow < nrows; irow++)
      for (Id icol = 0; icol < ncols; icol++)
        setValue(irow, icol, tab[lec++]);
  }
}

bool MatrixInt::_isNumbersValid(Id nrows, Id ncols)
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

  sstr << "- Number of rows    = " << _nRows << std::endl;
  sstr << "- Number of columns = " << _nCols << std::endl;

  sstr << toMatrix(String(), VectorString(), VectorString(), true, _nRows, _nCols,
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
MatrixInt* MatrixInt::createFromVVI(const VectorVectorInt& X)
{
  Id nrow = static_cast<Id>(X.size());
  Id ncol = static_cast<Id>(X[0].size());

  auto* mat = new MatrixInt(nrow, ncol);
  for (Id irow = 0; irow < nrow; irow++)
    for (Id icol = 0; icol < ncol; icol++)
      mat->setValue(irow, icol, X[irow][icol]);

  return mat;
}

MatrixInt* MatrixInt::createFromVI(const VectorInt& X,
                                   Id nrow,
                                   Id ncol,
                                   bool byCol)
{
  if (nrow * ncol != static_cast<Id>(X.size()))
  {
    messerr("Inconsistency between arguments 'nrow'(%d) and 'ncol'(%d)", nrow, ncol);
    messerr("and the dimension of the input Vector (%d)", static_cast<Id>(X.size()));
  }
  auto* mat = new MatrixInt(nrow, ncol);

  Id lec = 0;
  if (byCol)
  {
    for (Id irow = 0; irow < nrow; irow++)
      for (Id icol = 0; icol < ncol; icol++)
      {
        mat->setValue(irow, icol, X[lec++]);
      }
  }
  else
  {
    for (Id icol = 0; icol < ncol; icol++)
      for (Id irow = 0; irow < nrow; irow++)
      {
        mat->setValue(irow, icol, X[lec++]);
      }
  }
  return mat;
}

/*****************************************************************************/
/*!
 **  Transpose a (square or rectangular) matrix
 **
 ** \param[in]  n1 matrix dimension
 ** \param[in]  n2 matrix dimension
 ** \param[in]  v1 rectangular matrix (n1,n2)
 **
 ** \param[out] w1 rectangular matrix (n2,n1)
 **
 ** \remark  The matrix w1[] may NOT coincide with v1[]
 **
 *****************************************************************************/
void MatrixInt::_transposeInPlace(Id n1, Id n2, const Id* v1, Id* w1)
{
#define SQ(i, j, neq) ((j) * neq + (i))
#define V1(i, j)      v1[SQ(i, j, n1)]
  Id ecr = 0;
  for (Id i1 = 0; i1 < n1; i1++)
    for (Id i2 = 0; i2 < n2; i2++)
      w1[ecr++] = V1(i1, i2);
}

} // namespace gstlrn
