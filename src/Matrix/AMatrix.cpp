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
#include "Matrix/AMatrix.hpp"
#include "Basic/AException.hpp"
#include "Basic/Law.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Matrix/NF_Triplet.hpp"
#include "geoslib_define.h"

#include <iostream>

namespace gstlrn
{

static I32 _globalMultiThread = 0;
static bool _flagMatrixCheck  = false;

AMatrix::AMatrix(Id nrow, Id ncol)
  : AStringable()
  , _nRows(nrow)
  , _nCols(ncol)
  , _nullTerm(0.)
{
}

AMatrix::AMatrix(const AMatrix& m)
  : AStringable(m)
  , _nRows(m._nRows)
  , _nCols(m._nCols)
  , _nullTerm(m._nullTerm)
{
}

AMatrix& AMatrix::operator=(const AMatrix& m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _nRows    = m._nRows;
    _nCols    = m._nCols;
    _nullTerm = m._nullTerm;
  }
  return *this;
}

AMatrix::~AMatrix()
{
}

void AMatrix::reset(Id nrows, Id ncols)
{
  // Check if numbers are valid
  if (!_isNumbersValid(nrows, ncols)) return;

  // Reset memory
  _deallocate(); // virtual call
  _setNRows(nrows);
  _setNCols(ncols);
  _allocate(); // virtual call
}

/**
 * @brief Reset the matrix to new dimensions and fill with a new value
 *
 * @param nrows New number of rows
 * @param ncols New number of columns
 * @param value The new value used to fill the matrix
 */
void AMatrix::resetFromValue(Id nrows, Id ncols, double value)
{
  reset(nrows, ncols);
  fill(value);
}

/**
 * @brief Reset the matrix from an array of double values
 *
 * @param nrows New number of rows
 * @param ncols New number of columns
 * @param tab The array of values
 * @param byCol True if values are column-major in the array
 */
void AMatrix::resetFromArray(Id nrows, Id ncols, const double* tab, bool byCol)
{
  reset(nrows, ncols);
  if (byCol)
  {
    Id lec = 0;
    for (Id icol = 0; icol < ncols; icol++)
      for (Id irow = 0; irow < nrows; irow++)
        setValue(irow, icol, tab[lec++]);
  }
  else
  {
    Id lec = 0;
    for (Id irow = 0; irow < nrows; irow++)
      for (Id icol = 0; icol < ncols; icol++)
        setValue(irow, icol, tab[lec++]);
  }
}

/**
 * @brief Reset the matrix from a vector of double values
 *
 * @param nrows New number of rows
 * @param ncols New number of columns
 * @param tab The vector of values
 * @param byCol True if values are column-major in the vector
 */
void AMatrix::resetFromVD(Id nrows, Id ncols, const VectorDouble& tab, bool byCol)
{
  resetFromArray(nrows, ncols, tab.data(), byCol);
}

/**
 * @brief Reset the matrix from an array of double values
 *
 * @param tab The array of values
 * @param byCol True if values are column-major in the array
 */
void AMatrix::resetFromVVD(const VectorVectorDouble& tab, bool byCol)
{
  if (!byCol)
  {
    Id nrows = static_cast<Id>(tab.size());
    Id ncols = static_cast<Id>(tab[0].size());
    reset(nrows, ncols);
    for (Id icol = 0; icol < ncols; icol++)
      for (Id irow = 0; irow < nrows; irow++)
        setValue(irow, icol, tab[irow][icol]);
  }
  else
  {
    Id ncols = static_cast<Id>(tab.size());
    Id nrows = static_cast<Id>(tab[0].size());
    reset(nrows, ncols);
    for (Id icol = 0; icol < ncols; icol++)
      for (Id irow = 0; irow < nrows; irow++)
        setValue(irow, icol, tab[icol][irow]);
  }
}

/**
 * Indicate if the given matrce is a square matrix
 *
 * @param printWhyNot Print the message is the answer is false
 * @return true if the matrix is square
 */
bool AMatrix::isSquare(bool printWhyNot) const
{
  if (empty()) return false;
  if (_nRows != _nCols)
  {
    if (printWhyNot)
      messerr("The number of rows (%d) should match the number of columns (%d)",
              _nRows, _nCols);
    return false;
  }
  return true;
}

/**
 * Indicate if the given indices are valid for the current matrix size
 *
 * @param irow Row index
 * @param icol Column index
 * @param printWhyNot Print the message is the answer is false
 * @return true if indices are valid for the current matrix size
 */
bool AMatrix::isValid(Id irow, Id icol, bool printWhyNot) const
{
  if (!_flagMatrixCheck) return true;
  if (irow < 0 || irow >= getNRows())
  {
    if (printWhyNot)
      messerr("Argument 'irow' invalid: it should lie in [0;%d[",
              irow, getNRows());
    return false;
  }
  if (icol < 0 || icol >= getNCols())
  {
    if (printWhyNot)
      messerr("Argument 'icol' invalid: it should lie in [0;%d[",
              icol, getNCols());
    return false;
  }
  return true;
}

/**
 * Check that Matrix 'm' share the same dimensions as current one
 *
 * @param m Matrix to be compared to the current Matrix
 * @return true if 'm' has same dimensions as the current Matrix
 *
 * @remark:  message is issued if dimensions are different
 */
bool AMatrix::isSameSize(const AMatrix& m) const
{
  if (_nRows == m.getNRows() && _nCols == m.getNCols()) return true;

  messerr("Dimensions of matrices do not match");
  messerr("- current matrix (%d x %d)", getNRows(), getNCols());
  messerr("- tested matrix  (%d x %d)", m.getNRows(), m.getNCols());
  return false;
}

/**
 * Check that Matrix 'm' is equal to the current Matrix
 *
 * @param m Matrix to be compared to the current Matrix
 * @param eps Epsilon for double equality comparison
 * @param printWhyNot Print the message is the answer is false
 * @return true if 'm'  is equal to the current Matrix
 */
bool AMatrix::isSame(const AMatrix& m, double eps, bool printWhyNot)
{
  if (!isSameSize(m)) return false;

  auto ncols = getNCols();
  auto nrows = getNRows();
  for (Id icol = 0; icol < ncols; icol++)
    for (Id irow = 0; irow < nrows; irow++)
    {
      double v1 = getValue(irow, icol);
      double v2 = m.getValue(irow, icol);
      if (ABS(v1 - v2) > eps)
      {
        if (printWhyNot)
        {
          messerr("Element (%d;%d) are different between:\n", irow, icol);
          messerr("- First matrix");
          m.display();
          messerr("- Second matrix");
          display();
        }
        return false;
      }
    }
  return true;
}

/**
 * Indicate if the current matrix is symmetric
 *
 * @param eps Epsilon for double equality comparison
 * @param printWhyNot Print the message is the answer is false
 * @return true if the current matrix is symmetric
 */
bool AMatrix::isSymmetric(double eps, bool printWhyNot) const
{
  if (empty() || !isSquare()) return false;

  for (Id irow = 0; irow < _nRows; irow++)
    for (Id icol = 0; icol < _nCols; icol++)
    {
      if (ABS(getValue(irow, icol) - getValue(icol, irow)) > eps)
      {
        if (printWhyNot)
          messerr("Elements (%d;%d)=%lf and (%d;%d)=%kf should be equal",
                  irow, icol, getValue(irow, icol),
                  icol, irow, getValue(icol, irow));
        return false;
      }
    }
  return true;
}

/**
 * Indicate if the current matrix is the Identity
 *
 * @param printWhyNot Print the message is the answer is false
 * @return true if the current matrix is the Identity
 */
bool AMatrix::isIdentity(bool printWhyNot) const
{
  for (Id irow = 0; irow < getNRows(); irow++)
    for (Id icol = 0; icol < getNCols(); icol++)
    {
      double refval = (irow == icol) ? 1. : 0.;
      if (ABS(getValue(irow, icol) - refval) > EPSILON10)
      {
        if (printWhyNot)
          messerr("The term (%d,%d) should be equal to %lf (%lf)", irow + 1,
                  icol + 1, refval, getValue(irow, icol));
        return false;
      }
    }
  return true;
}

void AMatrix::transposeInPlace()
{
  _transposeInPlace();
}

AMatrix* AMatrix::transpose() const
{
  AMatrix* mat = dynamic_cast<AMatrix*>(clone());
  mat->transposeInPlace();
  return mat;
}

/**
 * Filling the matrix with an array of values
 * Note that this array is ALWAYS dimensioned to the total number
 * of elements in the matrix.
 * Kept for compatibility with old code where matrix contents was stored as
 * a double* array
 * @param values Input array (Dimension: nrow * ncol)
 * @param byCol true for Column major; false for Row Major
 */
#ifndef SWIG
void AMatrix::_setValues(const double* values, bool byCol)
{
  Id ecr = 0;
  if (byCol)
  {
    for (Id icol = 0; icol < getNCols(); icol++)
      for (Id irow = 0; irow < getNRows(); irow++, ecr++)
        setValue(irow, icol, values[ecr]);
  }
  else
  {
    for (Id irow = 0; irow < getNRows(); irow++)
      for (Id icol = 0; icol < getNCols(); icol++, ecr++)
        setValue(irow, icol, values[ecr]);
  }
}
#endif

void AMatrix::fillRandom(Id seed, double zeroPercent)
{
  law_set_random_seed(seed);

  double value = 0.;
  for (Id irow = 0; irow < _nRows; irow++)
    for (Id icol = 0; icol < _nCols; icol++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;
      if (law_uniform(0., 1.) < zeroPercent)
        value = 0.;
      else
        value = law_gaussian();
      setValue(irow, icol, value);
    }
}

/**
 * Filling the matrix with an array of values
 * Note that this array is ALWAYS dimensioned to the total number
 * of elements in the matrix.
 * Kept for compatibility with old code where matrix contents was stored as
 * a VectorDouble
 * @param values
 * @param byCol true for Column major; false for Row Major
 */
void AMatrix::setValues(const VectorDouble& values, bool byCol)
{
  if (static_cast<Id>(values.size()) != size())
  {
    messerr("Inconsistency between 'values' and Matrix Dimension");
    messerr("Operation cancelled");
    return;
  }
  _setValues(values.data(), byCol);
}

void AMatrix::setIdentity(double value)
{
  for (Id icol = 0; icol < _nCols; icol++)
    for (Id irow = 0; irow < _nRows; irow++)
      setValue(irow, icol, value * (irow == icol));
}

VectorDouble AMatrix::prodMatVec(const VectorDouble& x, bool transpose) const
{
  if (_flagMatrixCheck &&
      !_isMatrixCompatible("AMatrix::prodMatVec",
                           this, 0, transpose,
                           nullptr, static_cast<Id>(x.size()), false)) return VectorDouble();
  Id size = (transpose) ? _nCols : _nRows;
  VectorDouble y(size, 0.);
  _addProdMatVecInPlacePtr(x, y, transpose);
  return y;
}

/**
 * Returns 'y' = 'this' %*% 'x'
 * @param x Input vector
 * @param y Output vector
 * @param transpose True if the matrix 'this' must be transposed
 */
void AMatrix::prodMatVecInPlace(const VectorDouble& x,
                                VectorDouble& y,
                                bool transpose) const
{
  if (_flagMatrixCheck &&
      !_isMatrixCompatible("AMatrix::prodMatVecInPlace",
                           this, 0, transpose,
                           nullptr, static_cast<Id>(x.size()), false)) return;
  Id size = (transpose) ? _nCols : _nRows;
  y.fill(0., size);
  _addProdMatVecInPlacePtr(x, y, transpose);
}

void AMatrix::prodMatVecInPlaceC(const constvect x,
                                 vect y,
                                 bool transpose) const
{
  if (_flagMatrixCheck &&
      !_isMatrixCompatible("AMatrix::prodMatVecInPlaceC",
                           this, 0, transpose,
                           nullptr, static_cast<Id>(x.size()), false)) return;
  Id size = (transpose) ? _nCols : _nRows;
  std::fill(y.begin(), y.begin() + size, 0.0);
  _addProdMatVecInPlacePtr(x, y, transpose);
}

void AMatrix::addProdMatVecInPlaceC(const constvect x,
                                    vect y,
                                    bool transpose) const
{
  if (_flagMatrixCheck &&
      !_isMatrixCompatible("AMatrix::addProdMatVecInPlaceC",
                           this, 0, transpose,
                           nullptr, static_cast<Id>(x.size()), false)) return;
  _addProdMatVecInPlacePtr(x, y, transpose);
}

VectorDouble AMatrix::prodVecMat(const VectorDouble& x, bool transpose) const
{
  if (_flagMatrixCheck &&
      !_isMatrixCompatible("AMatrix::prodVecMat",
                           nullptr, static_cast<Id>(x.size()), true,
                           this, 0, transpose)) return VectorDouble();
  Id size = (transpose) ? _nRows : _nCols;
  VectorDouble y(size, 0.);
  _addProdVecMatInPlacePtr(x, y, transpose);
  return y;
}

void AMatrix::prodVecMatInPlace(const VectorDouble& x, VectorDouble& y, bool transpose) const
{
  if (_flagMatrixCheck &&
      !_isMatrixCompatible("AMatrix::prodVecMat",
                           nullptr, static_cast<Id>(x.size()), true,
                           this, 0, transpose)) return;
  Id size = (transpose) ? _nRows : _nCols;
  y.fill(0., size);
  _addProdVecMatInPlacePtr(x, y, transpose);
}

void AMatrix::prodVecMatInPlaceC(const constvect x, vect y, bool transpose) const
{
  if (_flagMatrixCheck &&
      !_isMatrixCompatible("AMatrix::prodVecMatInPlaceC",
                           this, 0, transpose,
                           nullptr, static_cast<Id>(x.size()), false)) return;
  Id size = (transpose) ? _nCols : _nRows;
  std::fill(y.begin(), y.begin() + size, 0.0);
  _addProdVecMatInPlacePtr(x, y, transpose);
}

void AMatrix::addProdVecMatInPlaceC(const constvect x,
                                    vect y,
                                    bool transpose) const
{
  if (_flagMatrixCheck &&
      !_isMatrixCompatible("AMatrix::addProdVecMatInPlaceC",
                           this, 0, transpose,
                           nullptr, static_cast<Id>(x.size()), false)) return;
  _addProdVecMatInPlacePtr(x, y, transpose);
}

bool AMatrix::_matrixNeedToReset(Id nrows, Id ncols)
{
  return nrows != getNRows() || ncols != getNCols() || _needToReset(nrows, ncols);
}

bool AMatrix::_needToReset(Id nrows, Id ncols)
{
  DECLARE_UNUSED(nrows, ncols)
  return false;
}
/**
 * @brief Resize the matrix to new dimensions
 *        (this method doesn't change the storage type)
 *
 * @param nrows New number of rows
 * @param ncols New number of columns
 */
void AMatrix::resize(Id nrows, Id ncols)
{
  // Check if nothing is to be done
  if (!_matrixNeedToReset(nrows, ncols))
    return;

  // Reset the sizes (clear values)
  reset(nrows, ncols);
}

/**
 * Add the matrix 'y' to the current Matrix
 * @param y Matrix to be added
 * @param cx Multiplicative parameter for this
 * @param cy Multiplicative parameter for y
 */
void AMatrix::addMat(const AMatrix& y, double cx, double cy)
{
  if (!isSameSize(y)) return;

  for (Id irow = 0; irow < _nRows; irow++)
    for (Id icol = 0; icol < _nCols; icol++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;
      setValue(irow, icol, cx * getValue(irow, icol) + cy * y.getValue(irow, icol));
    }
}

/**
 * Store the product of 'x'(or 't(x)') by 'y' (or 't(y') in this
 * @param x First Matrix
 * @param y Second matrix
 * @param transposeX True if first matrix must be transposed
 * @param transposeY True if second matrix must be transposed
 */
void AMatrix::prodMatMatInPlace(const AMatrix* x,
                                const AMatrix* y,
                                bool transposeX,
                                bool transposeY)
{
  if (_flagMatrixCheck &&
      !_isMatrixCompatible("AMatrix::prodMatMatInPlace",
                           x, 0, transposeX,
                           y, 0, transposeY)) return;

  Id ni1 = (transposeX) ? x->getNCols() : x->getNRows();
  Id nm1 = (transposeX) ? x->getNRows() : x->getNCols();
  Id ni2 = (transposeY) ? y->getNRows() : y->getNCols();
  Id nm2 = (transposeY) ? y->getNCols() : y->getNRows();

  if (nm1 != ni2)
  {
    messerr("Matrices 'x' and 'y' should have matching dimensions");
    return;
  }

  for (Id irow = 0; irow < ni1; irow++)
  {
    for (Id icol = 0; icol < nm2; icol++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;

      double value = 0.;
      for (Id k = 0; k < nm1; k++)
      {
        double v1 = (transposeX) ? x->getValue(k, irow) : x->getValue(irow, k);
        double v2 = (transposeY) ? y->getValue(icol, k) : y->getValue(k, icol);
        value += v1 * v2;
      }
      setValue(irow, icol, value);
    }
  }
}

/**
 * Product 't(A)' %*% 'M' %*% 'A' or 'A' %*% 'M' %*% 't(A)' stored in 'this'
 * @param a Matrix A
 * @param m Matrix M
 * @param transpose True for first implementation, False for the second
 */
void AMatrix::prodNormMatMatInPlace(const AMatrix* a,
                                    const AMatrix* m,
                                    bool transpose)
{
  if (_flagMatrixCheck &&
      !_isMatrixCompatible("AMatrix::prodNormMatMatInPlace",
                           a, 0, transpose,
                           m, 0, false,
                           a, 0, !transpose)) return;

  Id n1 = (transpose) ? a->getNCols() : a->getNRows();
  Id n2 = (transpose) ? a->getNRows() : a->getNCols();
  for (Id i = 0; i < n1; i++)
    for (Id j = 0; j < n1; j++)
    {
      if (!_isPhysicallyPresent(i, j)) continue;
      double value = 0;
      for (Id k = 0; k < n2; k++)
        for (Id l = 0; l < n2; l++)
        {
          double a_ik = (transpose) ? a->getValue(k, i) : a->getValue(i, k);
          double a_lj = (transpose) ? a->getValue(l, j) : a->getValue(j, l);
          value += a_ik * m->getValue(k, l) * a_lj;
        }
      setValue(i, j, value);
    }
}

void AMatrix::prodNormMatVecInPlace(const AMatrix* a, const VectorDouble& vec, bool transpose)
{
  if (_flagMatrixCheck &&
      !_isMatrixCompatible("AMatrix::prodNormMatVecInPlace",
                           a, 0, transpose,
                           a, 0, !transpose)) return;

  Id n1 = (transpose) ? a->getNCols() : a->getNRows();
  Id n2 = (transpose) ? a->getNRows() : a->getNCols();
  for (Id i = 0; i < n1; i++)
    for (Id j = 0; j < n1; j++)
    {
      if (!_isPhysicallyPresent(i, j)) continue;
      double value = 0;
      for (Id k = 0; k < n2; k++)
      {
        double a_ik     = (transpose) ? a->getValue(k, i) : a->getValue(i, k);
        double a_kj     = (transpose) ? a->getValue(j, k) : a->getValue(k, j);
        double vecvalue = vec[k];
        value += a_ik * vecvalue * a_kj;
      }
      setValue(i, j, value);
    }
}

void AMatrix::prodNormMatInPlace(const AMatrix* a, bool transpose)
{
  if (_flagMatrixCheck &&
      !_isMatrixCompatible("AMatrix::prodNormMatInPlace",
                           a, 0, transpose,
                           a, 0, !transpose)) return;

  Id n1 = (transpose) ? a->getNCols() : a->getNRows();
  Id n2 = (transpose) ? a->getNRows() : a->getNCols();
  for (Id i = 0; i < n1; i++)
    for (Id j = 0; j < n1; j++)
    {
      if (!_isPhysicallyPresent(i, j)) continue;
      double value = 0;
      for (Id k = 0; k < n2; k++)
      {
        double a_ik = (transpose) ? a->getValue(k, i) : a->getValue(i, k);
        double a_kj = (transpose) ? a->getValue(j, k) : a->getValue(k, j);
        value += a_ik * a_kj;
      }
      setValue(i, j, value);
    }
}

double AMatrix::prodVecMatVec(const VectorDouble& x, const VectorDouble& y) const
{
  if (_flagMatrixCheck &&
      !_isMatrixCompatible("AMatrix::quadraticMatrix",
                           nullptr, static_cast<Id>(x.size()), true,
                           this, 0, false,
                           nullptr, static_cast<Id>(y.size()), false)) return TEST;

  VectorDouble left(_nRows);
  prodMatVecInPlace(y, left, false);
  return VH::innerProduct(x, left);
}

Id AMatrix::invert()
{
  if (!isSquare())
  {
    messerr("'invert' method is restricted to Square Matrices");
    return 1;
  }
  return _invert();
}

Id AMatrix::solve(const VectorDouble& b, VectorDouble& x) const
{
  if (!isSquare())
  {
    messerr("'solve' method is limited to Square Matrices");
    return 1;
  }
  if (static_cast<Id>(b.size()) != _nRows || static_cast<Id>(x.size()) != _nRows)
  {
    messerr("b' and 'x' should have the same dimension as the Matrix");
    return 1;
  }
  return _solve(b, x);
}

String AMatrix::toString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt);
  std::stringstream sstr;

  sstr << "- Number of rows    = " << _nRows << std::endl;
  sstr << "- Number of columns = " << _nCols << std::endl;

  bool flagSkipZero = false;
  if (isSparse())
  {
    sstr << "- Sparse Format" << std::endl;
    flagSkipZero = true;
  }
  sstr << toMatrix(String(), VectorString(), VectorString(), true, _nRows, _nCols,
                   getValues(), false, flagSkipZero);

  return sstr.str();
}

void AMatrix::clear()
{
  _clear();
}

void AMatrix::_clear()
{
  _setNRows(0);
  _setNCols(0);
}

bool AMatrix::_isNumbersValid(Id nrows, Id ncols) const
{
  if (!_flagMatrixCheck) return true;
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

bool AMatrix::_isRowValid(Id irow) const
{
  if (!_flagMatrixCheck) return true;
  return checkArg("Row index invalid", irow, getNRows());
}

bool AMatrix::_isColumnValid(Id icol) const
{
  if (!_flagMatrixCheck) return true;
  return checkArg("Column index invalid", icol, getNCols());
}

bool AMatrix::_isIndexValid(Id irow, Id icol) const
{
  if (!_flagMatrixCheck) return true;
  if (!_isRowValid(irow)) return false;
  if (!_isColumnValid(icol)) return false;
  return true;
}

bool AMatrix::_isRowVectorConsistent(const VectorDouble& tab) const
{
  if (static_cast<Id>(tab.size()) != getNRows())
  {
    messerr("Argument vector size should match the number of rows");
    return false;
  }
  return true;
}

bool AMatrix::_isColVectorConsistent(const VectorDouble& tab) const
{
  if (static_cast<Id>(tab.size()) != getNCols())
  {
    messerr("Argument vector size should match the number of columns");
    return false;
  }
  return true;
}

bool AMatrix::_isVectorSizeConsistent(const VectorDouble& tab) const
{
  auto nrows = getNRows();
  auto ncols = getNCols();
  if (static_cast<Id>(tab.size()) != nrows * ncols)
  {
    messerr("The argument 'tab'(%d) does not have correct dimension (%d)",
            static_cast<Id>(tab.size()), nrows * ncols);
    return false;
  }
  return true;
}

bool AMatrix::_isColumnSizeConsistent(const VectorDouble& tab) const
{
  auto nrows = getNRows();
  if (static_cast<Id>(tab.size()) != nrows)
  {
    messerr("The argument 'tab'(%d) does not have correct dimension (%d)",
            static_cast<Id>(tab.size()), nrows);
    return false;
  }
  return true;
}

bool AMatrix::_isRowSizeConsistent(const VectorDouble& tab) const
{
  auto ncols = getNCols();
  if (static_cast<Id>(tab.size()) != ncols)
  {
    messerr("The argument 'tab'(%d) does not have correct dimension (%d)",
            static_cast<Id>(tab.size()), ncols);
    return false;
  }
  return true;
}

bool AMatrix::_isRankValid(Id rank) const
{
  if (!_flagMatrixCheck) return true;
  return (rank >= 0 && rank < _getMatrixPhysicalSize());
}

void AMatrix::dumpElements(const String& title, Id ifrom, Id ito) const
{
  mestitle(1, "%s", title.c_str());
  for (Id rank = ifrom; rank < ito; rank++)
  {
    if (_isRankValid(rank))
      message("Element %d = %lf\n", rank, _getValueByRank(rank));
  }
}

void AMatrix::dumpStatistics(const String& title) const
{
  message("%s : %d rows and %d columns\n", title.c_str(), _nRows, _nCols);
}

bool AMatrix::_identifyRowAndCol(const AMatrix* mat,
                                 Id vsize,
                                 bool transpose,
                                 Id* nrow,
                                 Id* ncol)
{
  Id nr;
  Id nc;
  if (mat != nullptr)
  {
    if (vsize > 0)
    {
      messerr("Both Matrix and VectorDouble are defined. Only one should be used");
      return false;
    }
    nr = mat->getNRows();
    nc = mat->getNCols();
  }
  else if (vsize > 0)
  {
    nr = vsize;
    nc = 1;
  }
  else
  {
    // Both 'mat' and 'vd' are undefined
    return false;
  }

  if (transpose)
  {
    *nrow = nc;
    *ncol = nr;
  }
  else
  {
    *nrow = nr;
    *ncol = nc;
  }
  return true;
}

/**
 * Check that a set of matrices (or vectors) has the correct linkage
 * The linkage is defined as follows:
 *        mat1 (ot vd1) * mat2 (or vd2) * mat3 (or vd3)
 * @param name        Name of the calling function
 * @param mat1        First matrix (optional)
 * @param vsize1      Dimension of the First  VectorDouble (optional)
 * @param transpose1  True if the first matrix must be transposed
 * @param mat2        Second matrix (optional)
 * @param vsize2      Dimension of the Second VectorDouble (optional)
 * @param transpose2  True if the second matrix must be transposed
 * @param mat3        Third matrix (optional)
 * @param vsize3      Dimension of the Third VectorDouble (optional)
 * @param transpose3  True if the third matrix must be transposed
 * @return True if the linkage is correct
 *
 * @remark A vector must be defined with its 'ncol' set to 1
 */
bool AMatrix::_isMatrixCompatible(const String& name,
                                  const AMatrix* mat1,
                                  Id vsize1,
                                  bool transpose1,
                                  const AMatrix* mat2,
                                  Id vsize2,
                                  bool transpose2,
                                  const AMatrix* mat3,
                                  Id vsize3,
                                  bool transpose3)
{
  // First element
  Id nrow1;
  Id ncol1;
  if (!_identifyRowAndCol(mat1, vsize1, transpose1, &nrow1, &ncol1)) return false;

  // Second element
  Id nrow2;
  Id ncol2;
  if (_identifyRowAndCol(mat2, vsize2, transpose2, &nrow2, &ncol2))
  {
    if (nrow2 != ncol1)
    {
      messerr("Linkage error called from %s: Level = 1 - Nrow = %d - Ncol = %d",
              name.c_str(), nrow2, ncol1);
      return false;
    }

    // Third element
    Id nrow3;
    Id ncol3;
    if (_identifyRowAndCol(mat3, vsize3, transpose3, &nrow3, &ncol3))
    {
      if (nrow3 != ncol2)
      {
        messerr("Linkage error called from %s: Level = 2 - Nrow = %d - Ncol = %d",
                name.c_str(), nrow3, ncol2);
        return false;
      }
    }
  }

  return true;
}

/**
 * From a matrix of any type, creates the three vectors of the triplet
 * (specific format for creating efficiently a Sparse matrix)
 * It only takes the only non-zero elements of the matrix
 */
NF_Triplet AMatrix::getMatrixToTriplet(Id shiftRow, Id shiftCol) const
{
  NF_Triplet NF_T;

  for (Id icol = 0; icol < _nCols; icol++)
    for (Id irow = 0; irow < _nRows; irow++)
    {
      if (!isValid(irow, icol)) continue;
      double value = getValue(irow, icol);
      if (isZero(value)) continue;
      NF_T.add(irow + shiftRow, icol + shiftCol, value);
    }
  return NF_T;
}

VectorDouble AMatrix::getValues(bool byCol) const
{
  VectorDouble vect(_nCols * _nRows);
  VectorDouble::iterator itvect(vect.begin());

  if (byCol)
  {
    for (Id icol = 0; icol < _nCols; icol++)
      for (Id irow = 0; irow < _nRows; irow++)
      {
        (*itvect) = getValue(irow, icol);
        itvect++;
      }
  }
  else
  {
    for (Id irow = 0; irow < _nRows; irow++)
      for (Id icol = 0; icol < _nCols; icol++)
      {
        (*itvect) = getValue(irow, icol);
        itvect++;
      }
  }
  return vect;
}

double AMatrix::compare(const AMatrix& mat) const
{
  if (mat.getNRows() != _nRows || mat.getNCols() != _nCols)
  {
    messerr("We can only compare two matrices with same dimensions");
    return TEST;
  }

  double diff = 0.;
  for (Id icol = 0; icol < _nCols; icol++)
    for (Id irow = 0; irow < _nRows; irow++)
    {
      double v1 = (isValid(irow, icol)) ? getValue(irow, icol) : 0.;
      double v2 = (mat.isValid(irow, icol)) ? mat.getValue(irow, icol) : 0.;
      diff += ABS(v1 - v2);
    }
  return diff;
}

const VectorDouble& AMatrix::getDiagonal(Id shift) const
{
  _diagonal.clear();
  if (!isSquare())
  {
    messerr("This function is only valid for Square matrices");
    return this->_diagonal;
  }

  _diagonal.reserve(getNRows());
  for (Id rank = 0; rank < getNRows(); rank++)
  {
    Id irow = rank;
    if (shift < 0) irow += shift;
    if (irow < 0 || irow >= getNRows()) continue;
    Id icol = rank;
    if (shift > 0) icol += shift;
    if (icol < 0 || icol >= getNCols()) continue;
    _diagonal.push_back(getValue(irow, icol));
  }
  return _diagonal;
}

/*! Extract a Row */
VectorDouble AMatrix::getRow(Id irow) const
{
  VectorDouble vect;
  if (!checkArg("Incorrect argument 'irow'", irow, getNRows())) return vect;

  for (Id icol = 0; icol < getNCols(); icol++)
    vect.push_back(getValue(irow, icol));
  return vect;
}

/*! Extract a Column */
VectorDouble AMatrix::getColumn(Id icol) const
{
  if (icol < 0 || icol >= getNCols())
    my_throw("Incorrect argument 'icol'");

  VectorDouble vect;
  for (Id irow = 0; irow < getNRows(); irow++)
    vect.push_back(getValue(irow, icol));
  return vect;
}

VectorDouble AMatrix::getColumnByRowRange(Id icol, Id rowFrom, Id rowTo) const
{
  if (icol < 0 || icol >= getNCols())
    my_throw("Incorrect argument 'icol'");

  VectorDouble vect;
  for (Id irow = rowFrom; irow < rowTo; irow++)
    vect.push_back(getValue(irow, icol));
  return vect;
}

/*! Checks if a Column is valid (contains a non TEST value) */
bool AMatrix::isColumnDefined(Id icol) const
{
  if (icol < 0 || icol >= getNCols())
    my_throw("Incorrect argument 'icol'");

  for (Id irow = 0; irow < getNRows(); irow++)
  {
    if (!FFFF(getValue(irow, icol))) return true;
  }
  return false;
}

/*! Checks if a Row is valid (contains a non TEST value) */
bool AMatrix::isRowDefined(Id irow) const
{
  if (irow < 0 || irow >= getNRows())
    my_throw("Incorrect argument 'irow'");

  for (Id icol = 0; icol < getNCols(); icol++)
  {
    if (!FFFF(getValue(irow, icol))) return true;
  }
  return false;
}

/*! Define the number of defined columns */
Id AMatrix::getNColDefined() const
{
  Id ncol = 0;
  for (Id icol = 0; icol < getNCols(); icol++)
  {
    if (isColumnDefined(icol)) ncol++;
  }
  return ncol;
}

/*! Define the number of defined rows */
Id AMatrix::getNRowDefined() const
{
  Id nrow = 0;
  for (Id irow = 0; irow < getNRows(); irow++)
  {
    if (isRowDefined(irow)) nrow++;
  }
  return nrow;
}

void AMatrix::addValue(Id irow, Id icol, double value)
{
  double oldval = getValue(irow, icol);
  if (FFFF(oldval)) return;
  setValue(irow, icol, oldval + value);
}

double AMatrix::getMeanByColumn(Id icol) const
{
  double cumul = 0.;
  double count = 0.;
  for (Id irow = 0; irow < getNRows(); irow++)
  {
    double value = getValue(irow, icol);
    if (FFFF(value)) continue;
    cumul += value;
    count += 1.;
  }

  if (count <= 0.) return TEST;
  return cumul / count;
}

double AMatrix::getMinimum() const
{
  double minimum = MAXIMUM_BIG;
  for (Id icol = 0; icol < getNCols(); icol++)
    for (Id irow = 0; irow < getNRows(); irow++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;
      double value = getValue(irow, icol);
      if (FFFF(value)) continue;
      if (value < minimum) minimum = value;
    }
  if (isEqual(minimum, MAXIMUM_BIG)) minimum = TEST;
  return minimum;
}

double AMatrix::getMaximum() const
{
  double maximum = MINIMUM_BIG;
  for (Id icol = 0; icol < getNCols(); icol++)
    for (Id irow = 0; irow < getNRows(); irow++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;
      double value = getValue(irow, icol);
      if (FFFF(value)) continue;
      if (value > maximum) maximum = value;
    }
  if (isEqual(maximum, MINIMUM_BIG)) maximum = TEST;
  return maximum;
}

double AMatrix::getNormInf() const
{
  double norminf = 0.;
  for (Id icol = 0; icol < getNCols(); icol++)
    for (Id irow = 0; irow < getNRows(); irow++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;
      double value = getValue(irow, icol);
      if (FFFF(value)) continue;
      value = ABS(value);
      if (value > norminf) norminf = value;
    }
  return norminf;
}

void AMatrix::copyReduce(const AMatrix* x,
                         const VectorInt& validRows,
                         const VectorInt& validCols)
{
  for (Id irow = 0; irow < static_cast<Id>(validRows.size()); irow++)
    for (Id icol = 0; icol < static_cast<Id>(validCols.size()); icol++)
      setValue(irow, icol, x->getValue(validRows[irow], validCols[icol]));
}

/**
 * Copy the contents of matrix 'm' into 'this'
 * Warning: matrices must have the same dimensions (not checked)
 * @param m Input matrix
 * @param factor Multiplicative factor (applied to each element)
 */
void AMatrix::copyElements(const AMatrix& m, double factor)
{
  for (Id icol = 0; icol < m.getNCols(); icol++)
    for (Id irow = 0; irow < m.getNRows(); irow++)
      setValue(irow, icol, factor * m.getValue(irow, icol));
}

/**
 * Converts a VectorVectorDouble into a Matrix (generic)
 * Note: the input argument is stored by row (if coming from [] specification)
 * @param X Input VectorVectorDouble argument
 */
void AMatrix::_fillFromVVD(const VectorVectorDouble& X)
{
  Id nrow = static_cast<Id>(X.size());
  Id ncol = static_cast<Id>(X[0].size());

  for (Id irow = 0; irow < nrow; irow++)
    for (Id icol = 0; icol < ncol; icol++)
      setValue(irow, icol, X[irow][icol]);
}

/****************************************************************************/
/*!
 **  Check if all the elements of a matrix are non-negative
 **
 ** \return  True if the matrix is non-negative; False otherwise
 **
 ** \param[in]  verbose  True for the verbose option
 **
 *****************************************************************************/
bool AMatrix::isNonNegative(bool verbose) const
{
  for (Id irow = 0; irow < _nRows; irow++)
    for (Id icol = 0; icol < _nCols; icol++)
    {
      double value = getValue(irow, icol);
      if (value < 0.)
      {
        if (verbose) messerr("The matrix term (%d,%d) is not non-negative (%lf)", irow, icol, value);
        return false;
      }
    }
  return true;
}

/**
 * Modify the contents of the matrix so that each column has a positive sum of elements.
 * If this is not the case, simply invert the sign of the column
 */
void AMatrix::makePositiveColumn()
{
  for (Id icol = 0, ncol = getNCols(); icol < ncol; icol++)
  {
    // Extract each column
    VectorDouble column = getColumn(icol);

    // Calculate the sum of the elements
    double sum = VH::cumul(column);
    if (sum >= 0) continue;

    // Invert the sign of all its elements
    VH::multiplyConstant(column, -1.);

    // Replace the column
    setColumn(icol, column);
  }
}

void AMatrix::prodMat(const AMatrix* matY, bool transposeY)
{
  prodMatMatInPlace(this, matY, false, transposeY);
}

/**
 * @brief Perfom the algebraic equation
 * this = val1 * mat1 + val2 * mat2 + val3 * mat3
 *
 * @param val1 Coefficient of first matrx
 * @param mat1 First matrix (optional)
 * @param val2 Coefficient of second matrix
 * @param mat2 Second matrix (optional)
 * @param val3 Coefficient of third matrix
 * @param mat3 Third matrix (optional)
 */
void AMatrix::linearCombination(double val1,
                                const AMatrix* mat1,
                                double val2,
                                const AMatrix* mat2,
                                double val3,
                                const AMatrix* mat3)
{
  // Check dimensions
  if (mat1 != nullptr && !isSameSize(*mat1)) return;
  if (mat2 != nullptr && !isSameSize(*mat2)) return;
  if (mat3 != nullptr && !isSameSize(*mat3)) return;

  // Calculations
  for (Id irow = 0; irow < getNRows(); irow++)
    for (Id icol = 0; icol < getNCols(); icol++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;
      double value = 0;
      if (mat1 != nullptr) value += val1 * mat1->getValue(irow, icol);
      if (mat2 != nullptr) value += val2 * mat2->getValue(irow, icol);
      if (mat3 != nullptr) value += val3 * mat3->getValue(irow, icol);
      setValue(irow, icol, value);
    }
}

void setMultiThread(I32 nthreads)
{
  if (nthreads > 0) _globalMultiThread = nthreads;
}

I32 getMultiThread()
{
  return _globalMultiThread;
}

bool isMultiThread()
{
  return _globalMultiThread > 0;
}

void setFlagMatrixCheck(bool flagMatrixCheck)
{
  _flagMatrixCheck = flagMatrixCheck;
}
bool getFlagMatrixCheck()
{
  return _flagMatrixCheck;
}

void AMatrix::dumpRange(const char* title)
{
  VectorDouble elements = getValues();

  /* Calculate the extreme values */

  StatResults stats = ut_statistics(static_cast<Id>(elements.size()), elements.data());

  /* Printout */

  if (title != nullptr)
    message("%s\n", title);
  else
    message("Sparse matrix\n");
  message(" Descr: m=%d n=%d nnz=%d\n", getNRows(), getNCols(), _getMatrixPhysicalSize());
  if (elements.size() > 0)
    message(" Range: [%lf ; %lf] (%d/%d)\n", stats.mini, stats.maxi, stats.nvalid, elements.size());
  else
    message(" All terms are set to zero\n");
}
} // namespace gstlrn
