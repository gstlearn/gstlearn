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
#include "Matrix/MatrixFactory.hpp"
#include "Matrix/LinkMatrixSparse.hpp"
#include "Matrix/NF_Triplet.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "geoslib_define.h"

#include <iostream>

static int  globalMultiThread = 0;

AMatrix::AMatrix(int nrow, int ncol)
  : AStringable(),
    _nRows(nrow),
    _nCols(ncol),
    _flagCheckAddress(false),
    _nullTerm(0.)
{
}

AMatrix::AMatrix(const AMatrix &m)
  : AStringable(m),
    _nRows(m._nRows),
    _nCols(m._nCols),
    _flagCheckAddress(m._flagCheckAddress),
    _nullTerm(m._nullTerm)
{
}

AMatrix& AMatrix::operator=(const AMatrix &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _nRows = m._nRows;
    _nCols = m._nCols;
    _flagCheckAddress = m._flagCheckAddress;
    _nullTerm = m._nullTerm;
  }
  return *this;
}

AMatrix::~AMatrix()
{
}

void AMatrix::reset(int nrows, int ncols)
{
  // Check if numbers are valid
  if (! _isNumbersValid(nrows, ncols)) return;

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
void AMatrix::resetFromValue(int nrows, int ncols, double value)
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
void AMatrix::resetFromArray(int nrows, int ncols, const double* tab, bool byCol)
{
  reset(nrows, ncols);
  if (byCol)
  {
    int lec = 0;
    for (int icol=0; icol<ncols; icol++)
      for (int irow=0; irow<nrows; irow++)
        setValue(irow,icol,tab[lec++]);
  }
  else
  {
    int lec = 0;
    for (int irow=0; irow<nrows; irow++)
      for (int icol=0; icol<ncols; icol++)
        setValue(irow,icol,tab[lec++]);
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
void AMatrix::resetFromVD(int nrows, int ncols, const VectorDouble& tab, bool byCol)
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
    int nrows = (int) tab.size();
    int ncols = (int) tab[0].size();
    reset(nrows, ncols);
    for (int icol = 0; icol < ncols; icol++)
      for (int irow = 0; irow < nrows; irow++)
        setValue(irow, icol, tab[irow][icol]);
  }
  else
  {
    int ncols = (int) tab.size();
    int nrows = (int) tab[0].size();
    reset(nrows, ncols);
    for (int icol = 0; icol < ncols; icol++)
      for (int irow = 0; irow < nrows; irow++)
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
bool AMatrix::isValid(int irow, int icol, bool printWhyNot) const
{
  if (! _flagCheckAddress) return true;
  if (irow < 0 || irow >= getNRows())
  {
    if (printWhyNot)
      messerr("Argument 'irow' invalid: it should lie in [0;%d[",
              irow,getNRows());
    return false;
  }
  if (icol < 0 || icol >= getNCols())
  {
    if (printWhyNot)
      messerr("Argument 'icol' invalid: it should lie in [0;%d[",
              icol,getNCols());
    return false;
  }
  return true;
}

/**
 * Check that Matrix 'm' share the same dimensions as current one
 *
 * @param m Matrix to be compared to the current Matrix
 * @return true if 'm' has same dimensions as the current Matrix
 */
bool AMatrix::isSameSize(const AMatrix& m) const
{
  return (_nRows == m.getNRows() && _nCols == m.getNCols());
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
  if (! isSameSize(m)) return false;

  int ncols = getNCols();
  int nrows = getNRows();
  for (int icol=0; icol<ncols; icol++)
    for (int irow=0; irow<nrows; irow++)
    {
      double v1 = getValue(irow, icol, false);
      double v2 = m.getValue(irow, icol, false);
      if (ABS(v1 - v2) > eps)
      {
        if (printWhyNot)
        {
          messerr("Element (%d;%d) are different between:\n",irow,icol);
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

  for (int irow = 0; irow <_nRows; irow++)
    for (int icol = 0; icol < _nCols; icol++)
    {
      if (ABS(getValue(irow,icol,false) - getValue(icol,irow,false)) > eps)
      {
        if (printWhyNot)
          messerr("Elements (%d;%d)=%lf and (%d;%d)=%kf should be equal",
                  irow,icol,getValue(irow,icol,false),
                  icol,irow,getValue(icol,irow,false));
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
  for (int irow = 0; irow < getNRows(); irow++)
    for (int icol = 0; icol < getNCols(); icol++)
    {
      double refval = (irow == icol) ? 1. : 0.;
      if (ABS(getValue(irow, icol, false) - refval) > EPSILON10)
      {
        if (printWhyNot)
          messerr("The term (%d,%d) should be equal to %lf (%lf)", irow + 1,
                  icol + 1, refval, getValue(irow, icol,false));
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
 * Fill 'this' with the constant 'value'
 * @param value Constant value used for filling 'this'
 */
void AMatrix::fill(double value)
{
  for (int rank = 0, n = _getMatrixPhysicalSize(); rank < n; rank++)
    _setValueByRank(rank, value);
}

int AMatrix::_getMatrixPhysicalSize() const
{
  return (getNRows() * getNCols());
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
void AMatrix::_setValues(const double *values, bool byCol)
{
  int ecr = 0;
  if (byCol)
  {
    for (int icol = 0; icol < getNCols(); icol++)
      for (int irow = 0; irow < getNRows(); irow++, ecr++)
        setValue(irow, icol, values[ecr]);
  }
  else
  {
    for (int irow = 0; irow < getNRows(); irow++)
      for (int icol = 0; icol < getNCols(); icol++, ecr++)
        setValue(irow, icol, values[ecr]);
  }
}
#endif

void AMatrix::fillRandom(int seed, double zeroPercent)
{
  law_set_random_seed(seed);

  double value = 0.;
  for (int irow = 0; irow < _nRows; irow++)
    for (int icol = 0; icol <_nCols; icol++)
    {
      if (! _isPhysicallyPresent(irow, icol)) continue;
      if (law_uniform(0.,1.) < zeroPercent)
        value = 0.;
      else
        value = law_gaussian();
      setValue(irow,icol,value);
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
  if ((int) values.size() != size())
  {
    messerr("Inconsistency between 'values' and Matrix Dimension");
    messerr("Operation cancelled");
    return;
  }
  _setValues(values.data(),byCol);
}

void AMatrix::setIdentity(double value)
{
  for (int icol = 0; icol < _nCols; icol++)
    for (int irow = 0; irow < _nRows; irow++)
      setValue(irow, icol, value * (irow == icol), false);
}

/**
 *
 * @param v Add a scalar value to all (valid) terms of the current matrix
 */
void AMatrix::addScalar(double v)
{
  if (isZero(v)) return;
  for (int rank = 0; rank < _getMatrixPhysicalSize(); rank++)
    _setValueByRank(rank, _getValueByRank(rank) + v);
}

/**
 *
 * @param v Add constant value to the diagonal of the current Matrix
 */
void AMatrix::addScalarDiag(double v)
{
  if (isZero(v)) return;
  for (int irow = 0; irow < _nRows; irow++)
  {
    for (int icol = 0; icol < _nCols; icol++)
    {
      if (irow == icol)
      {
        int rank = _getIndexToRank(irow, icol);
        _setValueByRank(rank, _getValueByRank(rank) + v);
      }
    }
  }
}

/**
 *
 * @param v Multiply all the terms of the matrix by the scalar 'v'
 */
void AMatrix::prodScalar(double v)
{
  if (isOne(v)) return;
  for (int rank = 0; rank < _getMatrixPhysicalSize(); rank++)
    _setValueByRank(rank, _getValueByRank(rank) * v);
}

/**
 * Returns 'y' = 'this' %*% 'x'
 * @param x Input vector
 * @param y Output vector obtained by multiplying 'inv' by current Matrix
 * @param transpose True if the matrix 'this' must be transposed
 */
void AMatrix::prodMatVecInPlace(const VectorDouble& x, VectorDouble& y, bool transpose) const
{
  y.fill(0.,y.size());
  _addProdMatVecInPlaceToDestPtr(x.data(),y.data(),transpose);
}

int AMatrix::addProdMatVecInPlace(const constvect x,
                                  vect y,
                                  bool transpose) const
{
  if (_flagCheckAddress)
  {
    bool error = false;
    if (!transpose)
    {
      error = ((int) x.size() != _nCols || (int) y.size() != _nRows);
    }
    else
    {
      error = ((int) x.size() != _nRows || (int) y.size() != _nCols);
    }
    if (error)
    {
      messerr("Inconsistency between:");
      messerr("- the dimension of 'x' = %d", (int) x.size());
      messerr("- the dimension of 'y' = %d", (int) y.size());
      messerr("- the matrix: number of rows (%d) and columns (%d)", _nRows, _nCols);
      return 1;
    }
  }
  _addProdMatVecInPlaceToDestPtr(x.data(), y.data(), transpose);
  return 0;
}

int AMatrix::prodMatVecInPlace(const constvect x, vect y, bool transpose) const
{
  if (_flagCheckAddress)
  {
    bool error = false;
    if (!transpose)
    {
      error = ((int) x.size() != _nCols || (int) y.size() != _nRows);
    }
    else
    {
      error = ((int) x.size() != _nRows || (int) y.size() != _nCols);
    }
    if (error)
    {
      messerr("Inconsistency between:");
      messerr("- the dimension of 'x' = %d", (int) x.size());
      messerr("- the dimension of 'y' = %d", (int) y.size());
      messerr("- the matrix: number of rows (%d) and columns (%d)", _nRows, _nCols);
      return 1;
    }
  }
  _prodMatVecInPlacePtr(x.data(), y.data(), transpose);
  return 0;
}

void AMatrix::prodMatVecInPlacePtr(const double* x, double* y, bool transpose) const
{
  _prodMatVecInPlacePtr(x, y, transpose);
}

void AMatrix::prodVecMatInPlace(const VectorDouble& x, VectorDouble& y, bool transpose) const
{
  if (_flagCheckAddress)
  {
    bool error = false;
    if (!transpose)
      error = ((int) x.size() != _nRows || (int) y.size() != _nCols);
    else
      error = ((int) x.size() != _nCols || (int) y.size() != _nRows);
    if (error)
    {
      messerr("Inconsistency between:");
      messerr("- the dimension of 'x' = %d", (int) x.size());
      messerr("- the dimension of 'y' = %d", (int) y.size());
      messerr("- the matrix: number of rows (%d) and columns (%d)", _nRows, _nCols);
      return;
    }
  }
  _prodVecMatInPlacePtr(x.data(), y.data(), transpose);
}

void AMatrix::prodVecMatInPlacePtr(const double* x, double* y, bool transpose) const
{
  _prodVecMatInPlacePtr(x, y, transpose);
}

bool AMatrix::needToReset(int nrows, int ncols)
{
  return nrows != getNRows() || ncols != getNCols() || _needToReset(nrows, ncols);
}

bool AMatrix::_needToReset(int nrows, int ncols)
{
  DECLARE_UNUSED(nrows,ncols)
  return false;
}
/**
 * @brief Resize the matrix to new dimensions
 *        (this method doesn't change the storage type)
 * 
 * @param nrows New number of rows
 * @param ncols New number of columns
 */
void AMatrix::resize(int nrows, int ncols)
{
  // Check if nothing is to be done
  if (!needToReset(nrows, ncols)) 
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
void AMatrix::addMatInPlace(const AMatrix& y, double cx, double cy)
{
  if (! isSameSize(y))
  {
    messerr("Matrices 'y' and 'this' should have the same size");
    return;
  }

  for (int irow = 0; irow < _nRows; irow++)
    for (int icol = 0; icol < _nCols; icol++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;
      setValue(irow, icol, cx * getValue(irow, icol) + cy * y.getValue(irow, icol), false);
    }
}

/**
 * Store the product of 'x'(or 't(x)') by 'y' (or 't(y') in this
 * @param x First Matrix
 * @param y Second matrix
 * @param transposeX True if first matrix must be transposed
 * @param transposeY True if second matrix must be transposed
 */
void AMatrix::prodMatMatInPlace(const AMatrix *x,
                                const AMatrix *y,
                                bool transposeX,
                                bool transposeY)
{
  int ni1 = (transposeX) ? x->getNCols() : x->getNRows();
  int nm1 = (transposeX) ? x->getNRows() : x->getNCols();
  int ni2 = (transposeY) ? y->getNRows() : y->getNCols();
  int nm2 = (transposeY) ? y->getNCols() : y->getNRows();

  if (nm1 != ni2)
  {
    messerr("Matrices 'x' and 'y' should have matching dimensions");
    return;
  }

  if (!_checkLink(x->getNRows(), x->getNCols(), transposeX,
                  y->getNRows(), y->getNCols(), transposeY)) return;

  for (int irow = 0; irow < ni1; irow++)
  {
    for (int icol = 0; icol < nm2; icol++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;

      double value = 0.;
      for (int k = 0; k < nm1; k++)
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
  int n1 = (transpose) ? a->getNCols() : a->getNRows();
  int n2 = (transpose) ? a->getNRows() : a->getNCols();

  if (!_checkLink(a->getNRows(), a->getNCols(), transpose,
                  m->getNRows(), m->getNCols(), false,
                  a->getNRows(), a->getNCols(), !transpose)) return;

  for (int i = 0; i < n1; i++)
    for (int j = 0; j < n1; j++)
    {
      if (!_isPhysicallyPresent(i, j)) continue;
      double value = 0;
      for (int k = 0; k < n2; k++)
        for (int l = 0; l < n2; l++)
        {
          double a_ik = (transpose) ? a->getValue(k,i) : a->getValue(i,k);
          double a_lj = (transpose) ? a->getValue(l,j) : a->getValue(j,l);
          value += a_ik * m->getValue(k,l) * a_lj;
        }
      setValue(i, j, value);
    }
}

void AMatrix::prodNormMatVecInPlace(const AMatrix &a, const VectorDouble& vec, bool transpose)
{
  int n1 = (transpose) ? a.getNCols() : a.getNRows();
  int n2 = (transpose) ? a.getNRows() : a.getNCols();

  if (!_checkLink(getNRows(), getNCols(), false, a.getNRows(), a.getNCols(), transpose,
                 (int) vec.size(), 1, false)) return;

  for (int i = 0; i < n1; i++)
    for (int j = 0; j < n1; j++)
    {
      if (!_isPhysicallyPresent(i, j)) continue;
      double value = 0;
      for (int k = 0; k < n2; k++)
      {
        double a_ik = (transpose) ? a.getValue(k,i) : a.getValue(i,k);
        double a_kj = (transpose) ? a.getValue(j,k) : a.getValue(k,j);
        double vecvalue = (vec.empty()) ? 1. : vec[k];
        value += a_ik * vecvalue * a_kj;
      }
      setValue(i, j, value);
    }
}

void AMatrix::multiplyRow(const VectorDouble& vec)
{
  if (_nRows != (int) vec.size())
  {
    messerr("The size of 'vec' must match the number of rows. Nothing is done");
    return;
  }
  for (int irow = 0; irow < _nRows; irow++)
    for (int icol = 0; icol < _nCols; icol++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;
      setValue(irow, icol, getValue(irow, icol, false) * vec[irow], false);
    }
}

void AMatrix::divideRow(const VectorDouble& vec)
{
  if (_nRows != (int) vec.size())
  {
    messerr("The size of 'vec' must match the number of rows. Nothing is done");
    return;
  }
  for (int irow = 0; irow < _nRows; irow++)
    for (int icol = 0; icol < _nCols; icol++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;
      setValue(irow, icol, getValue(irow, icol, false) / vec[irow], false);
    }
}

void AMatrix::multiplyColumn(const VectorDouble& vec)
{
  if (_nCols != (int) vec.size())
  {
    messerr("The size of 'vec' must match the number of columns. Nothing is done");
    return;
  }
  for (int irow = 0; irow < _nRows; irow++)
    for (int icol = 0; icol < _nCols; icol++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;
      setValue(irow, icol, getValue(irow, icol, false) * vec[icol], false);
    }
}
void AMatrix::divideColumn(const VectorDouble& vec)
{
  if (_nCols != (int) vec.size())
  {
    messerr("The size of 'vec' must match the number of columns. Nothing is done");
    return;
  }
  for (int irow = 0; irow < _nRows; irow++)
    for (int icol = 0; icol < _nCols; icol++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;
      setValue(irow, icol, getValue(irow, icol, false) / vec[icol], false);
    }
}

VectorDouble AMatrix::prodVecMat(const VectorDouble& x, bool transpose) const
{
  VectorDouble y;

  if (!_checkLink((int) x.size(), 1, false, getNRows(), getNCols(), transpose)) return y;

  int nval = (! transpose) ? _nCols : _nRows;
  y.resize(nval, 0.);
  prodVecMatInPlace(x, y, transpose);
  return y;
}

VectorDouble AMatrix::prodMatVec(const VectorDouble& x, bool transpose) const
{
  VectorDouble y;

  if (!_checkLink(getNRows(), getNCols(), transpose, (int) x.size(), 1, false)) return y;

  int nval = (transpose) ? _nCols : _nRows;
  y.resize(nval, 0.);
  prodMatVecInPlace(x, y, transpose);
  return y;
}

double AMatrix::quadraticMatrix(const VectorDouble& x, const VectorDouble& y)
{
  if (!_checkLink((int) x.size(), 1, true,
                  getNRows(), getNCols(), false,
                  (int) y.size(), 1, false)) return TEST;

  VectorDouble left(_nRows);
  prodMatVecInPlace(y, left, false);
  return VH::innerProduct(x, left);
}

int AMatrix::invert()
{
  if (! isSquare())
  {
    messerr("'invert' method is restricted to Square Matrices");
    return 1;
  }
  return _invert();
}

int AMatrix::solve(const VectorDouble& b, VectorDouble& x) const
{
  if (! isSquare())
  {
    messerr("'solve' method is limited to Square Matrices");
    return 1;
  }
  if ((int) b.size() != _nRows || (int) x.size() != _nRows)
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

  sstr << "- Number of rows    = " <<  _nRows << std::endl;
  sstr << "- Number of columns = " <<  _nCols << std::endl;

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

void AMatrix::_clear()
{
  _setNRows(0);
  _setNCols(0);
}

bool AMatrix::_isNumbersValid(int nrows, int ncols) const
{
  if (! _flagCheckAddress) return true;
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

bool AMatrix::_isRowValid(int irow) const
{
  if (!_flagCheckAddress) return true;
  return checkArg("Row index invalid", irow, getNRows());
}

bool AMatrix::_isColumnValid(int icol) const
{
  if (!_flagCheckAddress) return true;
  return checkArg("Column index invalid", icol, getNCols());
}

bool AMatrix::_isIndexValid(int irow, int icol) const
{
  if (! _flagCheckAddress) return true;
  if (! _isRowValid(irow)) return false;
  if (! _isColumnValid(icol)) return false;
  return true;
}

bool AMatrix::_isRowVectorConsistent(const VectorDouble& tab) const
{
  if ((int) tab.size() != getNRows())
  {
    messerr("Argument vector size should match the number of rows");
    return false;
  }
  return true;
}

bool AMatrix::_isColVectorConsistent(const VectorDouble& tab) const
{
  if ((int) tab.size() != getNCols())
  {
    messerr("Argument vector size should match the number of columns");
    return false;
  }
  return true;
}

bool AMatrix::_isVectorSizeConsistent(const VectorDouble &tab) const
{
  int nrows = getNRows();
  int ncols = getNCols();
  if ((int) tab.size() != nrows * ncols)
  {
    messerr("The argument 'tab'(%d) does not have correct dimension (%d)",
            (int) tab.size(), nrows * ncols);
    return false;
  }
  return true;
}

bool AMatrix::_isColumnSizeConsistent(const VectorDouble &tab) const
{
  int nrows = getNRows();
  if ((int) tab.size() != nrows)
  {
    messerr("The argument 'tab'(%d) does not have correct dimension (%d)",
            (int) tab.size(), nrows);
    return false;
  }
  return true;
}

bool AMatrix::_isRowSizeConsistent(const VectorDouble &tab) const
{
  int ncols = getNCols();
  if ((int) tab.size() != ncols)
  {
    messerr("The argument 'tab'(%d) does not have correct dimension (%d)",
            (int) tab.size(), ncols);
    return false;
  }
  return true;
}

bool AMatrix::_isRankValid(int rank) const
{
  if (! _flagCheckAddress) return true;
  return (rank >= 0 && rank < _getMatrixPhysicalSize());
}

void AMatrix::dumpElements(const String& title, int ifrom, int ito) const
{
  mestitle(1, "%s", title.c_str());
  for (int rank = ifrom; rank < ito; rank++)
  {
    if (_isRankValid(rank))
      message("Element %d = %lf\n", rank, _getValueByRank(rank));
  }
}

void AMatrix::dumpStatistics(const String& title) const
{
  message("%s : %d rows and %d columns\n", title.c_str(), _nRows, _nCols);
}
/**
 * Check that a set of matrices (or vectors) has the correct linkage
 * @param nrow1       Number of rows in the first matrix
 * @param ncol1       Number of columns in the first matrix
 * @param transpose1  True if the first matrix must be transposed
 * @param nrow2       Number of rows in the second matrix
 * @param ncol2       Number of columns in the second matrix
 * @param transpose2  True if the second matrix must be transposed
 * @param nrow3       Number of rows in the third matrix
 * @param ncol3       Number of Columns in the third matrix
 * @param transpose3  True if the third matrix must be transposed
 * @return True if the linkage is correct
 *
 * @remark A vector must be defined with its 'ncol' set to 1
 * @remark If a matrix is missing, set the corresponding 'nrow' to 0
 * @remark This test is not performed if _getFlagCheckAddress() returns false
 */
bool AMatrix::_checkLink(int nrow1,
                         int ncol1,
                         bool transpose1,
                         int nrow2,
                         int ncol2,
                         bool transpose2,
                         int nrow3,
                         int ncol3,
                         bool transpose3) const
{
  if (! _getFlagCheckAddress()) return true;
  int level = 0;
  int ncur = getNRows();

  if (nrow1 > 0)
  {
    int nr1 = (!transpose1) ? nrow1 : ncol1;
    int nc1 = (!transpose1) ? ncol1 : nrow1;
    if (ncur != nr1) level = 1;
    ncur = nc1;
  }

  if (nrow2 > 0)
  {
    int nr2 = (!transpose2) ? nrow2 : ncol2;
    int nc2 = (!transpose2) ? ncol2 : nrow2;
    if (ncur != nr2) level = 2;
    ncur = nc2;
  }

  if (nrow3 > 0)
  {
    int nr3 = (!transpose3) ? nrow3 : ncol3;
    int nc3 = (!transpose3) ? ncol3 : nrow3;
    if (ncur != nr3) level = 3;
    ncur = nc3;
  }

  if (ncur != getNCols()) level = -1;

  if (level != 0)
  {
    messerr("Error in the Linkage of matrices: Level = %d", level);
    messerr("Operation is cancelled");
  }
  return level == 0;
}

/**
 * From a matrix of any type, creates the three vectors of the triplet
 * (specific format for creating efficiently a Sparse matrix)
 * It only takes the only non-zero elements of the matrix
 */
NF_Triplet AMatrix::getMatrixToTriplet(int shiftRow, int shiftCol) const
{
  NF_Triplet NF_T;

  for (int icol = 0; icol < _nCols; icol++)
    for (int irow = 0; irow < _nRows; irow++)
    {
      if (!isValid(irow, icol)) continue;
      double value = getValue(irow, icol);
      if (isZero(value)) continue;
      NF_T.add(irow+shiftRow, icol+shiftCol, value);
    }

  return NF_T;
}

VectorDouble AMatrix::getValues(bool byCol) const
{
  VectorDouble vect(_nCols * _nRows);
  VectorDouble::iterator itvect(vect.begin());

  if (byCol)
  {
    for (int icol = 0; icol < _nCols; icol++)
      for (int irow = 0; irow < _nRows; irow++)
      {
        (*itvect) = getValue(irow, icol, false);
        itvect++;
      }
  }
  else
  {
    for (int irow = 0; irow < _nRows; irow++)
      for (int icol = 0; icol < _nCols; icol++)
      {
        (*itvect) = getValue(irow, icol, false);
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
  for (int icol = 0; icol < _nCols; icol++)
    for (int irow = 0; irow < _nRows; irow++)
    {
      double v1 = (    isValid(irow, icol)) ?     getValue(irow, icol) : 0.;
      double v2 = (mat.isValid(irow, icol)) ? mat.getValue(irow, icol) : 0.;
      diff += ABS(v1 - v2);
    }
  return diff;
}

VectorDouble AMatrix::getDiagonal(int shift) const
{
  if (! isSquare())
  {
    messerr("This function is only valid for Square matrices");
    return VectorDouble();
  }

  VectorDouble vect;
  for (int rank = 0; rank < getNRows(); rank++)
  {
    int irow = rank;
    if (shift < 0) irow += shift;
    if (irow < 0 || irow >= getNRows()) continue;
    int icol = rank;
    if (shift > 0) icol += shift;
    if (icol < 0 || icol >= getNCols()) continue;
    vect.push_back(getValue(irow,icol));
  }
  return vect;
}

/**
 * Reset the contents of a matrix by setting all terms to 0 and
 * update diagonal terms from the input argument 'tab'
 * @param tab Input vector to be copied to the diagonal of the output matrix
 * @param flagCheck When True, check the input arguments
 */
void AMatrix::setDiagonal(const VectorDouble& tab, bool flagCheck)
{
  int nrows = getNRows();
  if (! isSquare())
  {
    messerr("This function is only valid for Square matrices. Nothing is done");
    return;
  }
  if (flagCheck)
  {
    if (! _isRowSizeConsistent(tab)) return;
  }

  for (int irow = 0; irow < nrows; irow++)
  {
    if (irow >= getNCols()) continue;
    setValue(irow,irow,tab[irow]);
  }
}

void AMatrix::setDiagonalToConstant(double value)
{
  if (! isSquare())
  {
    messerr("This function is only valid for Square matrices. Nothing is done");
    return;
  }

  for (int irow = 0; irow < getNRows(); irow++)
  {
    int icol = irow;
    if (icol < 0 || icol >= getNCols()) continue;
    setValue(irow,icol,value);
  }
}

/*! Extract a Row */
VectorDouble AMatrix::getRow(int irow) const
{
  VectorDouble vect;
  if (!checkArg("Incorrect argument 'irow'", irow, getNRows())) return vect;

  for (int icol = 0; icol < getNCols(); icol++)
    vect.push_back(getValue(irow,icol));
  return vect;
}

/*! Set the contents of a Row */
void AMatrix::setRow(int irow, const VectorDouble& tab, bool flagCheck)
{
  if (irow < 0 || irow >= getNRows())
    my_throw("Incorrect argument 'irow'");
  if ((int) tab.size() != getNCols())
    my_throw("Incorrect dimension of 'tab'");

  for (int icol = 0; icol < getNCols(); icol++)
    setValue(irow,icol,tab[icol],flagCheck);
}

/*! Extract a Column */
VectorDouble AMatrix::getColumn(int icol) const
{
  if (icol < 0 || icol >= getNCols())
    my_throw("Incorrect argument 'icol'");

  VectorDouble vect;
  for (int irow = 0; irow < getNRows(); irow++)
    vect.push_back(getValue(irow,icol));
  return vect;
}

VectorDouble AMatrix::getColumnByRowRange(int icol, int rowFrom, int rowTo) const
{
  if (icol < 0 || icol >= getNCols())
    my_throw("Incorrect argument 'icol'");

  VectorDouble vect;
  for (int irow = rowFrom; irow < rowTo; irow++)
    vect.push_back(getValue(irow, icol));
  return vect;
}
/*! Set the contents of a Column */
void AMatrix::setColumn(int icol, const VectorDouble& tab, bool flagCheck)
{
  if (icol < 0 || icol >= getNCols())
    my_throw("Incorrect argument 'icol'");
  if ((int) tab.size() != getNRows())
    my_throw("Incorrect dimension of 'tab'");

  for (int irow = 0; irow < getNRows(); irow++)
    setValue(irow,icol,tab[irow], flagCheck);
}

/*! Checks if a Column is valid (contains a non TEST value) */
bool AMatrix::isColumnDefined(int icol) const
{
  if (icol < 0 || icol >= getNCols())
    my_throw("Incorrect argument 'icol'");

  for (int irow = 0; irow < getNRows(); irow++)
  {
    if (! FFFF(getValue(irow,icol))) return true;
  }
  return false;
}

/*! Checks if a Row is valid (contains a non TEST value) */
bool AMatrix::isRowDefined(int irow) const
{
  if (irow < 0 || irow >= getNRows())
    my_throw("Incorrect argument 'irow'");

  for (int icol = 0; icol < getNCols(); icol++)
  {
    if (! FFFF(getValue(irow,icol))) return true;
  }
  return false;
}

/*! Define the number of defined columns */
int AMatrix::getNColDefined() const
{
  int ncol = 0;
  for (int icol = 0; icol < getNCols(); icol++)
  {
    if (isColumnDefined(icol)) ncol++;
  }
  return ncol;
}

/*! Define the number of defined rows */
int AMatrix::getNRowDefined() const
{
  int nrow = 0;
  for (int irow = 0; irow < getNRows(); irow++)
  {
    if (isRowDefined(irow)) nrow++;
  }
  return nrow;
}

void AMatrix::addValue(int irow, int icol, double value)
{
  double oldval = getValue(irow, icol);
  if (FFFF(oldval)) return;
  setValue(irow, icol, oldval + value);
}

double AMatrix::getMeanByColumn(int icol) const
{
  double cumul = 0.;
  double count = 0.;
  for (int irow = 0; irow < getNRows(); irow++)
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
  double minimum = 1.e30;
  for (int icol = 0; icol < getNCols(); icol++)
    for (int irow = 0; irow < getNRows(); irow++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;
      double value = getValue(irow, icol, false);
      if (FFFF(value)) continue;
      if (value < minimum) minimum = value;
    }
  if (isEqual(minimum,1.e30)) minimum = TEST;
  return minimum;
}

double AMatrix::getMaximum() const
{
  double maximum = -1.e30;
  for (int icol = 0; icol < getNCols(); icol++)
    for (int irow = 0; irow < getNRows(); irow++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;
      double value = getValue(irow, icol, false);
      if (FFFF(value)) continue;
      if (value > maximum) maximum = value;
    }
  if (isEqual(maximum,-1.e30)) maximum = TEST;
  return maximum;
}

double AMatrix::getNormInf() const
{
  double norminf = 0.;
  for (int icol = 0; icol < getNCols(); icol++)
    for (int irow = 0; irow < getNRows(); irow++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;
      double value = getValue(irow, icol, false);
      if (FFFF(value)) continue;
      value = ABS(value);
      if (value > norminf) norminf = value;
  }
  return norminf;
}

double& AMatrix::_getValueRef(int irow, int icol)
{
  DECLARE_UNUSED(irow);
  DECLARE_UNUSED(icol);
  return _nullTerm;
}

void AMatrix::copyReduce(const AMatrix *x,
                         const VectorInt &validRows,
                         const VectorInt &validCols)
{
  for (int irow = 0; irow < (int) validRows.size(); irow++)
    for (int icol = 0; icol < (int) validCols.size(); icol++)
      setValue(irow, icol, x->getValue(validRows[irow], validCols[icol], false), false);
}

/**
 * Copy the contents of matrix 'm' into 'this'
 * Warning: matrices must have the same dimensions (not checked)
 * @param m Input matrix
 * @param factor Multiplicative factor (applied to each element)
 */
void AMatrix::copyElements(const AMatrix &m, double factor)
{
  for (int icol = 0; icol < m.getNCols(); icol++)
    for (int irow = 0; irow < m.getNRows(); irow++)
      setValue(irow, icol, factor * m.getValue(irow, icol, false), false);
}

/**
 * Converts a VectorVectorDouble into a Matrix (generic)
 * Note: the input argument is stored by row (if coming from [] specification)
 * @param X Input VectorVectorDouble argument
 */
void AMatrix::_fillFromVVD(const VectorVectorDouble& X)
{
  int nrow = (int) X.size();
  int ncol = (int) X[0].size();

  for (int irow = 0; irow < nrow; irow++)
    for (int icol = 0; icol < ncol; icol++)
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
  for (int irow = 0; irow < _nRows; irow++)
    for (int icol = 0; icol < _nCols; icol++)
    {
      double value = getValue(irow,icol);
      if (value < 0.)
      {
        if (verbose) messerr("The matrix term (%d,%d) is not non-negative (%lf)",irow,icol,value);
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
  for (int icol = 0, ncol = getNCols(); icol < ncol; icol++)
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

void AMatrix::prodMatInPlace(const AMatrix *matY, bool transposeY)
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
  if (mat1 != nullptr && ! isSameSize(*mat1))
  {
    messerr("AMatrix::linearCombination: Dimensions of 'mat1' do not match "
            "dimensions of current matrix. Nothing is done");
    return;
  }
  if (mat2 != nullptr && ! isSameSize(*mat2))
  {
    messerr("AMatrix::linearCombination: Dimensions of 'mat2' do not match "
            "dimensions of current matrix. Nothing is done");
    return;
  }
  if (mat3 != nullptr && !isSameSize(*mat3))
  {
    messerr("AMatrix::linearCombination: Dimensions of 'mat3' do not match "
            "dimensions of current matrix. Nothing is done");
    return;
  }

  // Calculations
  for (int irow = 0; irow < getNRows(); irow++)
    for (int icol = 0; icol < getNCols(); icol++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;
      double value = 0;
      if (mat1 != nullptr) value += val1 * mat1->getValue(irow, icol);
      if (mat2 != nullptr) value += val2 * mat2->getValue(irow, icol);
      if (mat3 != nullptr) value += val3 * mat3->getValue(irow, icol);
      setValue(irow, icol, value);
    }
}

void setMultiThread(int nthreads)
{
  if (nthreads > 0) globalMultiThread = nthreads;
}

int getMultiThread()
{
  return globalMultiThread;
}

bool isMultiThread()
{
  return globalMultiThread > 0;
}

