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
#include "Basic/VectorHelper.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"

#include <iostream>
#include <iomanip>

/**
 * This function switch ON/OFF the ability to use Eigen library for Algebra
 */
static bool globalFlagEigen   = true;
static int  globalMultiThread = 0;

AMatrix::AMatrix(int nrow, int ncol, int opt_eigen)
    : AStringable(),
      _nRows(nrow),
      _nCols(ncol),
      _flagEigen(),
      _flagCheckAddress(false),
      _nullTerm(0.)
{
  (void) _isNumbersValid(nrow, ncol);
  _flagEigen = _defineFlagEigen(opt_eigen);
}

AMatrix::AMatrix(const AMatrix &m)
    : AStringable(m),
      _nRows(m._nRows),
      _nCols(m._nCols),
      _flagEigen(m._flagEigen),
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
    _flagEigen = m ._flagEigen;
    _flagCheckAddress = m._flagCheckAddress;
    _nullTerm = m._nullTerm;
  }
  return *this;
}

AMatrix::~AMatrix()
{
}

void AMatrix::init(int nrows, int ncols, int opt_eigen)
{
  _nRows = nrows;
  _nCols = ncols;
  _flagEigen = _defineFlagEigen(opt_eigen);
  _allocate();
}

bool AMatrix::isSquare(bool printWhyNot) const
{
  if (empty()) return false;
  if (_nRows != _nCols)
  {
    if (printWhyNot)
      messerr("The number of rows (%d) should math the number of columns (%d)",
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
 * @param printWhyNot Print the message is the answer if false
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

bool AMatrix::isSame(const AMatrix& m, double eps)
{
  if (! isSameSize(m)) return false;

  int ncols = getNCols();
  int nrows = getNRows();
  for (int icol=0; icol<ncols; icol++)
    for (int irow=0; irow<nrows; irow++)
    {
      if (ABS(getValue(irow, icol) - m.getValue(irow, icol)) > eps) return false;
    }
  return true;
}

void AMatrix::reset(int nrows, int ncols, double value, int opt_eigen)
{
  if (! _isNumbersValid(nrows, ncols)) return;
  _deallocate();
  _nRows = nrows;
  _nCols = ncols;
  _flagEigen = _defineFlagEigen(opt_eigen);
  _allocate();
  fill(value);
  _clearDecoration();
}

void AMatrix::resetFromArray(int nrows, int ncols, const double* tab, bool byCol, int opt_eigen)
{
  if (! _isNumbersValid(nrows, ncols)) return;
  _deallocate();
  _nRows = nrows;
  _nCols = ncols;
  _flagEigen = _defineFlagEigen(opt_eigen);
  _allocate();
  int lec = 0;
  if (byCol)
  {
    for (int icol=0; icol<ncols; icol++)
      for (int irow=0; irow<nrows; irow++)
        _setValue(irow,icol,tab[lec++]);
  }
  else
  {
    for (int irow=0; irow<nrows; irow++)
      for (int icol=0; icol<ncols; icol++)
        _setValue(irow,icol,tab[lec++]);
  }
  _clearDecoration();
}

void AMatrix::resetFromVD(int nrows, int ncols, const VectorDouble& tab, bool byCol, int opt_eigen)
{
  if (! _isNumbersValid(nrows, ncols)) return;
  resetFromArray(nrows, ncols, tab.data(), byCol, opt_eigen);
}

void AMatrix::resetFromVVD(const VectorVectorDouble& tab, bool byCol, int opt_eigen)
{
  if (byCol)
  {
    _nRows = (int) tab.size();
    _nCols = (int) tab[0].size();
    _flagEigen = _defineFlagEigen(opt_eigen);
    _allocate();
    for (int icol = 0; icol < _nCols; icol++)
      for (int irow = 0; irow < _nRows; irow++)
        _setValue(irow, icol, tab[irow][icol]);
  }
  else
  {
    _nCols = (int) tab.size();
    _nRows = (int) tab[0].size();
    _flagEigen = _defineFlagEigen(opt_eigen);
    _allocate();
    for (int icol = 0; icol < _nCols; icol++)
      for (int irow = 0; irow < _nRows; irow++)
        _setValue(irow, icol, tab[icol][irow]);
  }
  _clearDecoration();
}

void AMatrix::fillRandom(int seed, double zeroPercent)
{
  law_set_random_seed(seed);

  double value = 0.;
  for (int irow = 0; irow < _nRows; irow++)
    for (int icol = 0; icol <_nCols; icol++)
    {
      if (! _isPhysicallyPresent(irow, icol)) continue;
      if (!mustBeDiagCst() && law_uniform(0.,1.) < zeroPercent)
        value = 0.;
      else
        value = law_gaussian();
      setValue(irow,icol,value);
    }
}

bool AMatrix::isSymmetric(bool printWhyNot) const
{
  if (empty() || ! isSquare()) return false;

  for (int irow = 0; irow <_nRows; irow++)
    for (int icol = 0; icol < _nCols; icol++)
    {
      if (getValue(irow,icol) != getValue(icol,irow))
      {
        if (printWhyNot)
          messerr("Elements (%d;%d)=%lf and (%d;%d)=%kf should be equal",
                  irow,icol,getValue(irow,icol),
                  icol,irow,getValue(icol,irow));
        return false;
      }
    }
  return true;
}

bool AMatrix::isIdentity(bool printWhyNot) const
{
  for (int irow = 0; irow < getNRows(); irow++)
    for (int icol = 0; icol < getNCols(); icol++)
    {
      double refval = (irow == icol) ? 1. : 0.;
      if (ABS(getValue_(irow, icol) - refval) > EPSILON10)
      {
        if (printWhyNot)
          messerr("The term (%d,%d) should be equal to %lf (%lf)", irow + 1,
                  icol + 1, refval, getValue(irow, icol));
        return false;
      }
    }
  return true;
}

bool AMatrix::isDiagonal(bool printWhyNot) const
{
  if (empty() || ! isSquare()) return false;

  for (int irow = 0; irow <_nRows; irow++)
    for (int icol = 0; icol < _nCols; icol++)
    {
      if (irow != icol)
      {
        if (ABS(getValue_(irow,icol)) > EPSILON10)
        {
          if (printWhyNot)
            messerr("The element (%d;%d)=%lf should be zero", irow, icol,
                    getValue_(irow, icol));
          return false;
        }
      }
    }
  return true;
}

bool AMatrix::isDiagCst(bool printWhyNot) const
{
  if (empty() || ! isSquare()) return false;

  double refval = TEST;
  for (int irow = 0; irow <_nRows; irow++)
    for (int icol = 0; icol < _nCols; icol++)
    {
      if (irow != icol)
      {
        if (ABS(getValue(irow, icol)) > EPSILON10)
        {
          if (printWhyNot)
            messerr("The element (%d,%d) is not zero (%lf)", irow, icol,
                    getValue(irow, icol));
          return false;
        }
      }
      else
      {
        if (FFFF(refval))
          refval = getValue(irow,icol);
        else
          if (ABS(refval - getValue(irow,icol)) > EPSILON10)
          {
            if (printWhyNot)
              messerr("The element(%d,%d) has value (%lf) different from a previous one (%lf)",
                      irow,icol,getValue(irow,icol),refval);
            return false;
          }
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

/*! Gets the value at row 'irow' and column 'icol' (no test) */
double AMatrix::getValue_(int irow, int icol) const
{
  return _getValue(irow, icol);
}

/*! Gets the value at row 'irow' and column 'icol' */
double AMatrix::getValue(int irow, int icol) const
{
  if (! _isIndexValid(irow, icol)) return TEST;
  return _getValue(irow, icol);
}

/*! Sets the value at row 'irow' and column 'icol' (no test) */
void AMatrix::setValue_(int irow, int icol, double value)
{
  return _setValue(irow, icol, value);
}

/*! Sets the value at row 'irow' and column 'icol' */
void AMatrix::setValue(int irow, int icol, double value)
{
  if (! _isIndexValid(irow, icol)) return;
  return _setValue(irow, icol, value);
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
#endif

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

void AMatrix::setValuesFromTriplet(const NF_Triplet& NF_T)
{
  for (int i = 0; i < NF_T.number; i++)
    setValue(NF_T.rows[i], NF_T.cols[i], NF_T.values[i]);
}

void AMatrix::setIdentity(double value)
{
  for (int icol = 0; icol < _nCols; icol++)
    for (int irow = 0; irow < _nRows; irow++)
      setValue_(irow, icol, value * (irow == icol));
}

/**
 *
 * @param v Add a scalar value to all (valid) terms of the current matrix
 */
void AMatrix::addScalar(double v)
{
  if (v == 0.) return;
  for (int rank = 0; rank < _getMatrixPhysicalSize(); rank++)
    _setValueByRank(rank, _getValueByRank(rank) + v);
}

/**
 *
 * @param v Add constant value to the diagonal of the current Matrix
 */
void AMatrix::addScalarDiag(double v)
{
  if (v == 0.) return;
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
  if (v == 1.) return;
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
      return;
    }
  }
  _prodMatVecInPlacePtr(x.data(), y.data(), transpose);
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

void AMatrix::resize(int nrows, int ncols)
{
  // Check if nothing to be done
  if (nrows == getNRows() && ncols == getNCols()) return;

  _deallocate();
  _setNRows(nrows);
  _setNCols(ncols);
  _allocate();
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

  if (!_isCompatible(y))
  {
    messerr("Matrices 'y' and 'this' are not compatible");
    return;
  }
  for (int irow = 0; irow < _nRows; irow++)
    for (int icol = 0; icol < _nCols; icol++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;
      _setValue(irow, icol, cx * _getValue(irow, icol) + cy * y.getValue(irow, icol));
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

  if (nm1 != nm2)
  {
    messerr("Matrices 'x' and 'y' should have matching dimensions");
    return;
  }

  if (!_checkLink(x->getNRows(), x->getNCols(), transposeX,
                  y->getNRows(), y->getNCols(), transposeY)) return;

  for (int irow = 0; irow < ni1; irow++)
  {
    for (int icol = 0; icol < ni2; icol++)
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
void AMatrix::prodNormMatMatInPlace(const AMatrix &a,
                                    const AMatrix &m,
                                    bool transpose)
{
  int n1 = (transpose) ? a.getNCols() : a.getNRows();
  int n2 = (transpose) ? a.getNRows() : a.getNCols();

  if (!_checkLink(a.getNRows(), a.getNCols(), transpose,
                  m.getNRows(), m.getNCols(), false,
                  a.getNRows(), a.getNCols(), !transpose)) return;

  for (int i = 0; i < n1; i++)
    for (int j = 0; j < n1; j++)
    {
      if (!_isPhysicallyPresent(i, j)) continue;
      double value = 0;
      for (int k = 0; k < n2; k++)
        for (int l = 0; l < n2; l++)
        {
          double a_ik = (transpose) ? a.getValue(k,i) : a.getValue(i,k);
          double a_lj = (transpose) ? a.getValue(l,j) : a.getValue(j,l);
          value += a_ik * m.getValue(k,l) * a_lj;
        }
      setValue(i, j, value);
    }
}

void AMatrix::prodNormMatInPlace(const AMatrix &a, const VectorDouble& vec, bool transpose)
{
  int n1 = (transpose) ? a.getNCols() : a.getNRows();
  int n2 = (transpose) ? a.getNRows() : a.getNCols();

  if (!_checkLink(getNRows(), getNCols(), a.getNRows(), a.getNCols(), transpose,
                 vec.size(), 1, false)) return;

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
    my_throw("The size of 'vec' must match the number of rows");
  for (int irow = 0; irow < _nRows; irow++)
    for (int icol = 0; icol < _nCols; icol++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;
      _setValue(irow, icol, _getValue(irow, icol) * vec[irow]);
    }
}

void AMatrix::divideRow(const VectorDouble& vec)
{
  if (_nRows != (int) vec.size())
    my_throw("The size of 'vec' must match the number of rows");
  for (int irow = 0; irow < _nRows; irow++)
    for (int icol = 0; icol < _nCols; icol++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;
      _setValue(irow, icol, _getValue(irow, icol) / vec[irow]);
    }
}

void AMatrix::multiplyColumn(const VectorDouble& vec)
{
  if (_nCols != (int) vec.size())
    my_throw("The size of 'vec' must match the number of columns");
  for (int irow = 0; irow < _nRows; irow++)
    for (int icol = 0; icol < _nCols; icol++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;
      _setValue(irow, icol, _getValue(irow, icol) * vec[icol]);
    }
}
void AMatrix::divideColumn(const VectorDouble& vec)
{
  if (_nCols != (int) vec.size())
    my_throw("The size of 'vec' must match the number of columns");
  for (int irow = 0; irow < _nRows; irow++)
    for (int icol = 0; icol < _nCols; icol++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;
      _setValue(irow, icol, _getValue(irow, icol) / vec[icol]);
    }
}

VectorDouble AMatrix::prodVecMat(const VectorDouble& x, bool transpose) const
{
  VectorDouble y;

  if (!_checkLink(x.size(), 1, false, getNRows(), getNCols(), transpose)) return y;

  int nval = (! transpose) ? _nCols : _nRows;
  y.resize(nval, 0.);
  prodVecMatInPlace(x, y, transpose);
  return y;
}

VectorDouble AMatrix::prodMatVec(const VectorDouble& x, bool transpose) const
{
  VectorDouble y;

  if (!_checkLink(getNRows(), getNCols(), transpose, x.size(), 1, false)) return y;

  int nval = (transpose) ? _nCols : _nRows;
  y.resize(nval, 0.);
  prodMatVecInPlace(x, y, transpose);
  return y;
}

double AMatrix::quadraticMatrix(const VectorDouble& x, const VectorDouble& y)
{
  if (!_checkLink(x.size(), 1, true,
                  getNRows(), getNCols(), false,
                  y.size(), 1, false)) return TEST;

  VectorDouble left(_nRows);
  prodMatVecInPlace(y, left, false);
  return VH::innerProduct(x, left);
}

int AMatrix::invert()
{
  if (! isSquare())
    my_throw("Invert method is restricted to Square matrices");
  return _invert();
}

int AMatrix::solve(const VectorDouble& b, VectorDouble& x) const
{
  if (! isSquare())
    my_throw("Invert method is limited to Square Matrices");
  if ((int) b.size() != _nRows || (int) x.size() != _nRows)
    my_throw("b' and 'x' should have the same dimension as the Matrix");
  return _solve(b, x);
}

String AMatrix::toString(const AStringFormat* /* strfmt*/) const
{
  std::stringstream sstr;

  sstr << "- Number of rows    = " <<  _nRows << std::endl;
  sstr << "- Number of columns = " <<  _nCols << std::endl;
  if (_flagEigen)
    sstr << "  (using Eigen Library)" << std::endl;
  else
    sstr << "  (not using Eigen Library)" << std::endl;

  bool flagSkipZero = false;
  if (isSparse())
  {
    sstr << "- Sparse Format" << std::endl;
    flagSkipZero = true;
  }
  sstr << toMatrix(String(), VectorString(), VectorString(), true, _nCols, _nRows,
                   getValues(), false, flagSkipZero);

  return sstr.str();
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
  if (! _flagCheckAddress) return true;
  if (irow < 0 || irow >= getNRows())
  {
    mesArg("Row index invalid",irow,getNRows());
    return false;
  }
  return true;
}

bool AMatrix::_isColumnValid(int icol) const
{
  if (! _flagCheckAddress) return true;
  if (icol < 0 || icol >= getNCols())
  {
    mesArg("Column index invalid",icol,getNCols());
    return false;
  }
  return true;
}

bool AMatrix::_isIndexValid(int irow, int icol) const
{
  if (! _flagCheckAddress) return true;
  if (! _isRowValid(irow)) return false;
  if (! _isColumnValid(icol)) return false;
  return true;
}

bool AMatrix::_isRowVectorConsistent(const VectorDouble& tab)
{
  if (tab.size() != (unsigned) getNRows())
  {
    messerr("Argument vector size should match the number of rows");
    return false;
  }
  return true;
}

bool AMatrix::_isColVectorConsistent(const VectorDouble& tab)
{
  if (tab.size() != (unsigned) getNCols())
  {
    messerr("Argument vector size should match the number of columns");
    return false;
  }
  return true;
}

bool AMatrix::_isVectorSizeConsistent(int nrows,
                                      int ncols,
                                      const VectorDouble &tab)
{
  if (tab.size() != (unsigned) (nrows * ncols))
  {
    messerr("The VectorDouble argument does not have correct dimension");
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
    int nr3 = (!transpose2) ? nrow3 : ncol3;
    int nc3 = (!transpose2) ? ncol3 : nrow3;
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
NF_Triplet AMatrix::getMatrixToTriplets(bool flag_from_1) const
{
  NF_Triplet NF_T = tripletInit(0);

  NF_T.flagFromOne = false;
  for (int icol = 0; icol < _nCols; icol++)
    for (int irow = 0; irow < _nRows; irow++)
    {
      if (!isValid(irow, icol)) continue;
      double value = getValue(irow, icol);
      if (ABS(value) < EPSILON10) continue;
      tripletAdd(NF_T, irow, icol, value);
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
        (*itvect) = getValue_(irow, icol);
        itvect++;
      }
  }
  else
  {
    for (int irow = 0; irow < _nRows; irow++)
      for (int icol = 0; icol < _nCols; icol++)
      {
        (*itvect) = getValue_(irow, icol);
        itvect++;
      }
  }
  return vect;
}

double AMatrix::compare(const AMatrix& mat) const
{
  if (mat.getNRows() != _nRows || mat.getNCols() != _nCols)
    my_throw("We can only compare two matrices with same dimensions");

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
    my_throw("This function is only valid for Square matrices");

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
 */
void AMatrix::setDiagonal(const VectorDouble& tab)
{
  if (! isSquare())
    my_throw("This function is only valid for Square matrices");

  fill(0.);
  for (int irow = 0; irow < getNRows(); irow++)
  {
    int icol = irow;
    if (icol < 0 || icol >= getNCols()) continue;
    setValue(irow,icol,tab[irow]);
  }
}

void AMatrix::setDiagonalToConstant(double value)
{
  if (! isSquare())
    my_throw("This function is only valid for Square matrices");

  fill(0.);
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
  if (irow < 0 || irow >= getNRows())
  {
    mesArg("Incorrect argument 'irow'",irow,getNRows());
    return vect;
  }

  for (int icol = 0; icol < getNCols(); icol++)
    vect.push_back(getValue(irow,icol));
  return vect;
}

/*! Set the contents of a Row */
void AMatrix::setRow(int irow, const VectorDouble& tab)
{
  if (irow < 0 || irow >= getNRows())
    my_throw("Incorrect argument 'irow'");
  if ((int) tab.size() != getNCols())
    my_throw("Incorrect dimension of 'tab'");

  for (int icol = 0; icol < getNCols(); icol++)
    setValue(irow,icol,tab[icol]);
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

/*! Set the contents of a Column */
void AMatrix::setColumn(int icol, const VectorDouble& tab)
{
  if (icol < 0 || icol >= getNCols())
    my_throw("Incorrect argument 'icol'");
  if ((int) tab.size() != getNRows())
    my_throw("Incorrect dimension of 'tab'");

  for (int irow = 0; irow < getNRows(); irow++)
    setValue(irow,icol,tab[irow]);
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
int AMatrix::getNumberColumnDefined() const
{
  int ncol = 0;
  for (int icol = 0; icol < getNCols(); icol++)
  {
    if (isColumnDefined(icol)) ncol++;
  }
  return ncol;
}

/*! Define the number of defined rows */
int AMatrix::getNumberRowDefined() const
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
  double oldval = _getValue(irow, icol);
  if (FFFF(oldval)) return;
  _setValue(irow, icol, oldval + value);
}

void AMatrix::_clear()
{
  _setNRows(0);
  _setNCols(0);
  _allocate();
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
      double value = _getValue(irow, icol);
      if (FFFF(value)) continue;
      if (value < minimum) minimum = value;
    }
  if (minimum == 1.e30) minimum = TEST;
  return minimum;
}

double AMatrix::getMaximum() const
{
  double maximum = -1.e30;
  for (int icol = 0; icol < getNCols(); icol++)
    for (int irow = 0; irow < getNRows(); irow++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;
      double value = _getValue(irow, icol);
      if (FFFF(value)) continue;
      if (value > maximum) maximum = value;
    }
  if (maximum == -1.e30) maximum = TEST;
  return maximum;
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
      setValue_(irow, icol, x->getValue_(validRows[irow], validCols[icol]));
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
      setValue_(irow, icol, factor * m.getValue_(irow, icol));
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

/**
 * Define the use of Eigen Library according to the value of input argulent 'opt_eigen'
 * @param opt_eigen Choice: 0: Do not use Eigen library; 1; Use Eigen library; -1: use global environment
 * @return Option for using the Eigen library
 */
bool AMatrix::_defineFlagEigen(int opt_eigen) const
{
  // Dispatch

  if (opt_eigen == 1)
    return true;
  else if (opt_eigen == 0)
    return false;
  else
    return globalFlagEigen;
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
bool AMatrix::isNonNegative(bool verbose)
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
 * Modify the parameter for using EIGEN library or not.
 * Warning: this must be performed very early in the script in order to forbid mixing two different styles.
 * @param flagEigen True if EIGEN library must be used; False otherwise (old style)
 */
void setGlobalFlagEigen(bool flagEigen)
{
  globalFlagEigen = flagEigen;
}

bool isGlobalFlagEigen()
{
  return globalFlagEigen;
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
