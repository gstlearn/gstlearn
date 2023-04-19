/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
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

// External library /// TODO : Dependency to csparse to be removed
#include "csparse_d.h"
#include "csparse_f.h"

AMatrix::AMatrix(int nrow, int ncol, bool sparse)
    : AStringable(),
      _nRows(nrow),
      _nCols(ncol),
      _sparse(sparse),
      _csMatrix(nullptr)
{
  (void) _isNumbersValid(nrow, ncol);
  if (sparse) _initiateSparse();
}
#ifndef SWIG
AMatrix::AMatrix(const cs* A)
    : AStringable(),
      _nRows(0),
      _nCols(0),
      _sparse(true),
      _csMatrix(nullptr)
{
  _recopySparse(A);
}
#endif
AMatrix::AMatrix(const AMatrix &m)
    : AStringable(m),
      _nRows(m._nRows),
      _nCols(m._nCols),
      _sparse(m._sparse),
      _csMatrix(nullptr)
{
  if (_sparse)
  {
    _recopySparse(m._csMatrix);
  }
}

AMatrix& AMatrix::operator=(const AMatrix &m)
{
  AStringable::operator=(m);
  _nRows = m._nRows;
  _nCols = m._nCols;
  _sparse = m._sparse;
  if (_sparse)
  {
    if (this != &m)
    {
      _deallocateSparse();
      _recopySparse(m._csMatrix);
    }
  }
  return *this;
}

AMatrix::~AMatrix()
{
  if (_sparse)
  {
    _deallocateSparse();
  }
}

void AMatrix::init(int nrows, int ncols, bool sparse)
{
  if (_sparse)
  {
    _deallocateSparse();
  }
  _nRows = nrows;
  _nCols = ncols;
  _sparse = sparse;
  _allocate();
}

bool AMatrix::isSquare(bool printWhyNot) const
{
  if (isEmpty()) return false;
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

void AMatrix::reset(int nrows, int ncols, bool sparse)
{
  if (! _isNumbersValid(nrows, ncols)) return;
  _sparse = sparse;
  if (_sparse)
  {
    _nRows = nrows;
    _nCols = ncols;
  }
  else
  {
    _deallocate();
    _nRows = nrows;
    _nCols = ncols;
    _allocate();
  }
  _clearContents();
}

void AMatrix::reset(int nrows, int ncols, double value, bool sparse)
{
  if (! _isNumbersValid(nrows, ncols)) return;
  _sparse = sparse;
  if (_sparse)
  {
    _forbiddenForSparse("reset");
  }
  else
  {
    _deallocate();
    _nRows = nrows;
    _nCols = ncols;
    _allocate();
    for (int i=0; i<_getMatrixSize(); i++)
      _setValue(i, value);
  }
  _clearContents();
}

void AMatrix::reset(int nrows, int ncols, const double* tab, bool sparse)
{
  if (! _isNumbersValid(nrows, ncols)) return;
  _sparse = sparse;
  if (_sparse)
  {
    _forbiddenForSparse("reset");
  }
  else
  {
    _deallocate();
    _nRows = nrows;
    _nCols = ncols;
    _allocate();
    int lec = 0;
    for (int icol=0; icol<ncols; icol++)
      for (int irow=0; irow<nrows; irow++)
        _setValue(irow,icol,tab[lec++]);
  }
  _clearContents();
}

void AMatrix::reset(int nrows, int ncols, const VectorDouble& tab, bool sparse, bool flagByRow)
{
  if (! _isNumbersValid(nrows, ncols)) return;
  _nRows = nrows;
  _nCols = ncols;
  _sparse = sparse;
  if (_sparse)
  {
    _forbiddenForSparse("reset");
  }
  else
  {
    _allocate();
    int lec = 0;
    if (flagByRow)
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
  }
  _clearContents();
}

void AMatrix::reset(const VectorVectorDouble& tab, bool flagByRow)
{
  if (flagByRow)
  {
    _nRows = (int) tab.size();
    _nCols = (int) tab[0].size();
    _allocate();
    for (int icol = 0; icol < _nCols; icol++)
      for (int irow = 0; irow < _nRows; irow++)
        _setValue(irow, icol, tab[irow][icol]);
  }
  else
  {
    _nCols = (int) tab.size();
    _nRows = (int) tab[0].size();
    _allocate();
    for (int icol = 0; icol < _nCols; icol++)
      for (int irow = 0; irow < _nRows; irow++)
        _setValue(irow, icol, tab[icol][irow]);
  }
  _clearContents();
}

void AMatrix::fillRandom(int seed, double zeroPercent)
{
  law_set_random_seed(seed);

  double value = 0.;
  if (_sparse)
  {
    cs *Atriplet;
    Atriplet = cs_spalloc(0, 0, 1, 1, 1);
    for (int irow = 0; irow < _nRows; irow++)
      for (int icol = 0; icol < _nCols; icol++)
      {
        if (! _isPhysicallyPresent(irow, icol)) continue;
        if (!mustBeDiagCst() && law_uniform(0., 1.) < zeroPercent) continue;
        value = law_gaussian();
        cs_entry(Atriplet, irow, icol, value);
      }
    _csMatrix = cs_triplet(Atriplet);
    Atriplet = cs_spfree(Atriplet);
  }
  else
  {
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
}

bool AMatrix::isSymmetric(bool printWhyNot) const
{
  if (isEmpty() || ! isSquare()) return false;

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

bool AMatrix::isDiagonal(bool printWhyNot) const
{
  if (isEmpty() || ! isSquare()) return false;

  for (int irow = 0; irow <_nRows; irow++)
    for (int icol = 0; icol < _nCols; icol++)
    {
      if (irow != icol)
      {
        if (ABS(getValue(irow,icol)) > EPSILON10)
        {
          if (printWhyNot)
            messerr("The element (%d;%d)=%lf should be zero",
                    irow,icol,getValue(irow,icol));
          return false;
        }
      }
    }
  return true;
}

bool AMatrix::isDiagCst(bool printWhyNot) const
{
  if (isEmpty() || ! isSquare()) return false;

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
  if (_sparse)
  {
    cs* old = _csMatrix;
    _csMatrix = cs_transpose(old, 0);
    cs_spfree(old);
  }
  else
  {
    _transposeInPlace();
  }
}

AMatrix* AMatrix::transpose() const
{
  AMatrix* mat = dynamic_cast<AMatrix*>(clone());
  mat->transposeInPlace();
  return mat;
}

AMatrix* AMatrix::toSparse() const
{
  AMatrix* mat = dynamic_cast<AMatrix*>(clone());

  // Create the triplet structure from the non-zero terms of source matrix

  VectorInt irows;
  VectorInt icols;
  VectorDouble values;
  getValuesAsTriplets(irows, icols, values);

  // Update the contents of the cloned matrix as sparse

  mat->reset(_nRows, _nCols, true);

  // Load the triplet information in the cloned matrix

  mat->setValuesByArrays(irows, icols, values);

  return mat;
}

void AMatrix::toSparseInPlace()
{
  // Create the triplet structure from the non-zero terms of source matrix

  VectorInt irows;
  VectorInt icols;
  VectorDouble    values;
  getValuesAsTriplets(irows, icols, values);

  // Update the sparse qualifier of this

  _deallocateSparse();
  _setSparse(true);

  // Load the triplet information in the cloned matrix

  setValuesByArrays(irows, icols, values);
}

/*! Gets the value at row 'irow' and column 'icol' */
double AMatrix::getValue(int irow, int icol) const
{
  if (! _isIndexValid(irow, icol)) return TEST;
  if (_sparse)
  {
    return cs_get_value(_csMatrix, irow, icol);
  }
  else
  {
    return _getValue(irow, icol);
  }
}

/*! Sets the value at row 'irow' and column 'icol' */
void AMatrix::setValue(int irow, int icol, double value)
{
  if (! _isIndexValid(irow, icol)) return;
  if (_sparse)
  {
    cs_set_value(_csMatrix, irow, icol, value);
  }
  else
  {
    return _setValue(irow, icol, value);
  }
}

/*! Gets the value at rank 'rank' */
double AMatrix::getValue(int rank) const
{
  if (! _isRankValid(rank)) return TEST;
  if (_sparse)
  {
    _forbiddenForSparse("getValue");
  }
  else
  {
    return _getValue(rank);
  }
  return 0.;
}

/*! Sets the value at rank 'rank' */
void AMatrix::setValue(int rank, double value)
{
  if (! _isRankValid(rank)) return;
  if (_sparse)
  {
    _forbiddenForSparse("getValue");
  }
  else
  {
    _setValue(rank, value);
  }
}

/*! Gets a reference to the value at row 'irow' and column 'icol' */
double& AMatrix::getValueRef(int irow, int icol)
{
  if (_sparse)
  {
    _forbiddenForSparse("getValueRef");
  }
  return _getValueRef(irow, icol);
}

/**
 * Fill 'this' with the constant 'value'
 * @param value Constant value used for filling 'this'
 */
void AMatrix::fill(double value)
{
  if (_sparse)
  {
    cs* Mtriplet = cs_spalloc(0,0,1,1,1);
    for (int icol = 0; icol < _nCols; icol++)
      for (int irow = 0; irow < _nRows; irow++)
      (void) cs_entry(Mtriplet, icol, irow, value);
    _csMatrix = cs_triplet(Mtriplet);
    Mtriplet = cs_spfree(Mtriplet);
  }
  else
  {
    for (int rank = 0; rank < _getMatrixSize(); rank++)
    {
      _setValue(rank, value);
    }
  }
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
void AMatrix::setValuesOldStyle(const double* values, bool byCol)
{
  if (_sparse)
  {
    cs* Mtriplet = cs_spalloc(0, 0, 1, 1, 1);
    if (byCol)
    {
      int lec = 0;
      for (int icol = 0; icol < _nCols; icol++)
        for (int irow = 0; irow < _nRows; irow++, lec++)
        {
          if (ABS(values[lec]) < EPSILON10) continue;
          (void) cs_entry(Mtriplet, irow, icol, values[lec]);
        }
    }
    else
    {
      int lec = 0;
      for (int irow = 0; irow < _nRows; irow++)
        for (int icol = 0; icol < _nCols; icol++, lec++)
        {
          if (ABS(values[lec]) < EPSILON10) continue;
          (void) cs_entry(Mtriplet, irow, icol, values[lec]);
        }
    }
    _csMatrix = cs_triplet(Mtriplet);
    Mtriplet = cs_spfree(Mtriplet);
  }
  else
  {
    _setValues(values, byCol);
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
  if ((int) values.size() != getNTotal())
  {
    messerr("Inconsistency between 'values' and Matrix Dimension");
    messerr("Operation cancelled");
    return;
  }
  setValuesOldStyle(values.data(),byCol);
}

void AMatrix::setValuesByArrays(const VectorInt &irows,
                                const VectorInt &icols,
                                const VectorDouble &values)
{
  int nelements = static_cast<int> (values.size());
  if (irows.size() != values.size() ||
      icols.size() != values.size())
  {
    messerr("The arguments 'icols', 'irows' and 'values' should share same positive dimension");
    messerr("Operation cancelled");
    return;
  }

  if (_sparse)
  {
    cs* Mtriplet = cs_spalloc(0,0,1,1,1);
    for (int i = 0; i < nelements; i++)
      (void) cs_entry(Mtriplet, irows[i], icols[i], values[i]);
    _csMatrix = cs_triplet(Mtriplet);
    Mtriplet = cs_spfree(Mtriplet);
  }
  else
  {

    for (int i = 0; i < nelements; i++)
      setValue(irows[i], icols[i], values[i]);
  }
}

void AMatrix::setIdentity(double value)
{
  for (int icol = 0; icol < _nCols; icol++)
    for (int irow = 0; irow < _nRows; irow++)
      setValue(irow, icol, value * (irow == icol));
}

/**
 *
 * @param v Add a scalar value to all terms of the current matrix
 */
void AMatrix::addScalar(double v)
{
  if (v == 0.) return;
  if (_sparse)
  {
    _forbiddenForSparse("addScalar");
  }
  else
  {
    for (int rank = 0; rank < _getMatrixSize(); rank++)
    {
      _setValue(rank, getValue(rank) + v);
    }
  }
}

/**
 *
 * @param v Add constant value to the diagonal of the current Matrix
 */
void AMatrix::addScalarDiag(double v)
{
  if (v == 0.) return;
  if (_sparse)
  {
    cs* csi = cs_eye(getNRows(), 1.);
    cs* res = cs_add(_csMatrix, csi, 1., v);
    cs_spfree(csi);
    cs_spfree(_csMatrix);
    _csMatrix = res;
  }
  else
  {
    for (int irow = 0; irow < _nRows; irow++)
    {
      for (int icol = 0; icol < _nCols; icol++)
      {
        if (irow == icol)
        {
          setValue(irow, icol, getValue(irow, icol) + v);
        }
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
  if (_sparse)
  {
    cs* res = cs_add(_csMatrix, _csMatrix, v, 0.);
    cs_spfree(_csMatrix);
    _csMatrix = res;
  }
  else
  {
    for (int rank = 0; rank < _getMatrixSize(); rank++)
    {
      _setValue(rank, getValue(rank) * v);
    }
  }
}

/**
 *
 * @param inv Input vector
 * @param outv Output vector obtained by multiplying 'inv' by current Matrix
 */
#ifndef SWIG
void AMatrix::prodVector(const double *inv, double *outv) const
{
  if (_sparse)
  {
    cs_vecmult(_csMatrix, _nRows, inv, outv);
  }
  else
  {
    _prodVector(inv, outv);
  }
}
#endif

void AMatrix::prodVector(const VectorDouble& inv, VectorDouble& outv) const
{
  // TODO : Check dimensions to avois SEGV
  prodVector(inv.data(), outv.data());
}

/**
 * Add the matrix 'y' to the current Matrix
 * @param y Matrix to be added
 */
void AMatrix::addMatrix(const AMatrix& y)
{
  if (! isSameSize(y))
  {
    messerr("Matrices 'y' and 'this' should have the same size");
    return;
  }

  if (_sparse)
  {
    cs* res = cs_add(_csMatrix, y._csMatrix, 1., 1.);
    cs_spfree(_csMatrix);
    _csMatrix = res;
  }
  else
  {
    if (! _isCompatible(y))
    {
      messerr("Matrices 'y' and 'this' are not compatible");
      return;
    }
    for (int irow = 0; irow < _nRows; irow++)
    {
      for (int icol = 0; icol < _nCols; icol++)
      {
        if (!_isPhysicallyPresent(irow, icol)) continue;
        setValue(irow, icol, getValue(irow, icol) + y.getValue(irow, icol));
      }
    }
  }
}

/**
 * Store the product of 'x' by 'y' in this
 * @param x First Matrix
 * @param y Second matrix
 */
void AMatrix::prodMatrix(const AMatrix& x, const AMatrix& y)
{
  if (x.getNCols() != y.getNRows() ||
      x.getNRows() != getNRows() ||
      y.getNCols() != getNCols())
  {
    messerr("Incompatible matrix dimensions for matrix product");
    messerr("- First matrix:  NRows = %d - NColumns = %d", x.getNRows(), x.getNCols());
    messerr("- Second matrix: NRows = %d - NColumns = %d", y.getNRows(), y.getNCols());
    messerr("- Result matrix: NRows = %d - NColumns = %d", getNRows(), getNCols());
    messerr("Operation is cancelled");
    return;
  }

  if (_sparse)
  {
    cs* res = cs_multiply(x._csMatrix, y._csMatrix);
    cs_spfree(_csMatrix);
    _csMatrix = res;
  }
  else
  {
    int n = x.getNCols();
    for (int irow = 0; irow < _nRows; irow++)
    {
      for (int icol = 0; icol < _nCols; icol++)
      {
        if (!_isPhysicallyPresent(irow, icol)) continue;

        double value = 0.;
        for (int k = 0; k < n; k++)
        {
          value += x.getValue(irow, k) * y.getValue(k, icol);
        }
        setValue(irow, icol, value);
      }
    }
  }
}

/*!
 * Updates the current Matrix as a linear combination of matrices as follows:
 *  this <- cx * this + cy * y
 * @param cx Coefficient applied to the current Matrix
 * @param cy Coefficient applied to the Matrix  'y'
 * @param y Second Matrix in the Linear combination
 */
void AMatrix::linearCombination(double cx, double cy, const AMatrix& y)
{
  if (! isSameSize(y))
  {
    my_throw("Matrices should have same size");
  }

  if (_sparse)
  {
    if (!y._sparse)
      my_throw("This function can only combine sparse matrices together");
    cs* res = cs_add(_csMatrix, y._csMatrix, cx, cy);
    cs_spfree(_csMatrix);
    _csMatrix = res;
  }
  else
  {
    if (! _isCompatible(y))
      my_throw("Matrix 'y' is not compatible with 'this'");
    for (int rank = 0; rank < _getMatrixSize(); rank++)
    {
      double value = getValue(rank) * cx + cy * y.getValue(rank);
      _setValue(rank, value);
    }
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

int AMatrix::invert()
{
  if (! isSquare())
    my_throw("Invert method is restricted to Square matrices");
  if (_sparse)
  {
    cs *inv = cs_invert(_csMatrix,0);
    _deallocateSparse();
    _csMatrix = inv;
  }
  else
  {
    return _invert();
  }
  return 0;
}

int AMatrix::solve(const VectorDouble& b, VectorDouble& x) const
{
  int error = 0;

  if (! isSymmetric())
    my_throw("Invert method is limited to Square Symmetrical Matrices");
  if ((int) b.size() != _nRows || (int) x.size() != _nRows)
    my_throw("b' and 'x' should have the same dimension as the Matrix");

  if (_sparse)
  {
    x = b;
    error = cs_cholsol(_csMatrix,x.data(), 0);
  }
  else
  {
    error = _solve(b, x);
  }
  return error;
}

String AMatrix::toString(const AStringFormat* /* strfmt*/) const
{
  std::stringstream sstr;

  sstr << "- Number of rows    = " <<  _nRows << std::endl;
  sstr << "- Number of columns = " <<  _nCols << std::endl;

   if (_sparse)
   {
     sstr << "- Sparse Format" << std::endl;
     sstr << toMatrix(String(), _csMatrix);
   }
   else
   {
     sstr << toMatrix(String(), VectorString(), VectorString(), true, _nCols, _nRows,
                     getValues());
   }
  return sstr.str();
}

bool AMatrix::_isNumbersValid(int nrows, int ncols) const
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

bool AMatrix::_isRowValid(int irow) const
{
  if (irow < 0 || irow >= getNRows())
  {
    mesArg("Row index invalid",irow,getNRows());
    return false;
  }
  return true;
}

bool AMatrix::_isColumnValid(int icol) const
{
  if (icol < 0 || icol >= getNCols())
  {
    mesArg("Column index invalid",icol,getNCols());
    return false;
  }
  return true;
}

bool AMatrix::_isIndexValid(int irow, int icol) const
{
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
                                       const VectorDouble& tab)
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
  return (rank >= 0 && rank < _getMatrixSize());
}

/**
 * This strange function instantiate a sparse matrix with given dimensions
 * filled with zeroes. It should be an empty matrix... But this does not make sense.
 * Therefore it is created by setting a single element at the lower bottom size of
 * the matrix ... filled with a zero.
 */
void AMatrix::_initiateSparse()
{
  cs *Atriplet;
  Atriplet = cs_spalloc(0, 0, 1, 1, 1);
  cs_entry(Atriplet, getNRows() - 1, getNCols() - 1, 0.);
  _csMatrix = cs_triplet(Atriplet);
  Atriplet = cs_spfree(Atriplet);
}

void AMatrix::_recopySparse(const cs* cs)
{
  int nrows, ncols, count;
  double percent;

  cs_rowcol(cs, &nrows, &ncols, &count, &percent);
  _setNRows(nrows);
  _setNCols(ncols);

  _csMatrix = cs_duplicate(cs);
}

void AMatrix::_deallocateSparse()
{
  _csMatrix = cs_spfree(_csMatrix);
}

void AMatrix::_forbiddenForSparse(const String& func) const
{
  messerr("Problem with Function: %s",func.c_str());
  my_throw("This function is not available in Sparse Matrix");
}

AMatrix* createIdentity(int nrow, bool sparse)
{
  return MatrixFactory::createIdentity(nrow, sparse);
}

AMatrix* transpose(const AMatrix* mat)
{
  return mat->transpose();
}

AMatrix* prodMatrix(const AMatrix *mat1, const AMatrix *mat2)
{
  return MatrixFactory::matProduct(mat1, mat2);
}

void prodMatrixInPlace(AMatrix* mat1, const AMatrix* mat2)
{
  AMatrix* res = prodMatrix(mat1, mat2);
  for (int i = 0; i < res->getNTotal(); i++)
    mat1->setValue(i, res->getValue(i));
  delete res;
}

void AMatrix::dumpElements(const String& title, int ifrom, int ito) const
{
  if (_sparse)
    messerr("This method is not implemented for Sparse Matrix");
  else
  {
    mestitle(1, "%s", title.c_str());
    for (int rank=ifrom; rank<ito; rank++)
    {
      if (_isRankValid(rank))
        message("Element %d = %lf\n",rank,getValue(rank));
    }
  }
}

/**
 * From a matrix of any type, creates the three vectors of the triplet
 * (specific format for creating efficiently a Sparse matrix)
 * It only takes the only non-zero elements of the matrix
 * @param irows Output array of row indices
 * @param icols Output array of column indices
 * @param values Output array of non-zero values
 */
void AMatrix::getValuesAsTriplets(VectorInt&    irows,
                                  VectorInt&    icols,
                                  VectorDouble& values) const
{
  /// TODO : use cs_sparce corresponding function
  for (int icol = 0; icol < _nCols; icol++)
    for (int irow = 0; irow < _nRows; irow++)
    {
      if (!isValid(irow, icol)) continue;
      double value = getValue(irow, icol);
      if (ABS(value) < EPSILON10) continue;
      irows.push_back(irow);
      icols.push_back(icol);
      values.push_back(value);
    }
}

VectorDouble AMatrix::getValues() const
{
  VectorDouble vect;
  for (int icol = 0; icol < _nCols; icol++)
    for (int irow = 0; irow < _nRows; irow++)
    {
      double value = 0.;
      if (isValid(irow, icol)) value = getValue(irow,icol);
      vect.push_back(value);
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

void AMatrix::setDiagonal(double value)
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

/**
 * Perform: this += 'value' * 'tab'
 * @param tab    Current matrix
 * @param value  Multiplicative coefficient  (default = 1)
 */
void AMatrix::add(const AMatrix& tab, double value)
{
  if (! isSameSize(tab))
    my_throw("Can only add matrices of same dimensions");

  for (int icol = 0; icol < getNCols(); icol++)
    for (int irow = 0; irow < getNRows(); irow++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;
      double oldval = getValue(irow, icol);
      setValue(irow, icol, oldval + value * tab.getValue(irow, icol));
    }
}

void AMatrix::add(int irow, int icol, double value)
{
  double oldval = getValue(irow, icol);
  if (FFFF(oldval)) return;
  setValue(irow, icol, oldval + value);
}

void AMatrix::subtract(const AMatrix& tab, double value)
{
  if (! isSameSize(tab))
    my_throw("Can only add matrices of same dimensions");

  for (int icol = 0; icol < getNCols(); icol++)
    for (int irow = 0; irow < getNRows(); irow++)
    {
      if (!_isPhysicallyPresent(irow, icol)) continue;
      double oldval = getValue(irow, icol);
      setValue(irow, icol, oldval - value * tab.getValue(irow, icol));
    }
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
  for (int i = 0; i < getNTotal(); i++)
  {
    double value = _getValue(i);
    if (FFFF(value)) continue;
    if (value < minimum) minimum = value;
  }
  if (minimum == 1.e30) minimum = TEST;
  return minimum;
}

double AMatrix::getMaximum() const
{
  double maximum = -1.e30;
  for (int i = 0; i < getNTotal(); i++)
  {
    double value = _getValue(i);
    if (FFFF(value)) continue;
    if (value > maximum) maximum = value;
  }
  if (maximum == -1.e30) maximum = TEST;
  return maximum;
}

void AMatrix::copyReduce(const AMatrix *x,
                         const VectorInt &validRows,
                         const VectorInt &validCols)
{
  VH::display("copyreduce validRows",validRows);
  for (int irow = 0; irow < (int) validRows.size(); irow++)
    for (int icol = 0; icol < (int) validCols.size(); icol++)
      setValue(irow, icol, x->getValue(validRows[irow], validCols[icol]));
}

Triplet AMatrix::getCsToTriplet(bool flag_from_1) const
{
  return csToTriplet(getCs(), flag_from_1);
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

