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
#include "Matrix/MatrixSparse.hpp"
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

MatrixSparse::MatrixSparse(int nrow, int ncol)
    : AMatrix(nrow, ncol),
      _csMatrix(nullptr)
{
  _allocate();
}

#ifndef SWIG
MatrixSparse::MatrixSparse(const cs *A)
    : AMatrix(0, 0),
      _csMatrix(nullptr)
{
  _recopySparse(A);
}
#endif

MatrixSparse::MatrixSparse(const MatrixSparse &m)
    : AMatrix(m),
      _csMatrix(nullptr)
{
  _recopySparse(m._csMatrix);
}

MatrixSparse& MatrixSparse::operator=(const MatrixSparse &m)
{
  if (this != &m)
  {
    AMatrix::operator=(m);
    _deallocate();
    _recopySparse(m._csMatrix);
  }
  return *this;
}

MatrixSparse::~MatrixSparse()
{
  _deallocate();
}

void MatrixSparse::init(int nrows, int ncols)
{
  _deallocate();
  _setNRows(nrows);
  _setNCols(ncols);
  _allocate();
}

void MatrixSparse::reset(int nrows, int ncols)
{
  if (! _isNumbersValid(nrows, ncols)) return;
  _setNRows(nrows);
  _setNCols(ncols);
  _clearContents();
}

void MatrixSparse::reset(int nrows, int ncols, double value)
{
  if (! _isNumbersValid(nrows, ncols)) return;
  _forbiddenForSparse("reset");
}

void MatrixSparse::reset(int nrows, int ncols, const double* tab, bool byCol)
{
  if (! _isNumbersValid(nrows, ncols)) return;
  _forbiddenForSparse("reset");
}

void MatrixSparse::reset(int nrows, int ncols, const VectorDouble& tab, bool byCol)
{
  if (! _isNumbersValid(nrows, ncols)) return;
  reset(nrows, ncols, tab.data(), byCol);
}

void MatrixSparse::reset(const VectorVectorDouble& tab, bool byCol)
{
  if (byCol)
  {
    _setNRows((int) tab.size());
    _setNCols((int) tab[0].size());
    _allocate();
    for (int icol = 0; icol < getNCols(); icol++)
      for (int irow = 0; irow < getNRows(); irow++)
        _setValue(irow, icol, tab[irow][icol]);
  }
  else
  {
    _setNCols((int) tab.size());
    _setNRows((int) tab[0].size());
    _allocate();
    for (int icol = 0; icol < getNCols(); icol++)
      for (int irow = 0; irow < getNRows(); irow++)
        _setValue(irow, icol, tab[icol][irow]);
  }
  _clearContents();
}

void MatrixSparse::fillRandom(int seed, double zeroPercent)
{
  law_set_random_seed(seed);

  double value = 0.;
  cs *Atriplet;
  Atriplet = cs_spalloc(0, 0, 1, 1, 1);
  for (int irow = 0; irow < getNRows(); irow++)
    for (int icol = 0; icol < getNCols(); icol++)
    {
      if (! _isPhysicallyPresent(irow, icol)) continue;
      if (!mustBeDiagCst() && law_uniform(0., 1.) < zeroPercent) continue;
      value = law_gaussian();
      cs_entry(Atriplet, irow, icol, value);
    }
  _csMatrix = cs_triplet(Atriplet);
  Atriplet = cs_spfree(Atriplet);
}

void MatrixSparse::_transposeInPlace()
{
  cs* old = _csMatrix;
  _csMatrix = cs_transpose(old, 0);
  cs_spfree(old);
}

MatrixSparse* MatrixSparse::transpose() const
{
  MatrixSparse* mat = dynamic_cast<MatrixSparse*>(clone());
  mat->transposeInPlace();
  return mat;
}

/*! Gets the value for rank 'rank' */
double MatrixSparse::_getValue(int rank) const
{
  _forbiddenForSparse("_setValue (by rank)");
  return TEST;
}

/*! Gets the value at row 'irow' and column 'icol' */
double MatrixSparse::_getValue(int irow, int icol) const
{
  if (! _isIndexValid(irow, icol)) return TEST;
  return cs_get_value(_csMatrix, irow, icol);
}

void MatrixSparse::_setValue(int rank, double value)
{
  _forbiddenForSparse("_setValue (by rank)");
}

/*! Sets the value at row 'irow' and column 'icol' */
void MatrixSparse::_setValue(int irow, int icol, double value)
{
  if (! _isIndexValid(irow, icol)) return;
  cs_set_value(_csMatrix, irow, icol, value);
}

/*! Gets a reference to the value at row 'irow' and column 'icol' */
double& MatrixSparse::_getValueRef(int irow, int icol)
{
  _forbiddenForSparse("_getValueRef");
  return _csMatrix->x[0]; // This is never used
}

int MatrixSparse::_getMatrixSize() const
{
  return _csMatrix->nz;
}

/**
 * Fill 'this' with the constant 'value'
 * @param value Constant value used for filling 'this'
 */
void MatrixSparse::fill(double value)
{
  cs* Mtriplet = cs_spalloc(0,0,1,1,1);
  for (int icol = 0; icol < getNCols(); icol++)
    for (int irow = 0; irow < getNRows(); irow++)
      (void) cs_entry(Mtriplet, icol, irow, value);
  _csMatrix = cs_triplet(Mtriplet);
  Mtriplet = cs_spfree(Mtriplet);
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
void MatrixSparse::_setValues(const double* values, bool byCol)
{
  cs* Mtriplet = cs_spalloc(0, 0, 1, 1, 1);
  if (byCol)
  {
    int lec = 0;
    for (int icol = 0; icol < getNCols(); icol++)
      for (int irow = 0; irow < getNRows(); irow++, lec++)
      {
        if (ABS(values[lec]) < EPSILON10) continue;
        (void) cs_entry(Mtriplet, irow, icol, values[lec]);
      }
  }
  else
  {
    int lec = 0;
    for (int irow = 0; irow < getNRows(); irow++)
      for (int icol = 0; icol < getNCols(); icol++, lec++)
      {
        if (ABS(values[lec]) < EPSILON10) continue;
        (void) cs_entry(Mtriplet, irow, icol, values[lec]);
      }
  }
  _csMatrix = cs_triplet(Mtriplet);
  Mtriplet = cs_spfree(Mtriplet);
}
#endif

void MatrixSparse::setValuesByArrays(const VectorInt &irows,
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

  cs* Mtriplet = cs_spalloc(0,0,1,1,1);
  for (int i = 0; i < nelements; i++)
    (void) cs_entry(Mtriplet, irows[i], icols[i], values[i]);
  _csMatrix = cs_triplet(Mtriplet);
  Mtriplet = cs_spfree(Mtriplet);
}

/**
 *
 * @param v Add a scalar value to all terms of the current matrix
 */
void MatrixSparse::addScalar(double v)
{
  if (v == 0.) return;
  _forbiddenForSparse("addScalar");
}

/**
 *
 * @param v Add constant value to the diagonal of the current Matrix
 */
void MatrixSparse::addScalarDiag(double v)
{
  if (v == 0.) return;
  cs* csi = cs_eye(getNRows(), 1.);
  cs* res = cs_add(_csMatrix, csi, 1., v);
  cs_spfree(csi);
  cs_spfree(_csMatrix);
  _csMatrix = res;
}

/**
 *
 * @param v Multiply all the terms of the matrix by the scalar 'v'
 */
void MatrixSparse::prodScalar(double v)
{
  if (v == 1.) return;
  cs* res = cs_add(_csMatrix, _csMatrix, v, 0.);
  cs_spfree(_csMatrix);
  _csMatrix = res;
}

/**
 *
 * @param inv Input vector
 * @param outv Output vector obtained by multiplying 'inv' by current Matrix
 */
#ifndef SWIG
void MatrixSparse::_prodVector(const double *inv, double *outv) const
{
  cs_vecmult(_csMatrix, getNRows(), inv, outv);
}
#endif

/**
 * Add the matrix 'y' to the current Matrix
 * @param y Matrix to be added
 */
void MatrixSparse::addMatrix(const MatrixSparse& y)
{
  if (! isSameSize(y))
  {
    messerr("Matrices 'y' and 'this' should have the same size");
    return;
  }

  cs* res = cs_add(_csMatrix, y._csMatrix, 1., 1.);
  cs_spfree(_csMatrix);
  _csMatrix = res;
}

/**
 * Store the product of 'x' by 'y' in this
 * @param x First Matrix
 * @param y Second matrix
 */
void MatrixSparse::prodMatrix(const MatrixSparse& x, const MatrixSparse& y)
{
  if (_getFlagCheckAddress())
  {
    if (x.getNCols() != y.getNRows() ||
        x.getNRows() != getNRows()   ||
        y.getNCols() != getNCols())
    {
      messerr("Incompatible matrix dimensions for matrix product");
      messerr("- First matrix:  NRows = %d - NColumns = %d", x.getNRows(),
              x.getNCols());
      messerr("- Second matrix: NRows = %d - NColumns = %d", y.getNRows(),
              y.getNCols());
      messerr("- Result matrix: NRows = %d - NColumns = %d", getNRows(),
              getNCols());
      messerr("Operation is cancelled");
      return;
    }
  }

  cs* res = cs_multiply(x._csMatrix, y._csMatrix);
  cs_spfree(_csMatrix);
  _csMatrix = res;
}

/*!
 * Updates the current Matrix as a linear combination of matrices as follows:
 *  this <- cx * this + cy * y
 * @param cx Coefficient applied to the current Matrix
 * @param cy Coefficient applied to the Matrix  'y'
 * @param y Second Matrix in the Linear combination
 */
void MatrixSparse::linearCombination(double cx, double cy, const MatrixSparse& y)
{
  if (! isSameSize(y))
  {
    my_throw("Matrices should have same size");
  }

  if (!y.isSparse())
    my_throw("This function can only combine sparse matrices together");
  cs* res = cs_add(_csMatrix, y._csMatrix, cx, cy);
  cs_spfree(_csMatrix);
  _csMatrix = res;
}

int MatrixSparse::_invert()
{
  if (! isSquare())
    my_throw("Invert method is restricted to Square matrices");
  cs *inv = cs_invert(_csMatrix,0);
  _deallocate();
  _csMatrix = inv;
  return 0;
}

int MatrixSparse::_solve(const VectorDouble& b, VectorDouble& x) const
{
  int error = 0;

  if (! isSymmetric())
    my_throw("Invert method is limited to Square Symmetrical Matrices");
  if ((int) b.size() != getNRows() || (int) x.size() != getNRows())
    my_throw("b' and 'x' should have the same dimension as the Matrix");

  x = b;
  error = cs_cholsol(_csMatrix,x.data(), 0);
  return error;
}

String MatrixSparse::toString(const AStringFormat* /* strfmt*/) const
{
  std::stringstream sstr;

  sstr << "- Number of rows    = " <<  getNRows() << std::endl;
  sstr << "- Number of columns = " <<  getNCols() << std::endl;

  sstr << "- Sparse Format" << std::endl;
  sstr << toMatrix(String(), _csMatrix);

  return sstr.str();
}

/**
 * This strange function instantiate a sparse matrix with given dimensions
 * filled with zeroes. It should be an empty matrix... But this does not make sense.
 * Therefore it is created by setting a single element at the lower bottom size of
 * the matrix ... filled with a zero.
 */
void MatrixSparse::_allocate()
{
  cs *Atriplet;
  Atriplet = cs_spalloc(0, 0, 1, 1, 1);
  cs_entry(Atriplet, getNRows() - 1, getNCols() - 1, 0.);
  _csMatrix = cs_triplet(Atriplet);
  Atriplet = cs_spfree(Atriplet);
}

void MatrixSparse::_recopySparse(const cs* cs)
{
  int nrows, ncols, count;
  double percent;

  cs_rowcol(cs, &nrows, &ncols, &count, &percent);
  _setNRows(nrows);
  _setNCols(ncols);

  _csMatrix = cs_duplicate(cs);
}

void MatrixSparse::_deallocate()
{
  _csMatrix = cs_spfree(_csMatrix);
}

void MatrixSparse::_forbiddenForSparse(const String& func) const
{
  messerr("Problem with Function: %s",func.c_str());
  my_throw("This function is not available in Sparse Matrix");
}

void MatrixSparse::dumpElements(const String& title, int ifrom, int ito) const
{
  messerr("This method is not implemented for Sparse Matrix");
}

/**
 * From a matrix of any type, creates the three vectors of the triplet
 * (specific format for creating efficiently a Sparse matrix)
 * It only takes the only non-zero elements of the matrix
 * @param irows Output array of row indices
 * @param icols Output array of column indices
 * @param values Output array of non-zero values
 */
void MatrixSparse::getValuesAsTriplets(VectorInt &irows,
                                       VectorInt &icols,
                                       VectorDouble &values) const
{
  /// TODO : use cs_sparce corresponding function
  for (int icol = 0; icol < getNCols(); icol++)
    for (int irow = 0; irow < getNRows(); irow++)
    {
      if (!isValid(irow, icol)) continue;
      double value = getValue(irow, icol);
      if (ABS(value) < EPSILON10) continue;
      irows.push_back(irow);
      icols.push_back(icol);
      values.push_back(value);
    }
}

void MatrixSparse::_clear()
{
  _setNRows(0);
  _setNCols(0);
  _allocate();
}

Triplet MatrixSparse::getCsToTriplet(bool flag_from_1) const
{
  return csToTriplet(getCs(), flag_from_1);
}

int MatrixSparse::_getIndexToRank(int irow,int icol) const
{
  _forbiddenForSparse("_getIndexToRank");
  return ITEST;
}

MatrixSparse* toSparse(const AMatrix* matin)
{
  MatrixSparse* matout = new MatrixSparse(matin->getNRows(), matin->getNCols());

  // Create the triplet structure from the non-zero terms of source matrix

  VectorInt irows;
  VectorInt icols;
  VectorDouble values;
  matin->getValuesAsTriplets(irows, icols, values);

  // Load the triplet information in the cloned matrix

  matout->setValuesByArrays(irows, icols, values);

  return matout;
}
