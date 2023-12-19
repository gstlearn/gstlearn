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
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Matrix/LinkMatrixSparse.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"

#include <iostream>
#include <iomanip>

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
  if (isFlagEigen())
    my_throw("Cannot copy a cs into an Eigen");
  else
    _csMatrix = cs_duplicate(A);
}
#endif

MatrixSparse::MatrixSparse(const MatrixSparse &m)
    : AMatrix(m),
      _csMatrix(nullptr)
{
  if (isFlagEigen())
    _eigenMatrix = m._eigenMatrix;
  else
    _csMatrix = cs_duplicate(m._csMatrix);
}

MatrixSparse& MatrixSparse::operator=(const MatrixSparse &m)
{
  if (this != &m)
  {
    AMatrix::operator=(m);
    _deallocate();
    if (isFlagEigen())
      _eigenMatrix = m._eigenMatrix;
    else
      _csMatrix = cs_duplicate(m._csMatrix);
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
  _clearDecoration();
}

void MatrixSparse::reset(int nrows, int ncols, double value)
{
  DECLARE_UNUSED(value);
  if (! _isNumbersValid(nrows, ncols)) return;
  _forbiddenForSparse("reset");
}

void MatrixSparse::reset(int nrows, int ncols, const double* tab, bool byCol)
{
  DECLARE_UNUSED(byCol);
  DECLARE_UNUSED(tab);
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
  _clearDecoration();
}

void MatrixSparse::fillRandom(int seed, double zeroPercent)
{
  law_set_random_seed(seed);

  if (isFlagEigen())
  {
    for (int k=0; k<_eigenMatrix.outerSize(); ++k)
      for (Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix,k); it; ++it)
        it.valueRef() = law_gaussian();
  }
  else
  {
    cs *Atriplet;
    Atriplet = cs_spalloc(0, 0, 1, 1, 1);
    for (int irow = 0; irow < getNRows(); irow++)
      for (int icol = 0; icol < getNCols(); icol++)
      {
        if (!_isPhysicallyPresent(irow, icol)) continue;
        if (!mustBeDiagCst() && law_uniform(0., 1.) < zeroPercent) continue;
        cs_entry(Atriplet, irow, icol, law_gaussian());
      }
    _csMatrix = cs_triplet(Atriplet);
    Atriplet = cs_spfree(Atriplet);
  }
}

void MatrixSparse::_transposeInPlace()
{
  if (isFlagEigen())
  {
    Eigen::SparseMatrix<double> temp;
    temp = _eigenMatrix.transpose();
    _eigenMatrix = temp;
  }
  else
  {
    cs *old = _csMatrix;
    _csMatrix = cs_transpose(old, 0);
    cs_spfree(old);
  }
}

MatrixSparse* MatrixSparse::transpose() const
{
  MatrixSparse* mat = dynamic_cast<MatrixSparse*>(clone());
  if (isFlagEigen())
  {
    mat->_eigenMatrix = _eigenMatrix.transpose();
  }
  else
  {
    mat->transposeInPlace();
  }
  return mat;
}

void MatrixSparse::setColumn(int icol, const VectorDouble& tab)
{
  if (isFlagEigen())
  {
    for (int k=0; k < _eigenMatrix.outerSize(); ++k)
    {
      if (k != icol) continue;
      for (Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix,k); it; ++it)
      {
        int irow = it.row();
        it.valueRef() = tab[irow];
      }
    }
  }
  else
    AMatrix::setColumn(icol, tab);
}

void MatrixSparse::setRow(int irow, const VectorDouble& tab)
{
  if (isFlagEigen())
  {
    for (int k=0; k < _eigenMatrix.outerSize(); ++k)
    {
      for (Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix,k); it; ++it)
      {
        if (it.row() != irow) continue;
        int icol = it.col();
        it.valueRef() = tab[icol];
      }
    }
  }
  else
  {
    AMatrix::setRow(irow, tab);
  }
}

void MatrixSparse::setDiagonal(const VectorDouble& tab)
{
  if (! isSquare())
    my_throw("This function is only valid for Square matrices");

  if (isFlagEigen())
  {
    fill(0);
    for (int k=0; k < _eigenMatrix.outerSize(); ++k)
    {
      for (Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix,k); it; ++it)
      {
        int icol = it.col();
        if (it.row() != icol) continue;
        it.valueRef() = tab[icol];
      }
    }
  }
  else
  {
    cs *Mtriplet = cs_spalloc(0, 0, 1, 1, 1);
    for (int icol = 0; icol < getNCols(); icol++)
    {
      if (ABS(tab[icol]) < EPSILON10) continue;
      (void) cs_entry(Mtriplet, icol, icol, tab[icol]);
    }
    _csMatrix = cs_triplet(Mtriplet);
    Mtriplet = cs_spfree(Mtriplet);
  }
}

void MatrixSparse::setDiagonalToConstant(double value)
{
  if (! isSquare())
    my_throw("This function is only valid for Square matrices");

  if (isFlagEigen())
  {
    fill(0);
    for (int k=0; k < _eigenMatrix.outerSize(); ++k)
    {
      for (Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix,k); it; ++it)
      {
        int icol = it.col();
        if (it.row() != icol) continue;
        it.valueRef() = value;
      }
    }
  }
  else
  {
    cs *Mtriplet = cs_spalloc(0, 0, 1, 1, 1);
    for (int icol = 0; icol < getNCols(); icol++)
    {
      (void) cs_entry(Mtriplet, icol, icol, value);
    }
    _csMatrix = cs_triplet(Mtriplet);
    Mtriplet = cs_spfree(Mtriplet);
  }
}

/*! Gets the value for rank 'rank' */
double MatrixSparse::_getValueByRank(int rank) const
{
  DECLARE_UNUSED(rank);
  _forbiddenForSparse("_getValue (by rank)");
  return TEST;
}

double MatrixSparse::_getValue(int irow, int icol) const
{
  if (isFlagEigen())
  {
    return _eigenMatrix.coeff(irow, icol);
  }
  else
  {
    if (!_isIndexValid(irow, icol)) return TEST;
    return cs_get_value(_csMatrix, irow, icol);
  }
}

double& MatrixSparse::_getValueRef(int irow, int icol)
{
  DECLARE_UNUSED(irow);
  DECLARE_UNUSED(icol);
  _forbiddenForSparse("_getValueRef");
  return AMatrix::_getValueRef(irow, icol);
}

void MatrixSparse::_setValueByRank(int rank, double value)
{
  DECLARE_UNUSED(rank);
  DECLARE_UNUSED(value);
  _forbiddenForSparse("_setValue (by rank)");
}

void MatrixSparse::_setValue(int irow, int icol, double value)
{
  if (isFlagEigen())
  {
    _eigenMatrix.coeffRef(irow, icol) = value;
  }
  else
  {
    if (! _isIndexValid(irow, icol)) return;
    cs_set_value(_csMatrix, irow, icol, value);
  }
}

int MatrixSparse::_getMatrixPhysicalSize() const
{
  if (isFlagEigen())
    return _eigenMatrix.nonZeros();
  else
    return cs_nnz(_csMatrix);
}

/**
 * Change any nonzero term to 'value'
 * @param value Constant value used for filling 'this'
 */
void MatrixSparse::fill(double value)
{
  if (isFlagEigen())
  {
    for (int k=0; k < _eigenMatrix.outerSize(); ++k)
      for (Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix,k); it; ++it)
        it.valueRef() = value;
  }
  else
  {
    cs_set_cste(_csMatrix, value);
  }
}

/*! Multiply a Matrix row-wise */
void MatrixSparse::multiplyRow(const VectorDouble& vec)
{
  if (isFlagEigen())
  {
    for (int k=0; k < _eigenMatrix.outerSize(); ++k)
      for (Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix,k); it; ++it)
        it.valueRef() *= vec[it.row()];
  }
  else
  {
    cs* temp = cs_matvecR(_csMatrix, vec.data(), 0);
    cs_spfree(_csMatrix);
    _csMatrix = temp;
  }
}

/*! Multiply a Matrix column-wise */
void MatrixSparse::multiplyColumn(const VectorDouble& vec)
{
  if (isFlagEigen())
  {
    for (int k=0; k < _eigenMatrix.outerSize(); ++k)
      for (Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix,k); it; ++it)
        it.valueRef() *= vec[it.col()];
  }
  else
  {
    cs* temp = cs_matvecL(_csMatrix, vec.data(), 0);
    cs_spfree(_csMatrix);
    _csMatrix = temp;
  }
}

/*! Divide a Matrix row-wise */
void MatrixSparse::divideRow(const VectorDouble& vec)
{
  if (isFlagEigen())
  {
    for (int k=0; k < _eigenMatrix.outerSize(); ++k)
      for (Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix,k); it; ++it)
        it.valueRef() /= vec[it.row()];
  }
  else
  {
    cs* temp = cs_matvecR(_csMatrix, vec.data(), 1);
    cs_spfree(_csMatrix);
    _csMatrix = temp;
  }
}

/*! Divide a Matrix column-wise */
void MatrixSparse::divideColumn(const VectorDouble& vec)
{
  if (isFlagEigen())
  {
    for (int k=0; k < _eigenMatrix.outerSize(); ++k)
      for (Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix,k); it; ++it)
        it.valueRef() /= vec[it.col()];
  }
  else
  {
    cs* temp = cs_matvecL(_csMatrix, vec.data(), 1);
    cs_spfree(_csMatrix);
    _csMatrix = temp;
  }
}

/*! Perform M * 'vec' */
VectorDouble MatrixSparse::prodVector(const VectorDouble& vec) const
{
  if (isFlagEigen())
  {
    Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), getNCols());
    Eigen::VectorXd resm = _eigenMatrix * vecm;
    VectorDouble res(resm.data(), resm.data() + resm.size());
    return res;
  }
  else
  {
    int nrow = getNRows();
    VectorDouble res(nrow);
    cs_mulvec(_csMatrix, nrow, vec.data(), res.data());
    return res;
  }
}

/*! Perform 'vec'^T * M */
VectorDouble MatrixSparse::prodTVector(const VectorDouble& vec) const
{
  if (isFlagEigen())
  {
    Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), getNRows());
    Eigen::VectorXd resm = vecm.transpose() * _eigenMatrix;
    VectorDouble res(resm.data(), resm.data() + resm.size());
    return res;
  }
  else
  {
    int ncol = getNCols();
    VectorDouble res(ncol);
    cs_tmulvec(_csMatrix, ncol, vec.data(), res.data());
    return res;
  }
}

/**
 * Filling the matrix with an array of values
 * Note that this array is ALWAYS dimensioned to the total number
 * of elements in the matrix.
 * Kept for compatibility with old code where matrix contents was stored as
 * a double* array
 * @param values Input array (Dimension: nrow * ncol)
 * @param byCol true for Column Major; false for Row Major
 */
#ifndef SWIG
void MatrixSparse::_setValues(const double* values, bool byCol)
{
  int lec = 0;
  if (isFlagEigen())
  {
    if (byCol)
    {
      Eigen::Map<const Eigen::MatrixXd> temp(values, getNRows(), getNCols());
      _eigenMatrix = temp.sparseView(1., EPSILON10);
    }
    else
    {
      Eigen::Map<const Eigen::MatrixXd> temp(values, getNCols(), getNRows());
      _eigenMatrix = temp.transpose().sparseView(1., EPSILON10);
    }
  }
  else
  {
    cs *Mtriplet = cs_spalloc(0, 0, 1, 1, 1);
    if (byCol)
    {
      for (int icol = 0; icol < getNCols(); icol++)
        for (int irow = 0; irow < getNRows(); irow++, lec++)
        {
          if (ABS(values[lec]) < EPSILON10) continue;
          (void) cs_entry(Mtriplet, irow, icol, values[lec]);
        }
    }
    else
    {
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

  if (isFlagEigen())
  {
    _eigenMatrix.reserve(nelements);
    for (int k = 0; k < nelements; k++)
      _eigenMatrix.insert(irows[k], icols[k]) = values[k];
  }
  else
  {
    cs* Mtriplet = cs_spalloc(0,0,1,1,1);
    for (int i = 0; i < nelements; i++)
      (void) cs_entry(Mtriplet, irows[i], icols[i], values[i]);
    _csMatrix = cs_triplet(Mtriplet);
    Mtriplet = cs_spfree(Mtriplet);
  }
}

bool MatrixSparse::_isElementPresent(int irow, int icol) const
{
  if (isFlagEigen())
  {
    for (Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix, icol); it; ++it)
    {
      if (it.row() == irow) return true;
    }
    return false;
  }
  else
  {
    return cs_exist(_csMatrix, irow, icol);
  }
}

/**
 *
 * @param v Add a scalar value to all terms of the current matrix
 */
void MatrixSparse::addScalar(double v)
{
  if (v == 0.) return;
  if (isFlagEigen())
  {
    for (int k=0; k<_eigenMatrix.outerSize(); ++k)
      for (Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix,k); it; ++it)
        it.valueRef() += v;
  }
  else
  {
    for (int irow = 0; irow < getNRows(); irow++)
      for (int icol = 0; icol < getNCols(); icol++)
      {
        if (_isElementPresent(irow, icol))
          _setValue(irow, icol, _getValue(irow, icol) + v);
      }
  }
}

/**
 *
 * @param v Add constant value to the diagonal of the current Matrix
 */
void MatrixSparse::addScalarDiag(double v)
{
  if (v == 0.) return;

  if (isFlagEigen())
  {
    for (int k=0; k<_eigenMatrix.outerSize(); ++k)
      for (Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix,k); it; ++it)
      {
        if (it.col() == it.row())
          it.valueRef() += v;
      }
  }
  else
  {
    cs *csi = cs_eye(getNRows(), 1.);
    cs *res = cs_add(_csMatrix, csi, 1., v);
    cs_spfree(csi);
    cs_spfree(_csMatrix);
    _csMatrix = res;
  }
}

/**
 *
 * @param v Multiply all the terms of the matrix by the scalar 'v'
 */
void MatrixSparse::prodScalar(double v)
{
  if (v == 1.) return;
  if (isFlagEigen())
  {
    for (int k=0; k<_eigenMatrix.outerSize(); ++k)
      for (Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix,k); it; ++it)
        it.valueRef() *= v;
  }
  else
  {
    cs *res = cs_add(_csMatrix, _csMatrix, v, 0.);
    cs_spfree(_csMatrix);
    _csMatrix = res;
  }
}

/**
 *
 * @param inv Input vector
 * @param outv Output vector obtained by multiplying 'inv' by current Matrix
 */
#ifndef SWIG
void MatrixSparse::_prodVectorInPlace(const double *inv, double *outv) const
{
  if (isFlagEigen())
  {
    Eigen::Map<const Eigen::VectorXd> inm(inv, getNCols());
    Eigen::Map<Eigen::VectorXd> outm(outv, getNRows());
    outm.noalias() = _eigenMatrix * inm;
  }
  else
  {
    cs_vecmult(_csMatrix, getNRows(), inv, outv);
  }
}
#endif

/**
 * Add the matrix 'y' to the current Matrix
 * @param y Matrix to be added
 * @param value Multiplicative coefficient
 */
void MatrixSparse::addMatrix(const MatrixSparse& y, double value)
{
  if (! isSameSize(y))
  {
    messerr("Matrices 'y' and 'this' should have the same size");
    return;
  }

  if (isFlagEigen())
  {
    _eigenMatrix += y._eigenMatrix * value;
  }
  else
  {
    cs *res = cs_add(_csMatrix, y._csMatrix, 1., value);
    cs_spfree(_csMatrix);
    _csMatrix = res;
  }
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
      messerr("- First matrix:  NRows = %d - NColumns = %d", x.getNRows(), x.getNCols());
      messerr("- Second matrix: NRows = %d - NColumns = %d", y.getNRows(), y.getNCols());
      messerr("- Result matrix: NRows = %d - NColumns = %d", getNRows(), getNCols());
      messerr("Operation is cancelled");
      return;
    }
  }

  if (isFlagEigen())
  {
    _eigenMatrix = x._eigenMatrix * y._eigenMatrix;
  }
  else
  {
    cs* res = cs_multiply(x._csMatrix, y._csMatrix);
    cs_spfree(_csMatrix);
    _csMatrix = res;
  }
}

/**
 * Store the product of 'transpose(x)' by 'y' in this
 * @param x First Matrix
 * @param y Second matrix
 */
void MatrixSparse::prodTMatrix(const MatrixSparse& x, const MatrixSparse& y)
{
  if (_getFlagCheckAddress())
  {
    if (x.getNRows() != y.getNRows() ||
        x.getNCols() != getNRows()   ||
        y.getNCols() != getNCols())
    {
      messerr("Incompatible matrix dimensions for matrix product");
      messerr("- First matrix:  NRows = %d - NColumns = %d", x.getNRows(), x.getNCols());
      messerr("- Second matrix: NRows = %d - NColumns = %d", y.getNRows(), y.getNCols());
      messerr("- Result matrix: NRows = %d - NColumns = %d", getNRows(), getNCols());
      messerr("Operation is cancelled");
      return;
    }
  }

  if (isFlagEigen())
  {
    _eigenMatrix = x._eigenMatrix.transpose() * y._eigenMatrix;
  }
  else
  {
    cs* xT = cs_transpose(x._csMatrix, 1);
    cs* res = cs_multiply(xT, y._csMatrix);
    cs_spfree(_csMatrix);
    _csMatrix = res;
    delete xT;
  }
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
    my_throw("Matrices should have same size");
  if (!y.isSparse())
    my_throw("This function can only combine sparse matrices together");

  if (isFlagEigen())
  {
    _eigenMatrix = cx * _eigenMatrix + cy * y._eigenMatrix;
  }
  else
  {
    cs *res = cs_add(_csMatrix, y._csMatrix, cx, cy);
    cs_spfree(_csMatrix);
    _csMatrix = res;
  }
}

int MatrixSparse::_invert()
{
  if (!isSquare())
    my_throw("Invert method is restricted to Square matrices");
  if (isFlagEigen())
  {
    int n = getNCols();
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver;
    solver.compute(_eigenMatrix);
    Eigen::SparseMatrix<double> I(n,n);
    I.setIdentity();
    _eigenMatrix = solver.solve(I);
  }
  else
  {
    cs *inv = cs_invert(_csMatrix, 0);
    _deallocate();
    _csMatrix = inv;
  }
  return 0;
}

int MatrixSparse::_solve(const VectorDouble& b, VectorDouble& x) const
{
  int error = 0;

  if (! isSymmetric())
    my_throw("Invert method is limited to Square Symmetrical Matrices");
  if ((int) b.size() != getNRows() || (int) x.size() != getNRows())
    my_throw("b' and 'x' should have the same dimension as the Matrix");

  if (isFlagEigen())
  {
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver;
    Eigen::Map<const Eigen::VectorXd> bm(b.data(), getNCols());
    Eigen::Map<Eigen::VectorXd> xm(x.data(), getNRows());
    xm = solver.compute(_eigenMatrix).solve(bm);
  }
  else
  {
    x = b;
    error = cs_cholsol(_csMatrix,x.data(), 0);
  }
  return error;
}

String MatrixSparse::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  if (isFlagEigen())
  {
    sstr << AMatrix::toString(strfmt) << std::endl;
  }
  else
  {
    sstr << "- Number of rows    = " << getNRows() << std::endl;
    sstr << "- Number of columns = " << getNCols() << std::endl;
    sstr << "- Sparse Format" << std::endl;
    sstr << toMatrix(String(), _csMatrix);
  }
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
  if (isFlagEigen())
  {
    if (isMultiThread()) omp_set_num_threads(getMultiThread());
    _eigenMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor>(getNRows(),getNCols());
  }
  else
  {
    cs *Atriplet;
    int nrow = getNRows();
    int ncol = getNCols();
    if (nrow > 0 && ncol > 0)
    {
      Atriplet = cs_spalloc(0, 0, 1, 1, 1);
      cs_entry(Atriplet, getNRows() - 1, getNCols() - 1, 0.);
      _csMatrix = cs_triplet(Atriplet);
      Atriplet = cs_spfree(Atriplet);
    }
  }
}

void MatrixSparse::_deallocate()
{
  if (isFlagEigen())
  {
    // This is where the specific code should take place
  }
  else
  {
    _csMatrix = cs_spfree(_csMatrix);
  }
}

void MatrixSparse::_forbiddenForSparse(const String& func) const
{
  messerr("Problem with Function: %s",func.c_str());
  messerr("This function is not available in Sparse Matrix");
  return;
}

void MatrixSparse::dumpElements(const String& title, int ifrom, int ito) const
{
  DECLARE_UNUSED(title);
  DECLARE_UNUSED(ifrom);
  DECLARE_UNUSED(ito);
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
  if (isFlagEigen())
  {
    for (int k = 0; k < _eigenMatrix.outerSize(); ++k)
    {
      for (Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix, k); it; ++it)
      {
        irows.push_back(it.row());
        icols.push_back(it.col());
        values.push_back(it.value());
      }
    }
  }
  else
  {
    int number = 0;
    int *cols = nullptr;
    int *rows = nullptr;
    double *vals = nullptr;
    cs_sparse_to_triplet(_csMatrix, 0, &number, &cols, &rows, &vals);
    irows = VH::initVInt(rows, number);
    icols = VH::initVInt(cols, number);
    values = VH::initVDouble(vals, number);
  }
}

void MatrixSparse::_clear()
{
  _setNRows(0);
  _setNCols(0);
  _allocate();
}

Triplet MatrixSparse::getSparseToTriplet(bool flag_from_1) const
{
  if (isFlagEigen())
  {
    Triplet trp;
    trp.flagFromOne = flag_from_1;
    trp.nrows = getNRows();
    trp.ncols = getNCols();
    trp.rows = VectorInt();
    trp.cols = VectorInt();
    trp.values = VectorDouble();
    for (int k = 0; k < _eigenMatrix.outerSize(); ++k)
    {
      for (Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix, k); it; ++it)
      {
        if (flag_from_1)
        {
          trp.rows.push_back(it.row() + 1);
          trp.cols.push_back(it.col() + 1);
        }
        else
        {
          trp.rows.push_back(it.row());
          trp.cols.push_back(it.col());
        }
        trp.values.push_back(it.value());
      }
    }
    return trp;
  }
  else
  {
    return csToTriplet(getCs(), flag_from_1);
  }
}

int MatrixSparse::_getIndexToRank(int irow,int icol) const
{
  DECLARE_UNUSED(irow);
  DECLARE_UNUSED(icol);
  _forbiddenForSparse("_getIndexToRank");
  return ITEST;
}

MatrixSparse* createFromAnyMatrix(const AMatrix* matin)
{
  // Create the output Sparse Matrix

  MatrixSparse* matout = new MatrixSparse(matin->getNRows(),matin->getNCols());

  // Create the triplet structure from the non-zero terms of source matrix

  VectorInt irows;
  VectorInt icols;
  VectorDouble values;
  matin->getValuesAsTriplets(irows, icols, values);

  // Load the triplet information in the cloned matrix

  matout->setValuesByArrays(irows, icols, values);
  return matout;
}

GSTLEARN_EXPORT void setUpdateNonZeroValue(int status)
{
  cs_set_status_update_nonzero_value(status);
}

GSTLEARN_EXPORT int getUpdateNonZeroValue()
{
  return cs_get_status_update_nonzero_value();
}
