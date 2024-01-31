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

#include <Eigen/SparseCholesky>

#include "csparse_d.h"
#include "csparse_f.h"

MatrixSparse::MatrixSparse(int nrow, int ncol, int opt_eigen)
    : AMatrix(nrow, ncol, opt_eigen),
      _csMatrix(nullptr),
      _eigenMatrix(),
      _flagDecomposeCholesky(false),
      _S(nullptr),
      _N(nullptr),
      _cholEigen()
{
  _allocate();
}

#ifndef SWIG
MatrixSparse::MatrixSparse(const cs *A, int opt_eigen)
    : AMatrix(cs_getnrow(A), cs_getncol(A), opt_eigen),
      _csMatrix(nullptr),
      _eigenMatrix(),
      _flagDecomposeCholesky(false), // Note: the class looses the Cholesky decomposition
      _S(nullptr),
      _N(nullptr),
      _cholEigen()
{
  if (_isFlagEigen())
    my_throw("Cannot copy a cs into an Eigen");
  else
    _csMatrix = cs_duplicate(A);
}
#endif

MatrixSparse::MatrixSparse(const MatrixSparse &m)
    : AMatrix(m),
      _csMatrix(nullptr),
      _eigenMatrix(),
      _flagDecomposeCholesky(false), // We loose the Cholesky decomposition
      _S(nullptr),
      _N(nullptr),
      _cholEigen()
{
  if (_isFlagEigen())
    _eigenMatrix = m._eigenMatrix;
  else
    _csMatrix = cs_duplicate(m._csMatrix);
}

MatrixSparse& MatrixSparse::operator=(const MatrixSparse &m)
{
  if (this != &m)
  {
    AMatrix::operator=(m);
    if (_isFlagEigen())
      _eigenMatrix = m._eigenMatrix;
    else
    {
      _csMatrix = cs_duplicate(m._csMatrix);

      // We loose the Cholesky decomposition when copying
      _S = nullptr;
      _N = nullptr;
    }
  }
  _flagDecomposeCholesky = false;
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

  if (_isFlagEigen())
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

int MatrixSparse::computeCholesky()
{
  if (! isSquare() || ! isSymmetric())
  {
    messerr("The 'Cholesky' decomposition is not possible as the matrix is not square symmetric");
    return 1;
  }

  if (_isFlagEigen())
  {
    _cholEigen.compute(_eigenMatrix);
  }
  else
  {
    _S = cs_schol(_csMatrix, 0);
    if (_S == nullptr)
    {
      messerr("Error in cs_schol function");
      return 1;
    }
    _N = cs_chol(_csMatrix, _S);
    if (_N == nullptr)
    {
      messerr("Error in cs_chol function");
      return 1;
    }
  }
  _flagDecomposeCholesky = true;
  return 0;
}

int MatrixSparse::solveCholesky(const VectorDouble& b, VectorDouble& x)
{
  if (! _flagDecomposeCholesky)
  {
    messerr("You must perform 'computeCholesky' beforehand");
    return 1;
  }
  int ncols = getNCols();
  if ((int) b.size() != ncols)
  {
    messerr("Dimension of input argument 'b' (%d) does not match",(int) b.size());
    messerr("the number of columns of the Matrix 'this' (%d)", ncols);
    return 1;
  }
  if ((int) x.size() != ncols)
  {
    messerr("Dimension of output argument 'x' (%d) does not match",(int) x.size());
    messerr("the number of columns of the Matrix 'this' (%d)", ncols);
    return 1;
  }

  if (_isFlagEigen())
  {
    Eigen::Map<const Eigen::VectorXd> bm(b.data(), getNCols());
    Eigen::Map<Eigen::VectorXd> xm(x.data(), getNRows());
    xm = _cholEigen.solve(bm);
  }
  else
  {
    VectorDouble work(ncols);
    cs_ipvec(ncols, _S->Pinv, b.data(), work.data());
    cs_lsolve(_N->L, work.data());
    cs_ltsolve(_N->L, work.data());
    cs_pvec(ncols, _S->Pinv, work.data(), x.data());
  }
  return 0;
}

void MatrixSparse::_transposeInPlace()
{
  if (_isFlagEigen())
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
  if (_isFlagEigen())
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
  if (_isFlagEigen())
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
  if (_isFlagEigen())
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

  if (_isFlagEigen())
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

  if (_isFlagEigen())
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
  if (_isFlagEigen())
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
  if (_isFlagEigen())
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
  if (_isFlagEigen())
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
  if (_isFlagEigen())
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
  if (_isFlagEigen())
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
  if (_isFlagEigen())
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
  if (_isFlagEigen())
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
  if (_isFlagEigen())
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
  if (_isFlagEigen())
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
  if (_isFlagEigen())
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
  if (_isFlagEigen())
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

MatrixSparse* MatrixSparse::createFromTriplet(const Triplet& T, int nrow, int ncol, int opt_eigen)
{
  MatrixSparse* mat = new MatrixSparse(nrow, ncol, opt_eigen);

  mat->setValuesFromTriplet(T);

  return mat;
}

void MatrixSparse::setValuesFromTriplet(const Triplet& T)
{
  if (_isFlagEigen())
  {
    _eigenMatrix.reserve(T.number);
    for (int k = 0; k < T.number; k++)
      _eigenMatrix.insert(T.rows[k], T.cols[k]) = T.values[k];
  }
  else
  {
    cs* Mtriplet = cs_spalloc(0,0,1,1,1);
    for (int i = 0; i < T.number; i++)
      (void) cs_entry(Mtriplet, T.rows[i], T.cols[i], T.values[i]);
    _csMatrix = cs_triplet(Mtriplet);
    Mtriplet = cs_spfree(Mtriplet);
  }
}

bool MatrixSparse::_isElementPresent(int irow, int icol) const
{
  if (_isFlagEigen())
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
  if (_isFlagEigen())
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

  if (_isFlagEigen())
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
  if (_isFlagEigen())
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
  if (_isFlagEigen())
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

  if (_isFlagEigen() && y._isFlagEigen())
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

  if (_isFlagEigen() && x._isFlagEigen() && y._isFlagEigen())
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

  if (_isFlagEigen() && y._isFlagEigen())
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
  if (_isFlagEigen())
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

  if (_isFlagEigen())
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

  if (_isFlagEigen())
  {
    sstr << AMatrix::toString(strfmt) << std::endl;
  }
  else
  {
    sstr << "- Number of rows    = " << getNRows() << std::endl;
    sstr << "- Number of columns = " << getNCols() << std::endl;
    sstr << "  (using Eigen Library)" << std::endl;
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
  if (_isFlagEigen())
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
  if (_isFlagEigen())
  {
    _eigenMatrix.data().squeeze();
  }
  else
  {
    _csMatrix = cs_spfree(_csMatrix);
    _flagDecomposeCholesky = false;
    _N = cs_nfree(_N);
    _S = cs_sfree(_S);
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
 */
Triplet MatrixSparse::getValuesAsTriplets() const
{
  Triplet T = triplet_init(0);
  if (_isFlagEigen())
  {
    T.nrows = _eigenMatrix.rows();
    T.ncols = _eigenMatrix.cols();
    int ecr = 0;
    for (int k = 0; k < _eigenMatrix.outerSize(); ++k)
    {
      for (Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix, k); it; ++it)
      {
        T.rows.push_back(it.row());
        T.cols.push_back(it.col());
        T.values.push_back(it.value());
        ecr++;
      }
    }
    T.number = ecr;
  }
  else
  {
    T = csToTriplet(_csMatrix, 0);
  }
  return T;
}

void MatrixSparse::_clear()
{
  _setNRows(0);
  _setNCols(0);
  _allocate();
}

Triplet MatrixSparse::getSparseToTriplet(bool flag_from_1) const
{
  if (_isFlagEigen())
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
    return csToTriplet(getCS(), flag_from_1);
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
  Triplet T = matin->getValuesAsTriplets();

  // Load the triplet information in the cloned matrix

  matout->setValuesFromTriplet(T);
  return matout;
}

void setUpdateNonZeroValue(int status)
{
  cs_set_status_update_nonzero_value(status);
}

int getUpdateNonZeroValue()
{
  return cs_get_status_update_nonzero_value();
}

const cs* _getCS(const MatrixSparse* A, bool optional)
{
  const cs* matCS = nullptr;
  if (A == nullptr)
  {
    if (optional)
      return matCS;
    else
      my_throw("Argument MatrixSparse is compulsory");
  }
  else
  {
    matCS = A->getCS();
    if (matCS == nullptr)
      my_throw("Argument should have a CS member");
  }
  return matCS;
}

cs* _getCSUnprotected(const MatrixSparse* A, bool optional)
{
  cs* matCS = nullptr;
  if (A == nullptr)
  {
    if (optional)
      return matCS;
    else
      my_throw("Argument MatrixSparse is compulsory");
  }
  else
  {
    matCS = A->getCSUnprotected();
    if (matCS == nullptr)
      my_throw("Argument should have a CS member");
  }
  return matCS;
}

MatrixSparse* matCS_glue(const MatrixSparse *A1,
                         const MatrixSparse *A2,
                         bool shiftRow,
                         bool shiftCol)
{
  const cs* A1cs = _getCS(A1);
  const cs* A2cs = _getCS(A2);
  cs* local = cs_glue(A1cs, A2cs, shiftRow, shiftCol);
  MatrixSparse* mat = new MatrixSparse(local, 0);
  local = cs_spfree(local);
  return mat;
}

MatrixSparse* matCS_matvecnorm(const MatrixSparse *A, const double *x, int oper)
{
  const cs* Acs = _getCS(A);
  cs* local = cs_matvecnorm(Acs, x, oper);
  MatrixSparse* mat = new MatrixSparse(local, 0);
  local = cs_spfree(local);
  return mat;
}

MatrixSparse* matCS_prod_norm(int mode, const MatrixSparse *A, const MatrixSparse *B)
{
  const cs* Acs = _getCS(A);
  const cs* Bcs = _getCS(B);
  cs* local = cs_prod_norm(mode, Acs, Bcs);
  MatrixSparse* mat = new MatrixSparse(local, 0);
  local = cs_spfree(local);
  return mat;
}

MatrixSparse* matCS_eye_tab(int number, double *values)
{
  cs* local = cs_eye_tab(number, values);
  MatrixSparse* mat = new MatrixSparse(local, 0);
  local = cs_spfree(local);
  return mat;
}

MatrixSparse* matCS_eye(int number, double value)
{
  cs* local = cs_eye(number, value);
  MatrixSparse* mat = new MatrixSparse(local, 0);
  local = cs_spfree(local);
  return mat;
}

MatrixSparse* matCS_triplet(const cs *T)
{
  cs* local = cs_triplet(T);
  MatrixSparse* mat = new MatrixSparse(local, 0);
  local = cs_spfree(local);
  return mat;
}

void matCS_tmulvec(const MatrixSparse *A, int nout, const double *x, double *y)
{
  const cs* Acs = _getCS(A);
  cs_tmulvec(Acs, nout, x, y);
}

void matCS_mulvec(const MatrixSparse *A, int nout, const double *x, double *y)
{
  const cs* Acs = _getCS(A);
  cs_mulvec(Acs, nout, x, y);
}

void matCS_vecmult(const MatrixSparse *A, int nout, const double *x, double *y)
{
  const cs* Acs = _getCS(A);
  cs_vecmult(Acs, nout, x, y);
}

MatrixSparse* matCS_prod_norm_diagonal(int mode, const MatrixSparse *B, VectorDouble diag)
{
  const cs* Bcs = _getCS(B);
  cs *local = cs_prod_norm_diagonal(mode, Bcs, diag);
  MatrixSparse *mat = new MatrixSparse(local, 0);
  local = cs_spfree(local);
  return mat;
}

MatrixSparse* matCS_transpose(const MatrixSparse *A, int values)
{
  const cs* Acs = _getCS(A);
  cs *local = cs_transpose(Acs, values);
  MatrixSparse *mat = new MatrixSparse(local, 0);
  local = cs_spfree(local);
  return mat;
}

MatrixSparse* matCS_multiply(const MatrixSparse *A, const MatrixSparse *B)
{
  const cs* Acs = _getCS(A);
  const cs* Bcs = _getCS(B);
  cs *local = cs_multiply(Acs, Bcs);
  MatrixSparse *mat = new MatrixSparse(local, 0);
  local = cs_spfree(local);
  return mat;
}

MatrixSparse* matCS_add(const MatrixSparse *A, const MatrixSparse *B, double alpha, double beta)
{
  const cs* Acs = _getCS(A);
  const cs* Bcs = _getCS(B, true);
  cs *local = cs_add(Acs, Bcs, alpha, beta);
  MatrixSparse *mat = new MatrixSparse(local, 0);
  local = cs_spfree(local);
  return mat;
}

int matCS_coarsening(MatrixSparse *Q, int type, int **indCo_ret, MatrixSparse **L_ret)
{
  const cs* Qcs = _getCS(Q);
  cs** LCS_ret = nullptr;
  int err = cs_coarsening(Qcs, type, indCo_ret, LCS_ret);
  if (err == 0) *L_ret = new MatrixSparse(*LCS_ret, 0);
  return err;
}

MatrixSparse* matCS_interpolate(MatrixSparse *AA, MatrixSparse *Lt, int *Co)
{
  const cs* AAcs = _getCS(AA);
  const cs* Ltcs = _getCS(Lt);
  cs *local = cs_interpolate(AAcs, Ltcs, Co);
  MatrixSparse *mat = new MatrixSparse(local, 0);
  local = cs_spfree(local);
  return mat;
}

MatrixSparse* matCS_extract_diag(MatrixSparse *C, int mode)
{
  const cs* Ccs = _getCS(C);
  cs* local = cs_extract_diag(Ccs, mode);
  MatrixSparse *mat = new MatrixSparse(local, 0);
  local = cs_spfree(local);
  return mat;
}

MatrixSparse* matCS_extract_submatrix_by_ranks(MatrixSparse *C, int *rank_rows, int *rank_cols)
{
  const cs* Ccs = _getCS(C);
  cs* local = cs_extract_submatrix_by_ranks(Ccs, rank_rows, rank_cols);
  MatrixSparse *mat = new MatrixSparse(local, 0);
  local = cs_spfree(local);
  return mat;
}

int matCS_gaxpy(const MatrixSparse *A, const double* x, double* y)
{
  const cs* Acs = _getCS(A);
  return cs_gaxpy(Acs, x, y);
}

void matCS_matvecnorm_inplace(MatrixSparse *A, const double* x, int oper)
{
  cs* Acs = _getCSUnprotected(A);
  cs_matvecnorm_inplace(Acs, x, oper);
  A->setCS(Acs);
}

double matCS_norm(const MatrixSparse *A)
{
  cs* Acs = _getCSUnprotected(A);
  return cs_norm(Acs);
}

MatrixSparse* matCS_prod_norm_single(int mode, MatrixSparse *B)
{
  const cs* Bcs = _getCS(B);
  cs* local = cs_prod_norm_single(mode, Bcs);
  MatrixSparse *mat = new MatrixSparse(local, 0);
  local = cs_spfree(local);
  return mat;
}

void matCS_add_value(const MatrixSparse *A, int row, int col, double value)
{
  const cs* Acs = _getCS(A);
  cs_add_value(Acs, row, col, value);
}

void matCS_set_cste(MatrixSparse *A, double value)
{
  const cs* Acs = _getCS(A);
  cs_set_cste(Acs, value);
}

VectorDouble  matCSD_extract_diag_VD(MatrixSparse *C, int mode)
{
  const cs* Ccs = _getCS(C);
  return csd_extract_diag_VD(Ccs, mode);
}

double* matCSD_extract_diag(const MatrixSparse *C, int mode)
{
  const cs* Ccs = _getCS(C);
  return csd_extract_diag(Ccs, mode);
}

MatrixSparse* matCS_diag(VectorDouble diag, double tol)
{
  cs* local = cs_diag(diag, tol);
  MatrixSparse *mat = new MatrixSparse(local, 0);
  local = cs_spfree(local);
  return mat;
}

int matCS_scale(MatrixSparse *A)
{
  const cs* Acs = _getCS(A);
  return cs_scale(Acs);
}

double matCS_get_value(const MatrixSparse *A, int row, int col)
{
  const cs* Acs = _getCS(A);
  return cs_get_value(Acs, row, col);
}

VectorInt matCS_color_coding(MatrixSparse *Q, int start, int *ncols)
{
  const cs* Qcs = _getCS(Q);
  return cs_color_coding(Qcs, start, ncols);
}

MatrixSparse* matCS_extract_submatrix_by_color(MatrixSparse *C,
                                               const VectorInt &colors,
                                               int ref_color,
                                               int row_ok,
                                               int col_ok)
{
  const cs* Ccs = _getCS(C);
  cs* local = cs_extract_submatrix_by_color(Ccs, colors, ref_color, row_ok, col_ok);
  MatrixSparse *mat = new MatrixSparse(local, 0);
  local = cs_spfree(local);
  return mat;
}

MatrixSparse* matCS_normalize_by_diag_and_release(MatrixSparse *Q, int flag_release)
{
  cs* Qcs = _getCSUnprotected(Q);
  cs *local = cs_normalize_by_diag_and_release(Qcs, flag_release);
  if (flag_release)
  {
    Q->setCS(local);
  }
  MatrixSparse *mat = new MatrixSparse(local, 0);
  local = cs_spfree(local);
  return mat;
}

MatrixSparse* matCS_add_and_release(MatrixSparse *b1,
                                    MatrixSparse *b2,
                                    double alpha,
                                    double beta,
                                    int flag_release)
{
  cs* b1cs = _getCSUnprotected(b1);
  const cs* b2cs = _getCS(b2, true);
  cs* local = cs_add_and_release(b1cs, b2cs, alpha, beta, flag_release);
  if (flag_release)
  {
    b1->setCS(local);
  }
  MatrixSparse *mat = new MatrixSparse(local, 0);
  local = cs_spfree(local);
  return mat;
}

MatrixSparse* matCS_multiply_and_release(MatrixSparse *b1, const MatrixSparse *b2,int flag_release)
{
  cs* b1cs = _getCSUnprotected(b1);
  const cs* b2cs = _getCS(b2, true);
  cs* local = cs_multiply_and_release(b1cs, b2cs, flag_release);
  if (flag_release)
  {
    b1->setCS(local);
  }
  MatrixSparse *mat = new MatrixSparse(local, 0);
  local = cs_spfree(local);
  return mat;
}
