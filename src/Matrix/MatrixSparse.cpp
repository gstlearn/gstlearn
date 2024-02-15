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
#include "LinearOp/Cholesky.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"

#include <iostream>
#include <iomanip>

#include <Eigen/SparseCholesky>

#include <csparse_f.h>

MatrixSparse::MatrixSparse(int nrow, int ncol, int opt_eigen)
    : AMatrix(nrow, ncol, opt_eigen),
      _csMatrix(nullptr),
      _eigenMatrix(),
      _factor(nullptr)
{
  _allocate();
}

#ifndef SWIG
MatrixSparse::MatrixSparse(const cs *A, int opt_eigen)
    : AMatrix(cs_get_nrow(A), cs_get_ncol(A), opt_eigen),
      _csMatrix(nullptr),
      _eigenMatrix(),
      _factor(nullptr)
{
  if (isFlagEigen())
    my_throw("Cannot copy a cs into an Eigen"); // TODO: improve convert cs to MatrixSparse
  else
    _csMatrix = cs_duplicate(A);
}
#endif

MatrixSparse::MatrixSparse(const MatrixSparse &m)
    : AMatrix(m),
      _csMatrix(nullptr),
      _eigenMatrix(),
      _factor(nullptr) // recompute cholesky (if needed)
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
    if (isFlagEigen())
      _eigenMatrix = m._eigenMatrix;
    else
    {
      _csMatrix = cs_duplicate(m._csMatrix);
    }
    _factor = nullptr; // TODO recompute Cholesky (if needed)
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
    cs *local = cs_spalloc2(0, 0, 1, 1, 1);
    for (int irow = 0; irow < getNRows(); irow++)
      for (int icol = 0; icol < getNCols(); icol++)
      {
        if (!_isPhysicallyPresent(irow, icol)) continue;
        if (!mustBeDiagCst() && law_uniform(0., 1.) < zeroPercent) continue;
        cs_entry2(local, irow, icol, law_gaussian());
      }
    _csMatrix = cs_triplet2(local);
    local = cs_spfree2(local);
  }
}

int MatrixSparse::computeCholesky()
{
  if (! isSquare())
  {
    messerr("The 'Cholesky' decomposition is not possible as the matrix is not square");
    return 1;
  }
  if (_factor != nullptr) return 0;
  _factor = new Cholesky(this);
  return _factor == nullptr;
}

int MatrixSparse::solveCholesky(const VectorDouble& b, VectorDouble& x)
{
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
  if (_factor == nullptr)
  {
    messerr("Use 'computeCholesky' beforehand");
    return 1;
  }
  return _factor->solve(b, x);
}

int MatrixSparse::simulateCholesky(const VectorDouble &b, VectorDouble &x)
{
  int ncols = getNCols();
  if ((int) b.size() != ncols)
  {
    messerr("Dimension of input argument 'b' (%d) does not match", (int) b.size());
    messerr("the number of columns of the Matrix 'this' (%d)", ncols);
    return 1;
  }
  if ((int) x.size() != ncols)
  {
    messerr("Dimension of output argument 'x' (%d) does not match", (int) x.size());
    messerr("the number of columns of the Matrix 'this' (%d)", ncols);
    return 1;
  }
  if (_factor == nullptr)
  {
    messerr("Use 'computeCholesky' beforehand");
    return 1;
  }

  return _factor->simulate(b, x);
}

double MatrixSparse::getCholeskyLogDeterminant()
{
  if (_factor == nullptr)
  {
    messerr("Use 'computeCholesky' beforehand");
    return TEST;
  }
  return _factor->getLogDeterminant();
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
    cs_spfree2(old);
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

/**
 * Fill a column of an already existing Sparse matrix, using 'tab' as entry
 * The input 'tab' corresponds to the whole column contents
 * @param icol Column rank
 * @param tab  Vector containing the information (Dimension: nrows)
 *
 * @warning: This method only copies the values at the non-zero existing entries
 */
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

/**
 * Fill a row of an already existing Sparse matrix, using 'tab' as entry
 * The input 'tab' corresponds to the whole row contents
 * @param irow Row rank
 * @param tab  Vector containing the information (Dimension: ncols)
 *
 * @warning: This method only copies the values at the non-zero existing entries
 */
void MatrixSparse::setRow(int irow, const VectorDouble& tab)
{
  if (isFlagEigen())
  {
    for (int k=0; k < _eigenMatrix.outerSize(); ++k)
      for (Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix,k); it; ++it)
      {
        if (it.row() != irow) continue;
        int icol = it.col();
        it.valueRef() = tab[icol];
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
    Eigen::Map<const Eigen::VectorXd> vecm(tab.data(), tab.size());
    _eigenMatrix = vecm.asDiagonal();
  }
  else
  {
    cs *local = cs_spalloc2(0, 0, 1, 1, 1);
    for (int icol = 0, ncol = getNCols(); icol < ncol; icol++)
    {
      if (ABS(tab[icol]) < EPSILON10) continue;
      (void) cs_entry2(local, icol, icol, tab[icol]);
    }
    _csMatrix = cs_triplet2(local);
    local = cs_spfree2(local);
  }
}

void MatrixSparse::setDiagonalToConstant(double value)
{
  if (! isSquare())
    my_throw("This function is only valid for Square matrices");

  if (isFlagEigen())
  {
    VectorDouble vec(getNRows(), value);
    Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), vec.size());
    _eigenMatrix = vecm.asDiagonal();
  }
  else
  {
    cs *local = cs_spalloc2(0, 0, 1, 1, 1);
    for (int icol = 0, ncol = getNCols(); icol < ncol; icol++)
    {
      (void) cs_entry2(local, icol, icol, value);
    }
    _csMatrix = cs_triplet2(local);
    local = cs_spfree2(local);
  }
}

/*! Gets the value for rank 'rank' */
double MatrixSparse::_getValueByRank(int rank) const
{
  DECLARE_UNUSED(rank);
  _forbiddenForSparse("_getValueByRank");
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
    cs* temp = cs_matvecR(_csMatrix, vec.data(), 1);
    _csMatrix = cs_spfree2(_csMatrix);
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
    cs* temp = cs_matvecL(_csMatrix, vec.data(), 1);
    _csMatrix = cs_spfree2(_csMatrix);
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
    cs* temp = cs_matvecR(_csMatrix, vec.data(), -1);
    _csMatrix = cs_spfree2(_csMatrix);
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
    cs* temp = cs_matvecL(_csMatrix, vec.data(), -1);
    _csMatrix = cs_spfree2(_csMatrix);
    _csMatrix = temp;
  }
}

/*! Perform y = x %*% 'this' */
VectorDouble MatrixSparse::prodVecMat(const VectorDouble& x, bool transpose) const
{
  if (isFlagEigen())
  {
    Eigen::Map<const Eigen::VectorXd> xm(x.data(), x.size());
    Eigen::VectorXd ym;
    if (transpose)
      ym = xm.transpose() * _eigenMatrix.transpose();
    else
      ym = xm.transpose() * _eigenMatrix;
    VectorDouble y(ym.data(), ym.data() + ym.size());
    return y;
  }
  else
  {
    VectorDouble y;
    if (transpose)
    {
      int ncol = getNCols();
      y.resize(ncol);
      cs_vector_xtM(_csMatrix, ncol, x.data(), y.data());
    }
    else
    {
      int nrow = getNRows();
      y.resize(nrow);
      cs_vector_xM(_csMatrix, nrow, x.data(), y.data());
    }
    return y;
  }
}

/*! Perform y = 'this' %*% x */
VectorDouble MatrixSparse::prodMatVec(const VectorDouble& x, bool transpose) const
{
  if (isFlagEigen())
  {
    Eigen::Map<const Eigen::VectorXd> xm(x.data(), x.size());
    Eigen::VectorXd ym;
    if (transpose)
      ym = _eigenMatrix.transpose() * xm;
    else
      ym = _eigenMatrix * xm;
    VectorDouble y(ym.data(), ym.data() + ym.size());
    return y;
  }
  else
  {
    VectorDouble y;
    if (transpose)
    {
      int ncol = getNCols();
      y.resize(ncol);
      cs_vector_tMx(_csMatrix, ncol, x.data(), y.data());
    }
    else
    {
      int nrow = getNRows();
      y.resize(nrow);
      cs_vector_Mx(_csMatrix, nrow, x.data(), y.data());
    }
    return y;
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
    cs *local = cs_spalloc2(0, 0, 1, 1, 1);
    if (byCol)
    {
      for (int icol = 0; icol < getNCols(); icol++)
        for (int irow = 0; irow < getNRows(); irow++, lec++)
        {
          if (ABS(values[lec]) < EPSILON10) continue;
          (void) cs_entry2(local, irow, icol, values[lec]);
        }
    }
    else
    {
      for (int irow = 0; irow < getNRows(); irow++)
        for (int icol = 0; icol < getNCols(); icol++, lec++)
        {
          if (ABS(values[lec]) < EPSILON10) continue;
          (void) cs_entry2(local, irow, icol, values[lec]);
        }
    }
    _csMatrix = cs_triplet2(local);
    local = cs_spfree2(local);
  }
}
#endif

MatrixSparse* MatrixSparse::createFromTriplet(const NF_Triplet &NF_T,
                                              int nrow,
                                              int ncol,
                                              int opt_eigen)
{
  // If 'nrow' and 'ncol' are not defined, derive them from NF_T
  if (nrow <= 0 || ncol <= 0)
  {
    nrow = VH::maximum(NF_T.rows) + 1;
    ncol = VH::maximum(NF_T.cols) + 1;
  }
  MatrixSparse* mat = new MatrixSparse(nrow, ncol, opt_eigen);

  mat->setValuesFromTriplet(NF_T);

  return mat;
}

MatrixSparse* MatrixSparse::addMatMat(const MatrixSparse *x, const MatrixSparse *y, double cx, double cy)
{
  MatrixSparse* mat = new MatrixSparse(x->getNRows(), x->getNCols(), x->isFlagEigen());
  if (x->isFlagEigen() && y->isFlagEigen())
  {
    mat->_eigenMatrix = cx * x->_eigenMatrix + cy * y->_eigenMatrix;
  }
  else
  {
    mat->_csMatrix = cs_spfree2(mat->_csMatrix);
    mat->_csMatrix = cs_add(x->_csMatrix, y->_csMatrix, cx, cy);
  }
  return mat;
}

MatrixSparse* MatrixSparse::diagVec(const VectorDouble& vec, int opt_eigen)
{
  int size = (int) vec.size();
  MatrixSparse *mat = new MatrixSparse(size, size, opt_eigen);

  if (mat->isFlagEigen())
  {
    mat->setDiagonal(vec);
  }
  else
  {
    mat->_csMatrix = cs_spfree2(mat->_csMatrix);
    mat->_csMatrix = cs_diag(vec, EPSILON10);
  }
  return mat;
}

MatrixSparse* MatrixSparse::diagConstant(int number, double value, int opt_eigen)
{
  MatrixSparse *mat = new MatrixSparse(number, number, opt_eigen);

  if (mat->isFlagEigen())
  {
    mat->setDiagonalToConstant(value);
  }
  else
  {
    mat->_csMatrix  = cs_spfree2(mat->_csMatrix);
    VectorDouble vec = VH::initVDouble(number, value);
    mat->_csMatrix = cs_diag(vec, EPSILON10);
  }
  return mat;
}

/**
 * Construct a sparse matrix with the diagonal of 'A', where each element is transformed
 * @param A    Input sparse matrix
 * @param oper_choice: Operation on the diagonal term (see Utilities::operate_XXX)
 * @param opt_eigen Option for choosing Eigen Library or not
 * @return
 */
MatrixSparse* MatrixSparse::diagMat(MatrixSparse *A, int oper_choice, int opt_eigen)
{
  if (! A->isSquare())
  {
    messerr("This method requires the matrix 'A' to be square");
    return nullptr;
  }

  VectorDouble diag = A->getDiagonal();
  VectorHelper::transformVD(diag, oper_choice);
  return MatrixSparse::diagVec(diag, opt_eigen);
}

void MatrixSparse::setValuesFromTriplet(const NF_Triplet& NF_T)
{
  if (isFlagEigen())
  {
    _eigenMatrix.reserve(NF_T.number);
    for (int k = 0; k < NF_T.number; k++)
      _eigenMatrix.insert(NF_T.rows[k], NF_T.cols[k]) = NF_T.values[k];
  }
  else
  {
    cs* local = cs_spalloc2(0,0,1,1,1);
    for (int i = 0; i < NF_T.number; i++)
      (void) cs_entry2(local, NF_T.rows[i], NF_T.cols[i], NF_T.values[i]);
    _csMatrix = cs_triplet2(local);
    local = cs_spfree2(local);
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

void MatrixSparse::addValue(int row, int col, double value)
{
  if (ABS(value) <= EPSILON10) return;
  if (isFlagEigen())
    _eigenMatrix.coeffRef(row, col) += value;
  else
    cs_add_value(_csMatrix, row, col, value);
}

double MatrixSparse::getValue(int row, int col) const
{
  if (isFlagEigen())
    return _eigenMatrix.coeff(row, col);
  else
    return cs_get_value(_csMatrix, row, col);
}

double MatrixSparse::L1Norm() const
{
  if (isFlagEigen())
    return (Eigen::RowVectorXd::Ones(_eigenMatrix.rows()) * _eigenMatrix.cwiseAbs()).maxCoeff();
  else
    return cs_norm(_csMatrix);
}

void MatrixSparse::getStats(int *nrows, int *ncols, int *count, double *percent) const
{
  if (isFlagEigen())
  {
    *nrows = getNRows();
    *ncols = getNCols();
    *count = _eigenMatrix.nonZeros();
    *percent = 0.;
    if ((*nrows) > 0 && (*ncols) > 0)
      (*percent) = ((100. * (double) (*count))
          / ((double) (*nrows) * (double) (*ncols)));
  }
  else
  {
    cs_rowcol(_csMatrix, nrows, ncols, count, percent);
  }
}

VectorDouble MatrixSparse::extractDiag(int oper_choice) const
{
  if (isFlagEigen())
  {
    Eigen::VectorXd ym = _eigenMatrix.diagonal();
    VectorDouble diag(ym.data(), ym.data() + ym.size());
    VH::transformVD(diag, oper_choice);
    return diag;
  }
  else
  {
    VectorDouble diag = csd_extract_diag_VD(_csMatrix, 1);
    VH::transformVD(diag, oper_choice);
    return diag;
  }
}

int MatrixSparse::addVecInPlace(const VectorDouble& x, VectorDouble& y)
{
  if (isFlagEigen())
  {
    Eigen::Map<const Eigen::VectorXd> xm(x.data(), x.size());
    Eigen::Map<Eigen::VectorXd> ym(y.data(), y.size());
    ym = _eigenMatrix * xm + ym;
    return 0;
  }
  else
  {
    return (!cs_gaxpy(_csMatrix, x.data(), y.data()));
  }
}

void MatrixSparse::setConstant(double value)
{
  if (isFlagEigen())
  {
    for (int k=0; k<_eigenMatrix.outerSize(); ++k)
      for (Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix,k); it; ++it)
        it.valueRef() = value;
  }
  else
  {
    cs_set_cste(_csMatrix, value);
  }
}

int MatrixSparse::scaleByDiag()
{
  if (isFlagEigen())
  {
    VectorDouble diag = extractDiag(-1);
    Eigen::Map<const Eigen::VectorXd> ym(diag.data(), diag.size());
    _eigenMatrix = ym.asDiagonal() * _eigenMatrix;
    return 0;
  }
  else
  {
    return cs_scale(_csMatrix);
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
    cs_spfree2(csi);
    cs_spfree2(_csMatrix);
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
    cs_spfree2(_csMatrix);
    _csMatrix = res;
  }
}

/**
 * Returns 'y' = 'this' %*% 'x'
 * @param x Input vector
 * @param y Output vector
 * @param transpose True if the matrix 'this' must be transposed
 */
void MatrixSparse::_prodMatVecInPlacePtr(const double *x, double *y, bool transpose) const
{
  if (isFlagEigen())
  {
    if (transpose)
    {
      Eigen::Map<const Eigen::VectorXd> xm(x, getNRows());
      Eigen::Map<Eigen::VectorXd> ym(y, getNCols());
      ym = _eigenMatrix.transpose() * xm;
    }
    else
    {
      Eigen::Map<const Eigen::VectorXd> xm(x, getNCols());
      Eigen::Map<Eigen::VectorXd> ym(y, getNRows());
      ym = _eigenMatrix * xm;
    }
  }
  else
  {
    if (transpose)
      cs_vector_tMx(_csMatrix, getNCols(), x, y);
    else
      cs_vector_Mx(_csMatrix, getNRows(), x, y);
  }
}

/**
 * Returns 'y' = 'x' %*% 'this'
 * @param x Input vector
 * @param y Output vector
 * @param transpose True if the matrix 'this' must be transposed
 */
void MatrixSparse::_prodVecMatInPlacePtr(const double *x, double *y, bool transpose) const
{
  if (isFlagEigen())
  {
    if (transpose)
    {
      Eigen::Map<const Eigen::VectorXd> xm(x, getNCols());
      Eigen::Map<Eigen::VectorXd> ym(y, getNRows());
      ym = xm.transpose() * _eigenMatrix.transpose() * xm;
    }
    else
    {
      Eigen::Map<const Eigen::VectorXd> xm(x, getNRows());
      Eigen::Map<Eigen::VectorXd> ym(y, getNCols());
      ym = xm.transpose() * _eigenMatrix;
    }
  }
  else
  {
    if (transpose)
      cs_vector_xtM(_csMatrix, getNRows(), x, y);
    else
      cs_vector_xM(_csMatrix, getNCols(), x, y);
  }
}

/**
 * Store the product of 'x' by 'y' in this
 * @param x First Matrix
 * @param y Second matrix
 * @param transposeX True if First matrix is transposed
 * @param transposeY True if Second matrix is transposed
 */
void MatrixSparse::prodMatMatInPlace(const AMatrix *x,
                                     const AMatrix *y,
                                     bool transposeX,
                                     bool transposeY)
{
  if (!_checkLink(x->getNRows(), x->getNCols(), transposeX,
                  y->getNRows(), y->getNCols(), transposeY)) return;

  const MatrixSparse* xm = dynamic_cast<const MatrixSparse*>(x);
  const MatrixSparse* ym = dynamic_cast<const MatrixSparse*>(y);
  if (xm == nullptr || ym == nullptr)
  {
    AMatrix::prodMatMatInPlace(x, y, transposeX, transposeY);
  }
  else
  {
    if (isFlagEigen() && xm->isFlagEigen() && ym->isFlagEigen())
    {
      if (transposeX)
      {
        if (transposeY)
        {
          _eigenMatrix = xm->_eigenMatrix.transpose() * ym->_eigenMatrix.transpose();
        }
        else
        {
          _eigenMatrix = xm->_eigenMatrix.transpose() * ym->_eigenMatrix;
        }
      }
      else
      {
        if (transposeY)
        {
          _eigenMatrix = xm->_eigenMatrix * ym->_eigenMatrix.transpose();
        }
        else
        {
          _eigenMatrix = xm->_eigenMatrix * ym->_eigenMatrix;
        }
      }
    }
    else
    {
      cs* csx = xm->_csMatrix;
      if (transposeX) csx = cs_transpose(csx, 1);
      cs *csy = ym->_csMatrix;
      if (transposeY) csy = cs_transpose(csy, 1);
      cs *res = cs_multiply(csx, csy);
      if (transposeX) csx = cs_spfree2(csx);
      if (transposeY) csy = cs_spfree2(csy);
      cs_spfree2(_csMatrix);
      _csMatrix = res;
    }
  }
}

MatrixSparse* prodNormMatMat(const MatrixSparse &a,
                             const MatrixSparse &m,
                             bool transpose)
{
  int nrow = (transpose) ? a.getNCols() : a.getNRows();
  int ncol = (transpose) ? a.getNRows() : a.getNCols();
  MatrixSparse *mat = new MatrixSparse(nrow, ncol, a.isFlagEigen());
  mat->prodNormMatMatInPlace(a, m, transpose);
  return mat;
}

MatrixSparse* prodNormMat(const MatrixSparse &a, const VectorDouble& vec, bool transpose)
{
  int nsym = (transpose) ? a.getNCols() : a.getNRows();
  MatrixSparse *mat = new MatrixSparse(nsym, nsym, a.isFlagEigen());
  mat->prodNormMatInPlace(a, vec, transpose);
  return mat;
}

MatrixSparse* prodNormDiagVec(const MatrixSparse &a,
                              const VectorDouble &vec,
                              int oper_choice)
{
  int nrow = a.getNRows();
  int ncol = a.getNCols();
  MatrixSparse *mat = new MatrixSparse(nrow, ncol, a.isFlagEigen());

  if (a.isFlagEigen())
  {
    // Perform the transformation of the input vector
    VectorDouble vecp = vec;
    VH::transformVD(vecp, oper_choice);

    Eigen::Map<const Eigen::VectorXd> vecm(vecp.data(), vecp.size());
    auto diag = vecm.asDiagonal();
    mat->setEigenMatrix(diag * a.getEigenMatrix() * diag);
  }
  else
  {
    cs* local = cs_matvecnorm(a.getCS(), vec.data(), oper_choice);
    mat->setCS(local);
    local = cs_spfree2(local);
  }
  return mat;
}

/**
 * Perform: 'this' = diag('vec') %*% 'A' %*% diag('vec')
 * @param vec  Input Vector
 * @param oper_choice Type of transformation
 */
void MatrixSparse::prodNormDiagVecInPlace(const VectorDouble &vec, int oper_choice)
{
  if (! isSquare())
  {
    messerr("This method is limited to square matrices");
    return;
  }
  if (getNRows() != (int) vec.size())
  {
    messerr("Matrix dimension (%d) does not math vector dimension (%d)",
            getNRows(), (int) vec.size());
    return;
  }

  if (isFlagEigen())
  {
    // Perform the transformation of the input vector
    VectorDouble vecp = vec;
    VH::transformVD(vecp, oper_choice);

    Eigen::Map<const Eigen::VectorXd> vecm(vecp.data(), vecp.size());
    auto diag = vecm.asDiagonal();
    _eigenMatrix = diag * _eigenMatrix * diag;
  }
  else
  {
    cs_matvecnorm_inplace(_csMatrix, vec.data(), oper_choice);
  }
}

void MatrixSparse::prodNormMatInPlace(const MatrixSparse &a, const VectorDouble& vec, bool transpose)
{
  if (!_checkLink(getNRows(), getNCols(), a.getNRows(), a.getNCols(), transpose,
                 vec.size(), 1, false)) return;

  if (isFlagEigen() && a.isFlagEigen())
  {
    if (transpose)
    {
      if (vec.empty())
        _eigenMatrix = a._eigenMatrix.transpose() * a._eigenMatrix;
      else
      {
        Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), vec.size());
        _eigenMatrix = a._eigenMatrix.transpose() * vecm.asDiagonal() * a._eigenMatrix;
      }
    }
    else
    {
      if (vec.empty())
        _eigenMatrix = a._eigenMatrix * a._eigenMatrix.transpose();
      else
      {
        Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), vec.size());
        _eigenMatrix = a._eigenMatrix * vecm.asDiagonal() * a._eigenMatrix.transpose();
      }
    }
  }
  else
  {
    cs* res = nullptr;
    if (vec.empty())
      res = cs_prod_norm_single((transpose) ? 1 : 2, a._csMatrix);
    else
      res = cs_prod_norm_diagonal((transpose) ? 1 : 2, a._csMatrix, vec);
    _csMatrix = cs_spfree2(_csMatrix);
    _csMatrix = res;
  }
}

void MatrixSparse::prodNormMatMatInPlace(const MatrixSparse &a,
                                         const MatrixSparse &m,
                                         bool transpose)
{
  if (!_checkLink(a.getNRows(), a.getNCols(), transpose,
                  m.getNRows(), m.getNCols(), false,
                  a.getNRows(), a.getNCols(), !transpose)) return;

  if (isFlagEigen() && a.isFlagEigen() && m.isFlagEigen())
  {
    if (transpose)
    {
      _eigenMatrix = a._eigenMatrix.transpose() * m._eigenMatrix * a._eigenMatrix;
    }
    else
    {
      _eigenMatrix = a._eigenMatrix * m._eigenMatrix * a._eigenMatrix.transpose();
    }
  }
  else
  {
    cs* res = cs_prod_norm((transpose) ? 1 : 2, m._csMatrix, a._csMatrix);
    _csMatrix = cs_spfree2(_csMatrix);
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
void MatrixSparse::addMatInPlace(const MatrixSparse& y, double cx, double cy)
{
  if (!_checkLink(y.getNRows(), y.getNCols(), false)) return;

  if (isFlagEigen() && y.isFlagEigen())
  {
    _eigenMatrix = cx * _eigenMatrix + cy * y._eigenMatrix;
  }
  else
  {
    cs *res = cs_add(_csMatrix, y._csMatrix, cx, cy);
    _csMatrix = cs_spfree2(_csMatrix);
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

  if (! isSquare())
    my_throw("Invert method is limited to Square Matrices");
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
    sstr << "  (not using Eigen Library)" << std::endl;
    sstr << "- Sparse Format" << std::endl;
    sstr << toMatrix(String(), *this);
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
    int nrow = getNRows();
    int ncol = getNCols();
    if (nrow > 0 && ncol > 0)
    {
      cs* local = cs_spalloc2(0, 0, 1, 1, 1);
      cs_entry2(local, nrow - 1, ncol - 1, 0.);
      _csMatrix = cs_triplet2(local);
      local = cs_spfree2(local);
    }
  }
}

void MatrixSparse::_deallocate()
{
  if (isFlagEigen())
  {
    _eigenMatrix.data().squeeze();
  }
  else
  {
    _csMatrix = cs_spfree2(_csMatrix);
  }
  delete _factor;
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
NF_Triplet MatrixSparse::getMatrixToTriplets(bool flag_from_1) const
{
  NF_Triplet NF_T = tripletInit(flag_from_1);
  if (isFlagEigen())
  {
    int ecr = 0;
    int shift = (flag_from_1) ? 1 : 0;
    for (int k = 0; k < _eigenMatrix.outerSize(); ++k)
      for (Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix, k); it; ++it)
      {
        double value = it.value();
        if (ABS(value) <= EPSILON10) continue;
        NF_T.rows.push_back(it.row() + shift);
        NF_T.cols.push_back(it.col() + shift);
        NF_T.values.push_back(value);
        ecr++;
      }
    NF_T.number = ecr;
  }
  else
  {
    NF_T = csToTriplet(_csMatrix, flag_from_1);
  }
  return NF_T;
}

void MatrixSparse::_clear()
{
  _setNRows(0);
  _setNCols(0);
  _allocate();
}

int MatrixSparse::_getIndexToRank(int irow,int icol) const
{
  DECLARE_UNUSED(irow);
  DECLARE_UNUSED(icol);
  _forbiddenForSparse("_getIndexToRank");
  return ITEST;
}

MatrixSparse* createFromAnyMatrix(const AMatrix* matin, bool flag_from_1)
{
  // Create the output Sparse Matrix

  MatrixSparse* matout = new MatrixSparse(matin->getNRows(),matin->getNCols());

  // Create the triplet structure from the non-zero terms of source matrix

  NF_Triplet T = matin->getMatrixToTriplets(flag_from_1);

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

/****************************************************************************/
/*!
 **  Creating the color coding for mesh nodes (in Eigen library)
 **
 ** \return  Error return code
 **
 ** \param[in] start    Starting value for colors ranking (usually 0 or 1)
 **
 *****************************************************************************/
VectorInt MatrixSparse::_eigen_color_coding(int start)
{
  int ncolor = 0;
  int nmesh = getNCols();

  /* Core allocation */

  VectorInt colors(nmesh, ITEST);
  VectorInt temp(nmesh);

  /* Loop on the nodes of the mesh */

  for (int imesh = 0; imesh < nmesh; imesh++)
  {
    temp.fill(0);
    for (Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix,imesh); it; ++it)
    {
      if (ABS(it.value()) < EPSILON10) continue;
      if (!IFFFF(colors[it.row()])) temp[colors[it.row()] - 1]++;
    }
    _updateColors(temp, colors, imesh, &ncolor);
  }

  /* Update the colors ranking fixing the starting value */

  VH::addConstant(colors, start-1);
  return colors;
}

VectorInt MatrixSparse::colorCoding(int start)
{
  if (isFlagEigen())
    return _eigen_color_coding(start);
  else
    return cs_color_coding(_csMatrix, start);
}

/**
 * Convert Sparse Matrix to Triplet (Eigen format)
 * @param shiftRow Shift the row number (if > 0)
 * @param shiftCol Shift the Column number (if > 0)
 * @return
 */
std::vector<EigT> MatrixSparse::toTriplet(int shiftRow, int shiftCol) const
{
  std::vector<EigT> tripletList;
  if (!isFlagEigen())
  {
    messerr("Method 'toTriplet' is restricted to the Eigen version");
    return tripletList;
  }

  int nnz = _eigenMatrix.nonZeros();
  tripletList.reserve(nnz);

  for (int i = 0; i < _eigenMatrix.outerSize(); i++)
    for (typename Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix, i); it; ++it)
      tripletList.emplace_back(it.row() + shiftRow, it.col() + shiftCol, it.value());
  return tripletList;
}

MatrixSparse* MatrixSparse::glue(const MatrixSparse *A1,
                                 const MatrixSparse *A2,
                                 bool flagShiftRow,
                                 bool flagShiftCol)
{
  if (A1->isFlagEigen())
  {
    int shiftRow =  (flagShiftRow) ? A1->getNRows() : 0;
    int shiftCol =  (flagShiftCol) ? A1->getNCols() : 0;

    // Create the two triplet lists
    std::vector<EigT> T1 = A1->toTriplet();
    std::vector<EigT> T2 = A2->toTriplet(shiftRow, shiftCol);

    // Concatenate the two triplet lists
    T1.insert( T1.end(), T2.begin(), T2.end() );

    // Create the new Eigen matrix from the resulting triplet list
    int nrow = (flagShiftRow) ? A1->getNRows() + A2->getNRows() : MAX(A1->getNRows(), A2->getNRows());
    int ncol = (flagShiftCol) ? A1->getNCols() + A2->getNCols() : MAX(A1->getNCols(), A2->getNCols());
    Eigen::SparseMatrix<double> local(nrow,ncol);
    local.setFromTriplets(T1.begin(), T1.end());

    MatrixSparse* mat = new MatrixSparse(nrow, ncol, 1);
    mat->setEigenMatrix(local);
    return mat;
  }
  else
  {
    cs* local = cs_glue(A1->getCS(), A2->getCS(), flagShiftRow, flagShiftCol);
    MatrixSparse* mat = new MatrixSparse(local, 0);
    local = cs_spfree2(local);
    return mat;
  }
}

/* Extract a sparse submatrix */
/* 'rank_rows' and 'rank_cols' must have same dimension as C */
/* The arrays 'rank_rows' and 'rank_cols' may be absent */
/* Their value gives the rank of the saved element or -1 */
MatrixSparse* MatrixSparse::extractSubmatrixByRanks(const VectorInt &rank_rows,
                                                    const VectorInt &rank_cols)
{
  int old_row, old_col, new_row, new_col;

  NF_Triplet NF_Tin = getMatrixToTriplets(false);
  NF_Triplet NF_Tout = tripletInit(0);

  /* Fill the new sparse triplet */

  for (int i = 0; i < NF_Tin.number; i++)
  {
    old_row = NF_Tin.rows[i];
    old_col = NF_Tin.cols[i];
    new_row = (!rank_rows.empty()) ? rank_rows[old_row] : old_row;
    new_col = (!rank_cols.empty()) ? rank_cols[old_col] : old_col;
    if (new_row < 0 || new_col < 0) continue;
    tripletAdd(NF_Tout, new_row, new_col, NF_Tin.values[i]);
  }

  return MatrixSparse::createFromTriplet(NF_Tout);
}

/* Extract a sparse submatrix */
/* The array 'colors' has the same dimension as C */
/* The element of 'C' must be kept if: */
/* - the color of its row number is equal to 'ref_color' if 'row_ok'==TRUE */
/*   or different if 'row_ok'== FALSE */
/* and if */
/* - the color of its column number is equal to 'ref-color' if 'col_ok'==TRUE*/
/*   or different if 'col_ok'== FALSE */
MatrixSparse* MatrixSparse::extractSubmatrixByColor(const VectorInt &colors,
                                                    int ref_color,
                                                    bool row_ok,
                                                    bool col_ok)
{
  /* Convert the contents of the sparse matrix into columns */

  NF_Triplet NF_Tin = getMatrixToTriplets(false);

  /* Initialize the output matrix */

  NF_Triplet NF_Tout = tripletInit(0);

  /* Core allocation */

  int n = getNCols();
  VectorInt u_row(n);
  VectorInt u_col(n);

  int ir = 0;
  for (int i = 0; i < n; i++)
  {
    u_row[i] = -1;
    if ( row_ok && colors[i] != ref_color) continue;
    if (!row_ok && colors[i] == ref_color) continue;
    u_row[i] = ir++;
  }

  int ic = 0;
  for (int i = 0; i < n; i++)
  {
    u_col[i] = -1;
    if ( col_ok && colors[i] != ref_color) continue;
    if (!col_ok && colors[i] == ref_color) continue;
    u_col[i] = ic++;
  }

  /* Fill the new sparse triplet */

  for (int i = 0; i < NF_Tin.number; i++)
  {
    ir = u_row[NF_Tin.rows[i]];
    ic = u_col[NF_Tin.cols[i]];
    if (ir < 0 || ic < 0) continue;
    tripletAdd(NF_Tout, ir, ic, NF_Tin.values[i]);
  }

  return MatrixSparse::createFromTriplet(NF_Tout);
}
