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
#include "Basic/AStringable.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Matrix/LinkMatrixSparse.hpp"
#include "Matrix/NF_Triplet.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Basic/WarningMacro.hpp"

#include <Eigen/src/Core/Map.h>
#include <Eigen/src/Core/Matrix.h>
#include <iostream>

DISABLE_WARNING_PUSH
DISABLE_WARNING_COND_EXPR_CONSTANT
DISABLE_WARNING_UNUSED_BUT_SET_VARIABLE
#include <Eigen/SparseCholesky>
#include <omp.h>
DISABLE_WARNING_POP

#include <csparse_f.h>

/**
 * This variable switches ON/OFF the ability to use Eigen library for Algebra
 */
static bool globalFlagEigen = true;

MatrixSparse::MatrixSparse(int nrow, int ncol, int opt_eigen)
  : AMatrix(nrow, ncol),
    _csMatrix(nullptr),
    _eigenMatrix(),
    _flagEigen(false)
{
  _flagEigen = _defineFlagEigen(opt_eigen);
  _allocate();
}

#ifndef SWIG
MatrixSparse::MatrixSparse(const cs *A)
  : AMatrix(cs_get_nrow(A), cs_get_ncol(A)),
    _csMatrix(nullptr),
    _eigenMatrix(),
    _flagEigen(false)
{
  _csMatrix = cs_duplicate(A);
}
#endif

MatrixSparse::MatrixSparse(const MatrixSparse &m)
    : AMatrix(m),
      _csMatrix(nullptr),
      _eigenMatrix(),
      _flagEigen(m._flagEigen)
{
  if (isFlagEigen())
    _eigenMatrix = m._eigenMatrix;
  else
    _csMatrix = cs_duplicate(m._csMatrix);
}

MatrixSparse& MatrixSparse::operator=(const MatrixSparse& m)
{
  if (this != &m)
  {
    AMatrix::operator=(m);
    if (!m.empty())
    {
      _flagEigen = m._flagEigen;
      if (_flagEigen)
        _eigenMatrix = m._eigenMatrix;
      else
      {
        _csMatrix = cs_duplicate(m._csMatrix);
      }
    }
  }
  return *this;
}

MatrixSparse::~MatrixSparse()
{
  _deallocate();
}

void MatrixSparse::resetFromValue(int nrows, int ncols, double value)
{
  DECLARE_UNUSED(nrows);
  DECLARE_UNUSED(ncols);
  DECLARE_UNUSED(value);
  _forbiddenForSparse("resetFromValue");
}

void MatrixSparse::resetFromArray(int nrows, int ncols, const double* tab, bool byCol)
{
  DECLARE_UNUSED(nrows);
  DECLARE_UNUSED(ncols);
  DECLARE_UNUSED(tab);
  DECLARE_UNUSED(byCol);
  _forbiddenForSparse("resetFromArray");
}

void MatrixSparse::resetFromVD(int nrows, int ncols, const VectorDouble &tab, bool byCol)
{
  DECLARE_UNUSED(nrows);
  DECLARE_UNUSED(ncols);
  DECLARE_UNUSED(tab);
  DECLARE_UNUSED(byCol);
  _forbiddenForSparse("resetFromVD");
}

void MatrixSparse::resetFromVVD(const VectorVectorDouble& tab, bool byCol)
{
  DECLARE_UNUSED(tab);
  DECLARE_UNUSED(byCol);
  _forbiddenForSparse("resetFromVVD");
}

void MatrixSparse::resetFromTriplet(const NF_Triplet& NF_T)
{
  delete _csMatrix;

  if (isFlagEigen())
  {
    _eigenMatrix = NF_T.buildEigenFromTriplet();
    _setNRows(_eigenMatrix.rows());
    _setNCols(_eigenMatrix.cols());
  }
  else
  {
    _csMatrix = NF_T.buildCsFromTriplet();
    _setNRows(_csMatrix->m);
    _setNCols(_csMatrix->n);
  }
}

void MatrixSparse::fillRandom(int seed, double zeroPercent)
{
  law_set_random_seed(seed);

  int nrow = getNRows();
  int ncol = getNCols();
  NF_Triplet NF_T;
  for (int irow = 0; irow < nrow; irow++)
    for (int icol = 0; icol < ncol; icol++)
    {
      if (law_uniform(0., 1.) < zeroPercent) continue;
      NF_T.add(irow, icol, law_gaussian());
    }
  NF_T.force(nrow,ncol);
  resetFromTriplet(NF_T);
}

void MatrixSparse::_transposeInPlace()
{
  if (isFlagEigen())
  {
    Eigen::SparseMatrix<double> temp;
    temp = _eigenMatrix.transpose();
    _eigenMatrix.swap(temp);
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
 * @param flagCheck When True, check the consistency of arguments
 */
void MatrixSparse::setColumn(int icol, const VectorDouble& tab, bool flagCheck)
{
  int nrows = getNRows();
  if (flagCheck)
  {
    if (! _isColumnValid(icol)) return;
    if (! _isColumnSizeConsistent(tab)) return;
  }
  if (isFlagEigen())
  {
    for (int irow = 0; irow < nrows; irow++)
      _eigenMatrix.coeffRef(irow, icol) = tab[irow];
  }
  else
    AMatrix::setColumn(icol, tab);
}

/**
 * Fill a row of an already existing Sparse matrix, using 'tab' as entry
 * The input 'tab' corresponds to the whole row contents
 * @param irow Row rank
 * @param tab  Vector containing the information (Dimension: ncols)
 * @param flagCheck True if the validity check must be performed
 *
 * @warning: This method only copies the values at the non-zero existing entries
 */
void MatrixSparse::setRow(int irow, const VectorDouble& tab, bool flagCheck)
{
  int ncols = getNCols();
  if (flagCheck)
  {
    if (! _isRowValid(irow)) return;
    if (! _isRowSizeConsistent(tab)) return;
  }
  if (isFlagEigen())
  {
    for (int icol = 0; icol < ncols; icol++)
      _eigenMatrix.coeffRef(irow, icol) = tab[icol];
  }
  else
  {
    AMatrix::setRow(irow, tab);
  }
}

void MatrixSparse::setDiagonal(const VectorDouble& tab, bool flagCheck)
{
  if (! isSquare())
    my_throw("This function is only valid for Square matrices");
  if (flagCheck)
  {
    if (! _isRowSizeConsistent(tab)) return;
  }

  if (isFlagEigen())
  {
    Eigen::Map<const Eigen::VectorXd> vecm(tab.data(), tab.size());
    _eigenMatrix = vecm.asDiagonal();
  }
  else
  {
    AMatrix::setDiagonal(tab, flagCheck);
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
    AMatrix::setDiagonalToConstant(value);
  }
}

/*! Gets the value for rank 'rank' */
double MatrixSparse::_getValueByRank(int rank) const
{
  DECLARE_UNUSED(rank);
  _forbiddenForSparse("_getValueByRank");
  return TEST;
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
  _forbiddenForSparse("_setValueByRank");
}

void MatrixSparse::setValue(int irow, int icol, double value, bool flagCheck)
{
  if (flagCheck && ! _isIndexValid(irow, icol)) return;
  if (isFlagEigen())
  {
    _eigenMatrix.coeffRef(irow, icol) = value;
  }
  else
  {
    cs_set_value(_csMatrix, irow, icol, value);
  }
}

void MatrixSparse::updValue(int irow,
                            int icol,
                            const EOperator &oper,
                            double value,
                            bool flagCheck)
{
  if (flagCheck && ! _isIndexValid(irow, icol)) return;
  if (isFlagEigen())
  {
    double newval = modifyOperator(oper, _eigenMatrix.coeff(irow, icol), value);
    _eigenMatrix.coeffRef(irow, icol) = newval;
  }
  else
  {
    if (! _isIndexValid(irow, icol)) return;
    double newval = modifyOperator(oper, cs_get_value(_csMatrix, irow, icol), value);
    cs_set_value(_csMatrix, irow, icol, newval);
  }
}

int MatrixSparse::_getMatrixPhysicalSize() const
{
  if (isFlagEigen())
    return _eigenMatrix.nonZeros();
  return cs_nnz(_csMatrix);
}

/**
 * @param value Constant value used for filling 'this'
 */
void MatrixSparse::fill(double value)
{
  int nrow = getNRows();
  int ncol = getNCols();
  NF_Triplet NF_T;
  for (int irow = 0; irow < nrow; irow++)
    for (int icol = 0; icol < ncol; icol++)
      NF_T.add(irow, icol, value);

  resetFromTriplet(NF_T);
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
    VectorDouble y(transpose ? getNRows() : getNCols());
    Eigen::Map<Eigen::VectorXd> ym(y.data(), y.size());
    if (transpose)
      ym = xm.transpose() * _eigenMatrix.transpose();
    else
      ym = xm.transpose() * _eigenMatrix;
    return y;
  }
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

/*! Perform y = 'this' %*% x */
VectorDouble MatrixSparse::prodMatVec(const VectorDouble& x, bool transpose) const
{
  if (isFlagEigen())
  {
    Eigen::Map<const Eigen::VectorXd> xm(x.data(), x.size());
    VectorDouble y(transpose ? getNCols() : getNRows());
    Eigen::Map<Eigen::VectorXd> ym(y.data(), y.size());
    if (transpose)
      ym = _eigenMatrix.transpose() * xm;
    else
      ym = _eigenMatrix * xm;
    return y;
  }
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

/*! Perform y += 'this' %*% x */
void MatrixSparse::addProdMatVecInPlaceToDest(const constvect in,
                                              vect out,
                                              bool transpose) const
{
  Eigen::Map<const Eigen::VectorXd> inm(in.data(),in.size());
  Eigen::Map<Eigen::VectorXd>       outm(out.data(),out.size());
  if (isFlagEigen())
  {
    if (transpose)
      outm += _eigenMatrix.transpose() * inm;
    else
      outm += _eigenMatrix * inm;
  }
  else
  {
    if (transpose)
    {
      int ncol = getNCols();
      cs_vector_addToDest_tMx(_csMatrix, ncol, in.data(), out.data());
    }
    else
    {
      int nrow = getNRows();
      cs_vector_addToDest_Mx(_csMatrix, nrow, in.data(), out.data());
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
          if (isZero(values[lec])) continue;
          (void) cs_entry2(local, irow, icol, values[lec]);
        }
    }
    else
    {
      for (int irow = 0; irow < getNRows(); irow++)
        for (int icol = 0; icol < getNCols(); icol++, lec++)
        {
          if (isZero(values[lec])) continue;
          (void) cs_entry2(local, irow, icol, values[lec]);
        }
    }
    _csMatrix = cs_triplet2(local);
    local = cs_spfree2(local);
  }
}
#endif

MatrixSparse* MatrixSparse::create(const MatrixSparse* mat)
{
  return new MatrixSparse(*mat);
}

MatrixSparse* MatrixSparse::create(int nrow, int ncol)
{
  return new MatrixSparse(nrow, ncol);
}

MatrixSparse* MatrixSparse::createFromTriplet(const NF_Triplet &NF_T,
                                              int nrow,
                                              int ncol,
                                              int opt_eigen)
{
  // If 'nrow' a  nd 'ncol' are not defined, derive them from NF_T
  if (nrow <= 0 || ncol <= 0)
  {
    nrow = NF_T.getNRows() + 1;
    ncol = NF_T.getNCols() + 1;
  }
  MatrixSparse* mat = new MatrixSparse(nrow, ncol, opt_eigen);

  mat->resetFromTriplet(NF_T);

  return mat;
}

MatrixSparse* MatrixSparse::addMatMat(const MatrixSparse *x,
                                      const MatrixSparse *y,
                                      double cx, double cy)
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
  return cs_exist(_csMatrix, irow, icol);
}

void MatrixSparse::addValue(int row, int col, double value)
{
  if (ABS(value) <= EPSILON10) return;
  if (isFlagEigen())
    _eigenMatrix.coeffRef(row, col) += value;
  else
    cs_add_value(_csMatrix, row, col, value);
}

double MatrixSparse::getValue(int row, int col, bool flagCheck) const
{
  if (flagCheck && ! _isIndexValid(row, col)) return TEST;
  if (isFlagEigen())
    return _eigenMatrix.coeff(row, col);
  return cs_get_value(_csMatrix, row, col);
}

double MatrixSparse::L1Norm() const
{
  if (isFlagEigen())
    return (Eigen::RowVectorXd::Ones(_eigenMatrix.rows()) * _eigenMatrix.cwiseAbs()).maxCoeff();
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
    VectorDouble diag(std::min(getNCols(), getNRows()));
    Eigen::Map<Eigen::VectorXd> ym(diag.data(), diag.size());
    ym = _eigenMatrix.diagonal();
    VH::transformVD(diag, oper_choice);
    return diag;
  }
  VectorDouble diag = csd_extract_diag_VD(_csMatrix, 1);
  VH::transformVD(diag, oper_choice);
  return diag;
}

int MatrixSparse::addVecInPlaceEigen(const Eigen::Map<const Eigen::VectorXd>& xm,
                                     Eigen::Map<Eigen::VectorXd>& ym) const
{
  if (isFlagEigen())
  {
    ym = _eigenMatrix * xm + ym;
    return 0;
  }
  return (!cs_gaxpy(_csMatrix, xm.data(), ym.data()));
}

int MatrixSparse::addVecInPlace(const constvect xm, vect ym) const
{
  if (isFlagEigen())
  {
    Eigen::Map<const Eigen::VectorXd> xmm(xm.data(), xm.size());
    Eigen::Map<Eigen::VectorXd> ymm(ym.data(), ym.size());
    ymm = _eigenMatrix * xmm + ymm;
    return 0;
  }
  return (!cs_gaxpy(_csMatrix, xm.data(), ym.data()));
}

int MatrixSparse::addVecInPlaceVD(const VectorDouble& x, VectorDouble& y) const
{
  if (isFlagEigen())
  {
    Eigen::Map<const Eigen::VectorXd> xm(x.data(), x.size());
    Eigen::Map<Eigen::VectorXd> ym(y.data(), y.size());
    ym = _eigenMatrix * xm + ym;
    return 0;
  }
  return (!cs_gaxpy(_csMatrix, x.data(), y.data()));
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
  return cs_scale(_csMatrix);
}

/**
 *
 * @param v Add a scalar value to all terms of the current matrix
 */
void MatrixSparse::addScalar(double v)
{
  if (isZero(v)) return;
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
          setValue(irow, icol, getValue(irow, icol, false) + v, false);
      }
  }
}

/**
 *
 * @param v Add constant value to the diagonal of the current Matrix
 */
void MatrixSparse::addScalarDiag(double v)
{
  if (isZero(v)) return;

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
  if (isOne(v)) return;
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

void MatrixSparse::_addProdMatVecInPlaceToDestPtr(const double *x, double *y, bool transpose) const
{
  if (isFlagEigen())
  {
    if (transpose)
    {
      Eigen::Map<const Eigen::VectorXd> xm(x, getNRows());
      Eigen::Map<Eigen::VectorXd> ym(y, getNCols());
      ym += _eigenMatrix.transpose() * xm;
    }
    else
    {
      Eigen::Map<const Eigen::VectorXd> xm(x, getNCols());
      Eigen::Map<Eigen::VectorXd> ym(y, getNRows());
      ym += _eigenMatrix * xm;
    }
  }
  else
  {
    if (transpose)
      cs_vector_addToDest_tMx(_csMatrix, getNCols(), x, y);
    else
      cs_vector_addToDest_Mx(_csMatrix, getNRows(), x, y);
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

MatrixSparse* prodNormMatMat(const MatrixSparse* a,
                             const MatrixSparse* m,
                             bool transpose)
{
  int nrow = (transpose) ? a->getNCols() : a->getNRows();
  int ncol = (transpose) ? a->getNRows() : a->getNCols();
  MatrixSparse *mat = new MatrixSparse(nrow, ncol, a->isFlagEigen() ? 1 : 0);
  mat->prodNormMatMatInPlace(a, m, transpose);
  return mat;
}

MatrixSparse* prodNormMat(const MatrixSparse* a, const VectorDouble& vec, bool transpose)
{
  int nsym = (transpose) ? a->getNCols() : a->getNRows();
  MatrixSparse *mat = new MatrixSparse(nsym, nsym, a->isFlagEigen() ? 1 : 0);
  mat->prodNormMatVecInPlace(a, vec, transpose);
  return mat;
}

MatrixSparse* prodNormDiagVec(const MatrixSparse* a,
                              const VectorDouble &vec,
                              int oper_choice)
{
  int nrow = a->getNRows();
  int ncol = a->getNCols();
  MatrixSparse *mat = new MatrixSparse(nrow, ncol, a->isFlagEigen() ? 1 : 0);

  if (a->isFlagEigen())
  {
    // Perform the transformation of the input vector
    VectorDouble vecp = vec;
    VH::transformVD(vecp, oper_choice);

    Eigen::Map<const Eigen::VectorXd> vecm(vecp.data(), vecp.size());
    auto diag = vecm.asDiagonal();
    mat->setEigenMatrix(diag * a->getEigenMatrix() * diag);
  }
  else
  {
    cs* local = cs_matvecnorm(a->getCS(), vec.data(), oper_choice);
    mat->setCS(local);
    cs_spfree2(local);
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
    messerr("Matrix dimension (%d) does not match vector dimension (%d)",
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

void MatrixSparse::prodNormMatVecInPlace(const MatrixSparse* a, const VectorDouble& vec, bool transpose)
{
  if (!_checkLink(getNRows(), getNCols(), transpose, a->getNRows(), a->getNCols(),
                  false, vec.size(), 1, false)) return;

  if (isFlagEigen() && a->isFlagEigen())
  {
    if (transpose)
    {
      if (vec.empty())
        _eigenMatrix = a->_eigenMatrix.transpose() * a->_eigenMatrix;
      else
      {
        Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), vec.size());
        _eigenMatrix = a->_eigenMatrix.transpose() * vecm.asDiagonal() * a->_eigenMatrix;
      }
    }
    else
    {
      if (vec.empty())
        _eigenMatrix = a->_eigenMatrix * a->_eigenMatrix.transpose();
      else
      {
        Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), vec.size());
        _eigenMatrix = a->_eigenMatrix * vecm.asDiagonal() * a->_eigenMatrix.transpose();
      }
    }
  }
  else
  {
    cs* res = nullptr;
    if (vec.empty())
      res = cs_prod_norm_single((transpose) ? 1 : 2, a->_csMatrix);
    else
      res = cs_prod_norm_diagonal((transpose) ? 1 : 2, a->_csMatrix, vec);
    _csMatrix = cs_spfree2(_csMatrix);
    _csMatrix = res;
  }
}


#ifndef SWIG
/*! Returns a pointer to the Sparse storage */
const cs* MatrixSparse::getCS() const
{
  return _csMatrix;
}
void MatrixSparse::setCS(cs* cs)
{
  _csMatrix = cs_duplicate(cs);
}
void MatrixSparse::freeCS()
{
  _csMatrix = cs_spfree2(_csMatrix);
}
/*! Temporary function to get the CS contents of Sparse Matrix */
cs* MatrixSparse::getCSUnprotected() const
{
  return _csMatrix;
}
#endif

void MatrixSparse::prodNormMatMatInPlace(const MatrixSparse* a,
                                         const MatrixSparse* m,
                                         bool transpose)
{
  if (!_checkLink(a->getNRows(), a->getNCols(), transpose,
                  m->getNRows(), m->getNCols(), false,
                  a->getNRows(), a->getNCols(), !transpose)) return;

  if (isFlagEigen() && a->isFlagEigen() && m->isFlagEigen())
  {
    if (transpose)
    {
      _eigenMatrix = (a->_eigenMatrix.transpose() * m->_eigenMatrix) * a->_eigenMatrix;
    }
    else
    {
      _eigenMatrix = (a->_eigenMatrix * m->_eigenMatrix) * a->_eigenMatrix.transpose();
    }
  }
  else
  {
    cs* res = cs_prod_norm((transpose) ? 1 : 2, m->_csMatrix, a->_csMatrix);
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
}

void MatrixSparse::_forbiddenForSparse(const String& func)
{
  messerr("Problem with Function: %s",func.c_str());
  messerr("This function is not available in Sparse Matrix");
}

void MatrixSparse::dumpElements(const String& title, int ifrom, int ito)
{
  DECLARE_UNUSED(title);
  DECLARE_UNUSED(ifrom);
  DECLARE_UNUSED(ito);
  messerr("This method is not implemented for Sparse Matrix");
}

/**
 * From a matrix of any type, creates the triplet
 * (specific format for creating efficiently a Sparse matrix)
 * It only takes the only non-zero elements of the matrix
 */
NF_Triplet MatrixSparse::getMatrixToTriplet(int shiftRow, int shiftCol) const
{
  if (isFlagEigen())
  {
    return NF_Triplet::createFromEigen(_eigenMatrix, shiftRow, shiftCol);
  }
  return NF_Triplet::createFromCs(_csMatrix, shiftRow, shiftCol);
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

MatrixSparse* createFromAnyMatrix(const AMatrix* matin, int opt_eigen)
{
  return MatrixSparse::createFromTriplet(matin->getMatrixToTriplet(),
                                         matin->getNRows(),
                                         matin->getNCols(),
                                         opt_eigen);
}

void setUpdateNonZeroValue(int status)
{
  cs_set_status_update_nonzero_value(status);
}

int getUpdateNonZeroValue()
{
  return cs_get_status_update_nonzero_value();
}

int MatrixSparse::_eigen_findColor(int imesh,
                                   int ncolor,
                                   VectorInt &colors,
                                   VectorInt &temp) const
{
  temp.fill(0);

  /* Checks the colors of the connected nodes */

  for (Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix,imesh); it; ++it)
  {
    if (isZero(it.value())) continue;
    int irow = it.row();
    if (!IFFFF(colors[irow])) temp[colors[irow] - 1]++;
  }

  /* Look for a free color */

  for (int j = 0; j < ncolor; j++)
  {
    if (temp[j] == 0) return (j + 1);
  }
  return (-1);
}

VectorInt MatrixSparse::colorCoding() const
{
  int next_col = 0;
  int ncol = 0;
  int nmesh = getNCols();

  /* Core allocation */

  VectorInt colors(nmesh, ITEST);
  VectorInt temp(nmesh);

  /* Loop on the nodes of the mesh */

  for (int imesh = 0; imesh < nmesh; imesh++)
  {
    if (isFlagEigen())
      next_col = _eigen_findColor(imesh, ncol, colors, temp);
    else
      next_col = _cs_findColor(_csMatrix, imesh, ncol, colors, temp);

    if (next_col < 0)
    {
      ncol++;
      colors[imesh] = ncol;
    }
    else
    {
      colors[imesh] = next_col;
    }
  }
  return colors;
}

void MatrixSparse::glueInPlace(MatrixSparse* A1,
                               const MatrixSparse* A2,
                               bool flagShiftRow,
                               bool flagShiftCol)
{
  int shiftRow = (flagShiftRow) ? A1->getNRows() : 0;
  int shiftCol = (flagShiftCol) ? A1->getNCols() : 0;
  NF_Triplet T1 = A1->getMatrixToTriplet();
  NF_Triplet T2 = A2->getMatrixToTriplet(shiftRow, shiftCol);

  // Concatenate the two triplet lists
  T1.appendInPlace(T2);

  // Create the new matrix from the resulting triplet list
  int nrow = (flagShiftRow) ? A1->getNRows() + A2->getNRows() : MAX(A1->getNRows(), A2->getNRows());
  int ncol = (flagShiftCol) ? A1->getNCols() + A2->getNCols() : MAX(A1->getNCols(), A2->getNCols());

  A1->resetFromTriplet(T1);
  A1->_setNRows(nrow);
  A1->_setNCols(ncol);
}

MatrixSparse* MatrixSparse::glue(const MatrixSparse *A1,
                                 const MatrixSparse *A2,
                                 bool flagShiftRow,
                                 bool flagShiftCol)
{
  int shiftRow = (flagShiftRow) ? A1->getNRows() : 0;
  int shiftCol = (flagShiftCol) ? A1->getNCols() : 0;

  // Create the two triplet lists
  NF_Triplet T1 = A1->getMatrixToTriplet();
  NF_Triplet T2 = A2->getMatrixToTriplet(shiftRow, shiftCol);

  // Concatenate the two triplet lists
  T1.appendInPlace(T2);

  // Create the new matrix from the resulting triplet list
  int nrow = (flagShiftRow) ? A1->getNRows() + A2->getNRows() : MAX(A1->getNRows(), A2->getNRows());
  int ncol = (flagShiftCol) ? A1->getNCols() + A2->getNCols() : MAX(A1->getNCols(), A2->getNCols());

  return MatrixSparse::createFromTriplet(T1, nrow, ncol, A1->isFlagEigen() ? 1 : 0);
}

/* Extract a sparse sub-matrix */
/* 'rank_rows' and 'rank_cols' must have same dimension as C */
/* The arrays 'rank_rows' and 'rank_cols' may be absent */
/* Their value gives the rank of the saved element or -1 */
MatrixSparse* MatrixSparse::extractSubmatrixByRanks(const VectorInt &rank_rows,
                                                    const VectorInt &rank_cols) const
{
  int old_row, old_col, new_row, new_col;

  NF_Triplet NF_Tin = getMatrixToTriplet();
  NF_Triplet NF_Tout;

  /* Fill the new sparse triplet */

  for (int i = 0; i < NF_Tin.getNumber(); i++)
  {
    old_row = NF_Tin.getRow(i);
    old_col = NF_Tin.getCol(i);
    new_row = (!rank_rows.empty()) ? rank_rows[old_row] : old_row;
    new_col = (!rank_cols.empty()) ? rank_cols[old_col] : old_col;
    if (new_row < 0 || new_col < 0) continue;
    NF_Tout.add(new_row, new_col, NF_Tin.getValue(i));
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

  NF_Triplet NF_Tin = getMatrixToTriplet();

  /* Initialize the output matrix */

  NF_Triplet NF_Tout;

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

  for (int i = 0; i < NF_Tin.getNumber(); i++)
  {
    ir = u_row[NF_Tin.getRow(i)];
    ic = u_col[NF_Tin.getCol(i)];
    if (ir < 0 || ic < 0) continue;
    NF_Tout.add(ir, ic, NF_Tin.getValue(i));
  }

  return MatrixSparse::createFromTriplet(NF_Tout,0,0, isFlagEigen() ? 1 : 0);
}

/**
 * Define the use of Eigen Library according to the value of input argulent 'opt_eigen'
 * @param opt_eigen Choice: 0: Do not use Eigen library; 1; Use Eigen library; -1: use global environment
 * @return Option for using the Eigen library
 */
bool MatrixSparse::_defineFlagEigen(int opt_eigen)
{
  // Dispatch

  if (opt_eigen == 1) return true;
  if (opt_eigen == 0) return false;
  return globalFlagEigen;
}

/**
 * Modify the parameter for using EIGEN library or not.
 * Warning: this must be performed very early in the script in order to forbid mixing two different styles.
 * @param flagEigen True if EIGEN library must be used; False otherwise (cs is used)
 */
void setGlobalFlagEigen(bool flagEigen)
{
  globalFlagEigen = flagEigen;
}

bool isGlobalFlagEigen()
{
  return globalFlagEigen;
}

void MatrixSparse::gibbs(int iech,
                         const VectorDouble& zcur,
                         double* yk,
                         double* sk)
{
  if (isFlagEigen())
  {
    *yk = 0.;
    for (Eigen::SparseMatrix<double>::InnerIterator it(_eigenMatrix, iech); it;
         ++it)
    {
      double coeff = it.valueRef();
      if (ABS(coeff) <= 0.) continue;
      int jech = it.row();

      if (iech == jech)
        *sk = coeff;
      else
        *yk -= coeff * zcur[jech];
    }
  }
  else
  {
    cs_gibbs(_csMatrix, iech, zcur, yk, sk);
  }

  // Returned arguments
  (*yk) /= (*sk);
  (*sk) = sqrt(1. / (*sk));
}

int MatrixSparse::_addToDest(const constvect inv, vect outv) const
{
  Eigen::Map<const Eigen::VectorXd> inm(inv.data(), inv.size());
  Eigen::Map<Eigen::VectorXd> outm(outv.data(), outv.size());
  outm += _eigenMatrix * inm;
  return 0;
}

void MatrixSparse::setDiagonal(const constvect tab)
{
  Eigen::Map<const Eigen::VectorXd> tabm(tab.data(),tab.size());
  setDiagonal(tabm);
}

void MatrixSparse::setDiagonal(const Eigen::Map<const Eigen::VectorXd>& tab)
{
  _eigenMatrix = tab.asDiagonal();
}


///////////////Not exported //////////

Eigen::SparseMatrix<double> AtMA(const Eigen::SparseMatrix<double>& A,
                                 const Eigen::SparseMatrix<double>& M)
{
  return A.transpose() * M * A;
}

