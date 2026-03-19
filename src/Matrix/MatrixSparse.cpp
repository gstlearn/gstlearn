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
#include "Basic/AException.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Law.hpp"
#include "Basic/OptCustom.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/WarningMacro.hpp"
#include "Matrix/AMatrix.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Matrix/NF_Triplet.hpp"

#include <Eigen/src/Core/Map.h>
#include <Eigen/src/Core/Matrix.h>
#include <iostream>

DISABLE_WARNING_PUSH
DISABLE_WARNING_COND_EXPR_CONSTANT
DISABLE_WARNING_UNUSED_BUT_SET_VARIABLE
#include <Eigen/SparseCholesky>
#include <omp.h>
DISABLE_WARNING_POP

/**
 * This variable switches ON/OFF the ability to use Eigen library for Algebra
 */

namespace gstlrn
{

MatrixSparse::MatrixSparse(Id nrow, Id ncol, Id ncolmax)
  : AMatrix(nrow, ncol)
  , _eigenMatrix()
  , _nColMax(ncolmax)
{
  _allocate();
}

MatrixSparse::MatrixSparse(const MatrixSparse& m)
  : AMatrix(m)
  // , ALinearOp(m)
  , _eigenMatrix(m._eigenMatrix)
{
}

MatrixSparse& MatrixSparse::operator=(const MatrixSparse& m)
{
  if (this != &m)
  {
    AMatrix::operator=(m);
    // ALinearOp::operator=(m);
    if (!m.empty())
    {
      _eigenMatrix = m._eigenMatrix;
    }
  }
  return *this;
}

MatrixSparse::~MatrixSparse()
{
  _deallocate();
}

MatrixSparse* MatrixSparse::createFromMatrix(const MatrixSparse& mat)
{
  return new MatrixSparse(mat);
}

void MatrixSparse::resetFromValue(Id nrows, Id ncols, double value)
{
  DECLARE_UNUSED(nrows);
  DECLARE_UNUSED(ncols);
  DECLARE_UNUSED(value);
  _forbiddenForSparse("resetFromValue");
}

void MatrixSparse::resetFromArray(Id nrows, Id ncols, const double* tab, bool byCol)
{
  DECLARE_UNUSED(nrows);
  DECLARE_UNUSED(ncols);
  DECLARE_UNUSED(tab);
  DECLARE_UNUSED(byCol);
  _forbiddenForSparse("resetFromArray");
}

void MatrixSparse::resetFromVD(Id nrows, Id ncols, const VectorDouble& tab, bool byCol)
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
  eigenMat() = NF_T.buildEigenFromTriplet();
  _setNRows(eigenMat().rows());
  _setNCols(eigenMat().cols());
}

/**
 * @brief Fill a sparse matrix with random values, using a given seed and a percentage of zero values
 *
 * @param zeroPercent Percentage of zero values (between 0 and 1)
 * @param seed Random seed for reproducibility
 *
 * @remarks The method also ensures that the last element of the matrix (at position [nrow-1, ncol-1])
 * is non-zero
 */
void MatrixSparse::fillRandom(double zeroPercent, Id seed)
{
  law_set_random_seed(seed);

  auto nrow = getNRows();
  auto ncol = getNCols();
  NF_Triplet NF_T;
  bool lastElementNonZero = false;
  for (Id irow = 0; irow < nrow; irow++)
    for (Id icol = 0; icol < ncol; icol++)
    {
      if (law_uniform(0., 1.) < zeroPercent) continue;
      NF_T.add(irow, icol, law_gaussian());
      if (irow == nrow - 1 && icol == ncol - 1) lastElementNonZero = true;
    }
  if (!lastElementNonZero)
  {
    NF_T.add(nrow - 1, ncol - 1, law_gaussian());
  }
  resetFromTriplet(NF_T);
}

void MatrixSparse::_transposeInPlace()
{
  EigenSparseMatrix temp;
  temp = eigenMat().transpose();
  eigenMat().swap(temp);
}

MatrixSparse* MatrixSparse::transpose() const
{
  auto* mat       = dynamic_cast<MatrixSparse*>(clone());
  mat->eigenMat() = eigenMat().transpose();
  return mat;
}

/**
 * Fill a column of an already existing Sparse matrix, using 'tab' as entry
 * The input 'tab' corresponds to the whole column contents
 * @param icol Column rank
 * @param tab  Vector containing the information (Dimension: nrows)
 */
void MatrixSparse::setColumn(Id icol, const VectorDouble& tab)
{
  auto nrows = getNRows();
  if (getFlagMatrixCheck())
  {
    if (!_isColumnValid(icol)) return;
    if (!_isRowConsistent(tab)) return;
  }
  for (Id irow = 0; irow < nrows; irow++)
    eigenMat().coeffRef(irow, icol) = tab[irow];
}

void MatrixSparse::setColumnToConstant(Id icol, double value)
{
  auto nrows = getNRows();
  if (getFlagMatrixCheck() && !_isColumnValid(icol)) return;
  for (Id irow = 0; irow < nrows; irow++)
    eigenMat().coeffRef(irow, icol) = value;
}

/**
 * Fill a row of an already existing Sparse matrix, using 'tab' as entry
 * The input 'tab' corresponds to the whole row contents
 * @param irow Row rank
 * @param tab  Vector containing the information (Dimension: ncols)
 *
 * @warning: This method only copies the values at the non-zero existing entries
 */
void MatrixSparse::setRow(Id irow, const VectorDouble& tab)
{
  auto ncols = getNCols();
  if (getFlagMatrixCheck())
  {
    if (!_isRowValid(irow)) return;
    if (!_isColumnConsistent(tab)) return;
  }
  for (Id icol = 0; icol < ncols; icol++)
    eigenMat().coeffRef(irow, icol) = tab[icol];
}

void MatrixSparse::setRowToConstant(Id irow, double value)
{
  auto ncols = getNCols();
  if (getFlagMatrixCheck() && !_isRowValid(irow)) return;
  for (Id icol = 0; icol < ncols; icol++)
    eigenMat().coeffRef(irow, icol) = value;
}

void MatrixSparse::setDiagonal(const VectorDouble& tab)
{
  if (getFlagMatrixCheck() && !_isColumnConsistent(tab)) return;
  Eigen::Map<const Eigen::VectorXd> vecm(tab.data(), tab.size());
  eigenMat() = vecm.asDiagonal();
}

void MatrixSparse::setDiagonalToConstant(double value)
{
  if (!isSquare())
    my_throw("This function is only valid for Square matrices");

  VectorDouble vec(getNRows(), value);
  Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), vec.size());
  eigenMat() = vecm.asDiagonal();
}

/*! Gets the value for rank 'rank' */
double MatrixSparse::_getValueByRank(Id rank) const
{
  DECLARE_UNUSED(rank);
  _forbiddenForSparse("_getValueByRank");
  return TEST;
}

double& MatrixSparse::_getValueRef(Id irow, Id icol)
{
  DECLARE_UNUSED(irow);
  DECLARE_UNUSED(icol);
  _forbiddenForSparse("_getValueRef");
  return _eigenMatrix.coeffRef(irow, icol);
}

void MatrixSparse::_setValueByRank(Id rank, double value)
{
  DECLARE_UNUSED(rank);
  DECLARE_UNUSED(value);
  _forbiddenForSparse("_setValueByRank");
}

void MatrixSparse::setValue(Id irow, Id icol, double value)
{
  if (getFlagMatrixCheck() && !_isIndexValid(irow, icol)) return;
  eigenMat().coeffRef(irow, icol) = value;
}

void MatrixSparse::updValue(Id irow,
                            Id icol,
                            const EOperator& oper,
                            double value)
{
  if (getFlagMatrixCheck() && !_isIndexValid(irow, icol)) return;
  double newval                   = modifyOperator(oper, eigenMat().coeff(irow, icol), value);
  eigenMat().coeffRef(irow, icol) = newval;
}

Id MatrixSparse::_getMatrixPhysicalSize() const
{
  return eigenMat().nonZeros();
}

/**
 * @param value Constant value used for filling 'this'
 */
void MatrixSparse::fill(double value)
{
  auto nrow = getNRows();
  auto ncol = getNCols();
  NF_Triplet NF_T;
  for (Id irow = 0; irow < nrow; irow++)
    for (Id icol = 0; icol < ncol; icol++)
      NF_T.add(irow, icol, value);

  resetFromTriplet(NF_T);
}

/*! Multiply a Matrix row-wise */
void MatrixSparse::multiplyRow(const VectorDouble& vec)
{
  if (getFlagMatrixCheck() && !_isRowConsistent(vec)) return;
  for (Id k = 0; k < eigenMat().outerSize(); ++k)
    for (EigenSparseMatrix::InnerIterator it(eigenMat(), k); it; ++it)
      it.valueRef() *= vec[it.row()];
}

/*! Multiply a Matrix column-wise */
void MatrixSparse::multiplyColumn(const VectorDouble& vec)
{
  if (getFlagMatrixCheck() && !_isColumnConsistent(vec)) return;
  for (Id k = 0; k < eigenMat().outerSize(); ++k)
    for (EigenSparseMatrix::InnerIterator it(eigenMat(), k); it; ++it)
      it.valueRef() *= vec[it.col()];
}

/*! Divide a Matrix row-wise */
void MatrixSparse::divideRow(const VectorDouble& vec)
{
  if (getFlagMatrixCheck() && !_isRowConsistent(vec)) return;
  for (Id k = 0; k < eigenMat().outerSize(); ++k)
    for (EigenSparseMatrix::InnerIterator it(eigenMat(), k); it; ++it)
      it.valueRef() /= vec[it.row()];
}

/*! Divide a Matrix column-wise */
void MatrixSparse::divideColumn(const VectorDouble& vec)
{
  if (getFlagMatrixCheck() && !_isColumnConsistent(vec)) return;
  for (Id k = 0; k < eigenMat().outerSize(); ++k)
    for (EigenSparseMatrix::InnerIterator it(eigenMat(), k); it; ++it)
      it.valueRef() /= vec[it.col()];
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
  if (byCol)
  {
    Eigen::Map<const Eigen::MatrixXd> temp(values, getNRows(), getNCols());
    eigenMat() = temp.sparseView(1., EPSILON10);
  }
  else
  {
    Eigen::Map<const Eigen::MatrixXd> temp(values, getNCols(), getNRows());
    eigenMat() = temp.transpose().sparseView(1., EPSILON10);
  }
}
#endif

MatrixSparse* MatrixSparse::createFromCopy(const MatrixSparse* mat)
{
  return new MatrixSparse(*mat);
}

MatrixSparse* MatrixSparse::create(Id nrow, Id ncol)
{
  return new MatrixSparse(nrow, ncol);
}

MatrixSparse* MatrixSparse::createFromTriplet(const NF_Triplet& NF_T,
                                              Id nrow,
                                              Id ncol,
                                              Id nrowmax)
{
  // If 'nrow' a  nd 'ncol' are not defined, derive them from NF_T
  if (nrow <= 0 || ncol <= 0)
  {
    nrow = NF_T.getNRows() + 1;
    ncol = NF_T.getNCols() + 1;
  }
  auto* mat = new MatrixSparse(nrow, ncol, nrowmax);
  mat->resetFromTriplet(NF_T);

  return mat;
}

MatrixSparse* MatrixSparse::createFillRandom(Id nrow,
                                             Id ncol,
                                             double zeroPercent,
                                             Id seed)
{
  auto* mat = new MatrixSparse(nrow, ncol);
  mat->fillRandom(zeroPercent, seed);
  return mat;
}

MatrixSparse* MatrixSparse::Identity(Id nrow, double value)
{
  auto* mat = new MatrixSparse(nrow, nrow);
  for (Id i = 0; i < nrow; i++)
  {
    mat->eigenMat().coeffRef(i, i) += value;
  }
  return mat;
}

MatrixDense* MatrixSparse::createFromSparse(const MatrixSparse& mat)
{
  auto nrow = mat.getNRows();
  auto ncol = mat.getNCols();
  auto* res = new MatrixDense(nrow, ncol);
  for (Id irow = 0; irow < nrow; irow++)
    for (Id icol = 0; icol < ncol; icol++)
      res->setValue(irow, icol, mat.getValue(irow, icol));
  return res;
}

MatrixSparse* MatrixSparse::addMatMat(const MatrixSparse* x,
                                      const MatrixSparse* y,
                                      double cx,
                                      double cy)
{
  auto* mat       = new MatrixSparse(x->getNRows(), x->getNCols(), -1);
  mat->eigenMat() = cx * x->eigenMat() + cy * y->eigenMat();
  return mat;
}

MatrixSparse* MatrixSparse::diagVec(const VectorDouble& vec)
{
  Id size   = static_cast<Id>(vec.size());
  auto* mat = new MatrixSparse(size, size);

  mat->setDiagonal(vec);
  return mat;
}

MatrixSparse* MatrixSparse::diagConstant(Id number, double value)
{
  auto* mat = new MatrixSparse(number, number);
  mat->setDiagonalToConstant(value);
  return mat;
}

/**
 * Construct a sparse matrix with the diagonal of 'A', where each element is transformed
 * @param A    Input sparse matrix
 * @param oper_choice: Operation on the diagonal term (see Utilities::operate_XXX)
 * @return
 */
MatrixSparse* MatrixSparse::diagMat(MatrixSparse* A, Id oper_choice)
{
  if (!A->isSquare())
  {
    messerr("This method requires the matrix 'A' to be square");
    return nullptr;
  }

  VectorDouble diag = A->getDiagonal();
  VectorHelper::transformVD(diag, oper_choice);
  return MatrixSparse::diagVec(diag);
}

bool MatrixSparse::_isElementPresent(Id irow, Id icol) const
{
  for (EigenSparseMatrix::InnerIterator it(eigenMat(), icol); it; ++it)
  {
    if (it.row() == irow) return true;
  }
  return false;
}

void MatrixSparse::addValue(Id row, Id col, double value)
{
  if (ABS(value) <= EPSILON10) return;
  eigenMat().coeffRef(row, col) += value;
}

double MatrixSparse::getValue(Id row, Id col) const
{
  if (getFlagMatrixCheck() && !_isIndexValid(row, col)) return TEST;
  return eigenMat().coeff(row, col);
}

double MatrixSparse::L1Norm() const
{
  return (Eigen::RowVectorXd::Ones(eigenMat().rows()) * eigenMat().cwiseAbs()).maxCoeff();
}

void MatrixSparse::getStats(Id* nrows, Id* ncols, Id* count, double* percent) const
{
  *nrows   = getNRows();
  *ncols   = getNCols();
  *count   = eigenMat().nonZeros();
  *percent = 0.;
  if ((*nrows) > 0 && (*ncols) > 0)
    (*percent) = ((100. * static_cast<double>(*count)) / (static_cast<double>(*nrows) * static_cast<double>(*ncols)));
}

VectorDouble MatrixSparse::extractDiag(Id oper_choice) const
{
  VectorDouble diag(std::min(getNCols(), getNRows()));
  Eigen::Map<Eigen::VectorXd> ym(diag.data(), diag.size());
  ym = eigenMat().diagonal();
  VH::transformVD(diag, oper_choice);
  return diag;
}

Id MatrixSparse::addVecInPlaceEigen(const Eigen::Map<const Eigen::VectorXd>& xm,
                                    Eigen::Map<Eigen::VectorXd>& ym) const
{
  ym = eigenMat() * xm + ym;
  return 0;
}

Id MatrixSparse::addVecInPlace(const constvect xm, vect ym) const
{
  Eigen::Map<const Eigen::VectorXd> xmm(xm.data(), xm.size());
  Eigen::Map<Eigen::VectorXd> ymm(ym.data(), ym.size());
  ymm = eigenMat() * xmm + ymm;
  return 0;
}

Id MatrixSparse::addVecInPlaceVD(const VectorDouble& x, VectorDouble& y) const
{
  Eigen::Map<const Eigen::VectorXd> xm(x.data(), x.size());
  Eigen::Map<Eigen::VectorXd> ym(y.data(), y.size());
  ym = eigenMat() * xm + ym;
  return 0;
}

void MatrixSparse::setConstant(double value)
{
  for (Id k = 0; k < eigenMat().outerSize(); ++k)
    for (EigenSparseMatrix::InnerIterator it(eigenMat(), k); it; ++it)
      it.valueRef() = value;
}

Id MatrixSparse::scaleByDiag()
{
  VectorDouble diag = extractDiag(-1);
  Eigen::Map<const Eigen::VectorXd> ym(diag.data(), diag.size());
  eigenMat() = ym.asDiagonal() * eigenMat();
  return 0;
}

void MatrixSparse::_productGeneral(MatrixSparse& res,
                                   const MatrixSparse& other1,
                                   const MatrixSparse& other2,
                                   bool transpose1,
                                   bool transpose2)
{
  if (transpose1)
  {
    if (transpose2)
      res.eigenMat() = other1.eigenMat().transpose() * other2.eigenMat().transpose();
    else
      res.eigenMat() = other1.eigenMat().transpose() * other2.eigenMat();
  }
  else
  {
    if (transpose2)
      res.eigenMat() = other1.eigenMat() * other2.eigenMat().transpose();
    else
      res.eigenMat() = other1.eigenMat() * other2.eigenMat();
  }
}

/**
 * @brief Add the Product of a matrix by a vector
 *
 * @param res Resulting vector
 * @param other Input matrix
 * @param vec  Input vector
 * @param transpose Should the matrix 'other' be transposed
 * @param flagInvert Product 'other' * 'vec' (F) or 'vec' * 'other' (T)
 */
void MatrixSparse::_prodVecAddGeneral(VectorDouble& res,
                                      const MatrixSparse& other,
                                      const VectorDouble& vec,
                                      bool transpose,
                                      bool flagInvert)
{
  Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), vec.size());
  Eigen::Map<Eigen::VectorXd> resm(res.data(), res.size());

  if (flagInvert)
  {
    // In order to obtain the correct result as a vector, the product is coded in a transposed manner
    if (transpose)
    {
      resm.noalias() += other.eigenMat() * vecm;
    }
    else
    {
      resm.noalias() += other.eigenMat().transpose() * vecm;
    }
  }
  else
  {
    if (transpose)
    {
      resm.noalias() += other.eigenMat().transpose() * vecm;
    }
    else
    {
      resm.noalias() += other.eigenMat() * vecm;
    }
  }
}

/**
 * @brief Add the Product of a matrix by a vector
 *
 * @param res Resulting vector
 * @param other Input matrix
 * @param vec  Input vector
 * @param transpose Should the matrix 'other' be transposed
 * @param flagInvert Product 'other' * 'vec' (F) or 'vec' * 'other' (T)
 */
void MatrixSparse::_prodVecAddGeneral(vect& res,
                                      const MatrixSparse& other,
                                      const constvect& vec,
                                      bool transpose,
                                      bool flagInvert)
{
  Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), vec.size());
  Eigen::Map<Eigen::VectorXd> resm(res.data(), res.size());

  if (flagInvert)
  {
    // In order to obtain the correct result as a vector, the product is coded in a transposed manner
    if (transpose)
    {
      resm.noalias() += other.eigenMat() * vecm;
    }
    else
    {
      resm.noalias() += other.eigenMat().transpose() * vecm;
    }
  }
  else
  {
    if (transpose)
    {
      resm.noalias() += other.eigenMat().transpose() * vecm;
    }
    else
    {
      resm.noalias() += other.eigenMat() * vecm;
    }
  }
}

/**
 * Store the product of 'x' by 'y' in this
 * @param x First Matrix
 * @param y Second matrix
 * @param transposeX True if First matrix is transposed
 * @param transposeY True if Second matrix is transposed
 */
void MatrixSparse::prodMatMatInPlace(const AMatrix* x,
                                     const AMatrix* y,
                                     bool transposeX,
                                     bool transposeY)
{
  if (getFlagMatrixCheck() &&
      !_isMatrixCompatible("MatrixSparse::prodMatMatInPlace",
                           x, 0, transposeX,
                           y, 0, transposeY)) return;

  const auto* xm = dynamic_cast<const MatrixSparse*>(x);
  const auto* ym = dynamic_cast<const MatrixSparse*>(y);
  if (xm == nullptr || ym == nullptr)
  {
    AMatrix::prodMatMatInPlace(x, y, transposeX, transposeY);
  }
  else
  {
    MatrixSparse temp;
    AMatrix::productInPlace(temp, *xm, *ym, transposeX, transposeY);

    if (transposeX)
    {
      if (transposeY)
        eigenMat() = xm->eigenMat().transpose() * ym->eigenMat().transpose();
      else
        eigenMat() = xm->eigenMat().transpose() * ym->eigenMat();
    }
    else
    {
      if (transposeY)
        eigenMat() = xm->eigenMat() * ym->eigenMat().transpose();
      else
        eigenMat() = xm->eigenMat() * ym->eigenMat();
    }
    (void)areIdentical(temp, *this);
  }
}

MatrixSparse* prodNormMatMat(const MatrixSparse* a,
                             const MatrixSparse* m,
                             bool transpose)
{
  Id nrow   = (transpose) ? a->getNCols() : a->getNRows();
  Id ncol   = (transpose) ? a->getNRows() : a->getNCols();
  auto* mat = new MatrixSparse(nrow, ncol);
  mat->prodNormMatMatInPlace(a, m, transpose);
  return mat;
}

MatrixSparse* prodNormMatVec(const MatrixSparse* a, const VectorDouble& vec, bool transpose)
{
  Id nsym   = (transpose) ? a->getNCols() : a->getNRows();
  auto* mat = new MatrixSparse(nsym, nsym);
  mat->prodNormMatVecInPlace(a, vec, transpose);
  return mat;
}

MatrixSparse* prodNormMat(const MatrixSparse* a, bool transpose)
{
  Id nsym   = (transpose) ? a->getNCols() : a->getNRows();
  auto* mat = new MatrixSparse(nsym, nsym);
  mat->prodNormMatInPlace(a, transpose);
  return mat;
}

MatrixSparse* prodNormDiagVec(const MatrixSparse* a,
                              const VectorDouble& vec,
                              Id oper_choice)
{
  Id nrow   = a->getNRows();
  Id ncol   = a->getNCols();
  auto* mat = new MatrixSparse(nrow, ncol, -1);

  // Perform the transformation of the input vector
  VectorDouble vecp = vec;
  VH::transformVD(vecp, oper_choice);

  Eigen::Map<const Eigen::VectorXd> vecm(vecp.data(), vecp.size());
  auto diag       = vecm.asDiagonal();
  mat->eigenMat() = (diag * a->eigenMat() * diag);
  return mat;
}

/**
 * Perform: 'this' = diag('vec') %*% 'A' %*% diag('vec')
 * @param vec  Input Vector
 * @param oper_choice Type of transformation
 */
void MatrixSparse::prodNormDiagVecInPlace(const VectorDouble& vec, Id oper_choice)
{
  if (!isSquare())
  {
    messerr("This method is limited to square matrices");
    return;
  }
  if (getNRows() != static_cast<Id>(vec.size()))
  {
    messerr("Matrix dimension (%d) does not match vector dimension (%d)",
            getNRows(), static_cast<Id>(vec.size()));
    return;
  }

  // Perform the transformation of the input vector
  VectorDouble vecp = vec;
  VH::transformVD(vecp, oper_choice);

  Eigen::Map<const Eigen::VectorXd> vecm(vecp.data(), vecp.size());
  auto diag  = vecm.asDiagonal();
  eigenMat() = diag * eigenMat() * diag;
}

void MatrixSparse::prodNormMatVecInPlace(const AMatrix* a,
                                         const VectorDouble& vec,
                                         bool transpose)
{
  if (getFlagMatrixCheck())
  {
    // Note that 'vec' is not tested for compatibility:
    // it is a vector but used as a diagonal of a square matrix
    if (!_isMatrixCompatible("MatrixSparse::prodNormMatVecInPlace",
                             a, 0, transpose,
                             a, 0, !transpose)) return;
  }

  const auto* am = dynamic_cast<const MatrixSparse*>(a);
  if (am == nullptr)
  {
    AMatrix::prodNormMatVecInPlace(a, vec, transpose);
  }
  else
  {
    if (transpose)
    {
      Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), vec.size());
      eigenMat() = am->eigenMat().transpose() * vecm.asDiagonal() * am->eigenMat();
    }
    else
    {
      Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), vec.size());
      eigenMat() = am->eigenMat() * vecm.asDiagonal() * am->eigenMat().transpose();
    }
  }

  MatrixSparse temp;
  AMatrix::prodnormInPlace(temp, *am, vec, transpose);
  (void)areIdentical(temp, *this);
}

void MatrixSparse::prodNormMatInPlace(const AMatrix* a, bool transpose)
{
  if (getFlagMatrixCheck() &&
      !_isMatrixCompatible("MatrixSparse::prodNormMatInPlace",
                           a, 0, transpose,
                           a, 0, !transpose)) return;

  const auto* am = dynamic_cast<const MatrixSparse*>(a);
  if (am == nullptr)
  {
    AMatrix::prodNormMatInPlace(a, transpose);
  }
  else
  {
    if (transpose)
    {
      eigenMat() = am->eigenMat().transpose() * am->eigenMat();
    }
    else
    {
      eigenMat() = am->eigenMat() * am->eigenMat().transpose();
    }
  }

  MatrixSparse temp;
  AMatrix::prodnormInPlace(temp, *am, MatrixSparse(), transpose);
  (void)areIdentical(temp, *this);
}

void MatrixSparse::prodNormMatMatInPlace(const AMatrix* a,
                                         const AMatrix* m,
                                         bool transpose)
{
  if (getFlagMatrixCheck() &&
      !_isMatrixCompatible("MatrixSparse::prodNormMatMatInPlace",
                           a, 0, transpose,
                           m, 0, false,
                           a, 0, !transpose)) return;

  const auto* am = dynamic_cast<const MatrixSparse*>(a);
  const auto* mm = dynamic_cast<const MatrixSparse*>(m);
  if (am == nullptr || mm == nullptr)
  {
    AMatrix::prodNormMatMatInPlace(a, m, transpose);
  }
  else
  {
    if (transpose)
    {
      eigenMat() = (am->eigenMat().transpose() * mm->eigenMat()) * am->eigenMat();
    }
    else
    {
      eigenMat() = (am->eigenMat() * mm->eigenMat()) * am->eigenMat().transpose();
    }
  }

  MatrixSparse temp;
  AMatrix::prodnormInPlace(temp, *am, *mm, transpose);
  (void)areIdentical(temp, *this);
}

Id MatrixSparse::_invert()
{
  if (!isSquare())
    my_throw("Invert method is restricted to Square matrices");
  auto n = getNCols();
  Eigen::SimplicialLLT<EigenSparseMatrix> solver;
  solver.compute(eigenMat());
  EigenSparseMatrix I(n, n);
  I.setIdentity();
  eigenMat() = solver.solve(I);

  return 0;
}

Id MatrixSparse::_solve(const VectorDouble& b, VectorDouble& x) const
{
  if (!isSquare())
    my_throw("Invert method is limited to Square Matrices");
  if (static_cast<Id>(b.size()) != getNRows() || static_cast<Id>(x.size()) != getNRows())
    my_throw("b' and 'x' should have the same dimension as the Matrix");

  Eigen::SimplicialLLT<EigenSparseMatrix> solver;
  Eigen::Map<const Eigen::VectorXd> bm(b.data(), getNCols());
  Eigen::Map<Eigen::VectorXd> xm(x.data(), getNRows());
  xm = solver.compute(eigenMat()).solve(bm);

  return 0;
}

String MatrixSparse::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  sstr << AMatrix::toString(strfmt) << std::endl;
  return sstr.str();
}

void MatrixSparse::_allocate(Id nrow, Id ncol, Id ncolmax)
{
  eigenMat() = EigenSparseMatrix(nrow, ncol);
  _setNCols(ncol);
  _setNRows(nrow);

  if (ncolmax > 0)
  {
    eigenMat().reserve(Eigen::VectorXi::Constant(nrow, static_cast<I32>(ncolmax)));
  }
  _nColMax = ncolmax;
}

/**
 * This strange function instantiate a sparse matrix with given dimensions
 * filled with zeroes. It should be an empty matrix... But this does not make sense.
 * Therefore it is created by setting a single element at the lower bottom size of
 * the matrix ... filled with a zero.
 */
void MatrixSparse::_allocate()
{
  eigenMat() = EigenSparseMatrix(getNRows(), getNCols());
  {
    eigenMat() = EigenSparseMatrix(getNRows(), getNCols());

    if (_nColMax > 0)
    {
      eigenMat().reserve(Eigen::VectorXi::Constant(getNCols(), static_cast<I32>(_nColMax)));
    }
    const auto nbthread = static_cast<I32>(OptCustom::query("ompthreads", 1));
    Eigen::setNbThreads(nbthread);
  }
}

void MatrixSparse::_deallocate()
{
  eigenMat().data().squeeze();
  {
    eigenMat().data().squeeze();
  }
}

void MatrixSparse::_forbiddenForSparse(const String& func)
{
  messerr("Problem with Function: %s", func.c_str());
  messerr("This function is not available in Sparse Matrix");
}

void MatrixSparse::dumpElements(const String& title, Id ifrom, Id ito)
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
NF_Triplet MatrixSparse::getMatrixToTriplet(Id shiftRow, Id shiftCol) const
{
  return NF_Triplet::createFromEigen(eigenMat(), shiftRow, shiftCol);
}

void MatrixSparse::_clear()
{
  _setNRows(0);
  _setNCols(0);
  _allocate();
}

Id MatrixSparse::_getIndexToRank(Id irow, Id icol) const
{
  DECLARE_UNUSED(irow);
  DECLARE_UNUSED(icol);
  _forbiddenForSparse("_getIndexToRank");
  return ITEST;
}

VectorInt MatrixSparse::getNonZeroCols(Id irow) const
{
  VectorInt cols;
  for (EigenSparseMatrix::InnerIterator it(eigenMat(), irow); it; ++it)
  {
    cols.push_back(it.index());
  }
  return cols;
}

VectorInt MatrixSparse::getNonZeroRows(Id icol) const
{
  if (!_isColumnValid(icol))
  {
    return VectorInt();
  }
  VectorInt rows;

  // On force la création d'une matrice RowMajor (C'est ICI que la copie a lieu)
  // Eigen va réorganiser les données en mémoire.
  Eigen::SparseMatrix<double, Eigen::RowMajor> rowMat = eigenMat();

  // Maintenant, l'itérateur fonctionne car les lignes sont contiguës
  for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(rowMat, icol); it; ++it)
  {
    rows.push_back(it.index());
  }

  return rows;
}

MatrixSparse* createFromAnyMatrix(const AMatrix* matin)
{
  return MatrixSparse::createFromTriplet(matin->getMatrixToTriplet(),
                                         matin->getNRows(),
                                         matin->getNCols(),
                                         -1);
}

Id MatrixSparse::_eigen_findColor(Id imesh,
                                  Id ncolor,
                                  VectorInt& colors,
                                  VectorInt& temp) const
{
  temp.fill(0);

  /* Checks the colors of the connected nodes */

  for (EigenSparseMatrix::InnerIterator it(eigenMat(), imesh); it; ++it)
  {
    if (isZero(it.value())) continue;
    Id irow = it.row();
    if (!isNA(colors[irow])) temp[colors[irow] - 1]++;
  }

  /* Look for a free color */

  for (Id j = 0; j < ncolor; j++)
  {
    if (temp[j] == 0) return (j + 1);
  }
  return (-1);
}

VectorInt MatrixSparse::colorCoding() const
{
  Id next_col = 0;
  Id ncol     = 0;
  auto nmesh  = getNCols();

  /* Core allocation */

  VectorInt colors(nmesh, ITEST);
  VectorInt temp(nmesh);

  /* Loop on the nodes of the mesh */

  for (Id imesh = 0; imesh < nmesh; imesh++)
  {
    next_col = _eigen_findColor(imesh, ncol, colors, temp);
    next_col = _eigen_findColor(imesh, ncol, colors, temp);

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
  Id shiftRow   = (flagShiftRow) ? A1->getNRows() : 0;
  Id shiftCol   = (flagShiftCol) ? A1->getNCols() : 0;
  NF_Triplet T1 = A1->getMatrixToTriplet();
  NF_Triplet T2 = A2->getMatrixToTriplet(shiftRow, shiftCol);

  // Concatenate the two triplet lists
  T1.appendInPlace(T2);

  // Create the new matrix from the resulting triplet list
  Id nrow = (flagShiftRow) ? A1->getNRows() + A2->getNRows() : MAX(A1->getNRows(), A2->getNRows());
  Id ncol = (flagShiftCol) ? A1->getNCols() + A2->getNCols() : MAX(A1->getNCols(), A2->getNCols());

  A1->resetFromTriplet(T1);
  A1->_setNRows(nrow);
  A1->_setNCols(ncol);
}

MatrixSparse* MatrixSparse::glue(const MatrixSparse* A1,
                                 const MatrixSparse* A2,
                                 bool flagShiftRow,
                                 bool flagShiftCol)
{
  Id shiftRow = (flagShiftRow) ? A1->getNRows() : 0;
  Id shiftCol = (flagShiftCol) ? A1->getNCols() : 0;

  // Create the two triplet lists
  NF_Triplet T1 = A1->getMatrixToTriplet();
  NF_Triplet T2 = A2->getMatrixToTriplet(shiftRow, shiftCol);

  // Concatenate the two triplet lists
  T1.appendInPlace(T2);

  // Create the new matrix from the resulting triplet list
  Id nrow = (flagShiftRow) ? A1->getNRows() + A2->getNRows() : MAX(A1->getNRows(), A2->getNRows());
  Id ncol = (flagShiftCol) ? A1->getNCols() + A2->getNCols() : MAX(A1->getNCols(), A2->getNCols());

  return MatrixSparse::createFromTriplet(T1, nrow, ncol, -1);
}

/* Extract a sparse sub-matrix */
/* 'rank_rows' and 'rank_cols' must have same dimension as C */
/* The arrays 'rank_rows' and 'rank_cols' may be absent */
/* Their value gives the rank of the saved element or -1 */
MatrixSparse* MatrixSparse::extractSubmatrixByRanks(const VectorInt& rank_rows,
                                                    const VectorInt& rank_cols) const
{
  Id old_row, old_col, new_row, new_col;

  NF_Triplet NF_Tin = getMatrixToTriplet();
  NF_Triplet NF_Tout;

  /* Fill the new sparse triplet */

  for (Id i = 0; i < NF_Tin.getNElements(); i++)
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
MatrixSparse* MatrixSparse::extractSubmatrixByColor(const VectorInt& colors,
                                                    Id ref_color,
                                                    bool row_ok,
                                                    bool col_ok)
{
  /* Convert the contents of the sparse matrix into columns */

  NF_Triplet NF_Tin = getMatrixToTriplet();

  /* Initialize the output matrix */

  NF_Triplet NF_Tout;

  /* Core allocation */

  auto n = getNCols();
  VectorInt u_row(n);
  VectorInt u_col(n);

  Id ir = 0;
  for (Id i = 0; i < n; i++)
  {
    u_row[i] = -1;
    if (row_ok && colors[i] != ref_color) continue;
    if (!row_ok && colors[i] == ref_color) continue;
    u_row[i] = ir++;
  }

  Id ic = 0;
  for (Id i = 0; i < n; i++)
  {
    u_col[i] = -1;
    if (col_ok && colors[i] != ref_color) continue;
    if (!col_ok && colors[i] == ref_color) continue;
    u_col[i] = ic++;
  }

  /* Fill the new sparse triplet */

  for (Id i = 0; i < NF_Tin.getNElements(); i++)
  {
    ir = u_row[NF_Tin.getRow(i)];
    ic = u_col[NF_Tin.getCol(i)];
    if (ir < 0 || ic < 0) continue;
    NF_Tout.add(ir, ic, NF_Tin.getValue(i));
  }

  return MatrixSparse::createFromTriplet(NF_Tout, 0, 0, -1);
}

void MatrixSparse::gibbs(Id iech,
                         const VectorDouble& zcur,
                         double* yk,
                         double* sk)
{
  *yk = 0.;
  for (EigenSparseMatrix::InnerIterator it(eigenMat(), iech); it;
       ++it)
  {
    double coeff = it.valueRef();
    if (ABS(coeff) <= 0.) continue;
    Id jech = it.row();
    {
      if (iech == jech)
        *sk = coeff;
      else
        *yk -= coeff * zcur[jech];
    }
  }

  // Returned arguments
  (*yk) /= (*sk);
  (*sk) = sqrt(1. / (*sk));
}

Id MatrixSparse::_addToDest(const constvect inv, vect outv) const
{
  Eigen::Map<const Eigen::VectorXd> inm(inv.data(), inv.size());
  Eigen::Map<Eigen::VectorXd> outm(outv.data(), outv.size());
  outm += eigenMat() * inm;
  return 0;
}

void MatrixSparse::setDiagonal(const constvect tab)
{
  Eigen::Map<const Eigen::VectorXd> tabm(tab.data(), tab.size());
  setDiagonal(tabm);
}

void MatrixSparse::setDiagonal(const Eigen::Map<const Eigen::VectorXd>& tab)
{
  eigenMat() = tab.asDiagonal();
}

/*! Extract a Row */
MatrixSparse* MatrixSparse::getRowAsMatrixSparse(Id irow, double coeff) const
{
  auto ncols = getNCols();
  auto* res  = new MatrixSparse(1, ncols);

  // The input sparse matrix being symmetrical, we benefit from its
  // column-major storage (setting icol = irow)
  Id icol = irow;
  for (EigenSparseMatrix::InnerIterator it(eigenMat(), icol); it; ++it)
    res->eigenMat().coeffRef(0, it.row()) = coeff * it.value();

  return res;
}

/*! Extract a Column */
MatrixSparse* MatrixSparse::getColumnAsMatrixSparse(Id icol, double coeff) const
{
  auto nrows = getNRows();
  auto* res  = new MatrixSparse(nrows, 1);

  for (EigenSparseMatrix::InnerIterator it(eigenMat(), icol); it; ++it)
    res->eigenMat().coeffRef(it.row(), 0) = coeff * it.value();

  return res;
}

///////////////Not exported //////////

EigenSparseMatrix AtMA(const EigenSparseMatrix& A,
                       const EigenSparseMatrix& M)
{
  return A.transpose() * M * A;
}

Id MatrixSparse::forwardLU(const VectorDouble& b, VectorDouble& x, bool flagLower) const
{
  Eigen::Map<const Eigen::VectorXd> bm(b.data(), b.size());
  Eigen::Map<Eigen::VectorXd> xm(x.data(), x.size());

  if (!flagLower)
  {
    const EigenSparseMatrix& Lx = eigenMat().transpose();
    xm                          = Lx.triangularView<Eigen::Upper>().solve(bm);
  }
  else
  {
    const EigenSparseMatrix& Lx = eigenMat();
    xm                          = Lx.triangularView<Eigen::Lower>().solve(bm);
  }
  return 0;
}

void MatrixSparse::forceDimension(Id maxRows, Id maxCols)
{
  // Redimensionner les matrices si nécessaire
  if (eigenMat().rows() < maxRows || eigenMat().cols() < maxCols)
  {
    eigenMat().conservativeResize(maxRows, maxCols);
    eigenMat().insert(maxRows - 1, maxCols - 1) = 0.0; // Élément fictif
    _setNRows(maxRows);
    _setNCols(maxCols);
  }
}

/*! New methodes for operator overloading */
MatrixSparse& MatrixSparse::addCst(double a)
{
  _addGeneral(*this, *this, a);
  return *this;
}

MatrixSparse& MatrixSparse::prodCst(double a)
{
  _prodGeneral(*this, *this, a);
  return *this;
}

void MatrixSparse::_addGeneral(MatrixSparse& res,
                               const MatrixSparse& other,
                               double cst)
{
  if (isZero(cst))
  {
    res = other;
    return;
  }

  MatrixSparse temp(other.getNRows(), other.getNCols());
  temp.fill(cst);

  auto result    = other.eigenMat() + temp.eigenMat();
  res.eigenMat() = result;
}

/**
 * @brief Add matrices 'other1' and 'other2' and store the result in 'res'
 * This code is alais-safe, as it does not modify 'other1' and 'other2' and only writes to 'res'.
 *
 * @param res Resulting matrix
 * @param other1 First matrix
 * @param other2 Second matrix
 */
void MatrixSparse::_addGeneral(MatrixSparse& res,
                               const MatrixSparse& other1,
                               const MatrixSparse& other2)
{
  res.eigenMat() = other1.eigenMat() + other2.eigenMat();
}

void MatrixSparse::_prodGeneral(MatrixSparse& res,
                                const MatrixSparse& other,
                                double cst)
{
  if (isOne(cst))
  {
    res = other;
    return;
  }

  res.eigenMat() = other.eigenMat() * cst;
}

/**
 * @brief Make the elementwise product of matrices 'other1' and 'other2' and store the result in 'res'
 * This code is alais-safe, as it does not modify 'other1' and 'other2' and only writes to 'res'.
 *
 * @param res Resulting matrix
 * @param other1 First matrix
 * @param other2 Second matrix
 */
void MatrixSparse::_prodGeneral(MatrixSparse& res,
                                const MatrixSparse& other1,
                                const MatrixSparse& other2)
{
  res.eigenMat() = other1.eigenMat().cwiseProduct(other2.eigenMat());
}

void MatrixSparse::_prodnormGeneral(MatrixSparse& res,
                                    const MatrixSparse& a,
                                    const MatrixSparse& m,
                                    bool transpose)
{
  if (m.empty())
  {
    if (transpose)
    {
      res.eigenMat() = a.eigenMat().transpose() * a.eigenMat();
    }
    else
    {
      res.eigenMat() = a.eigenMat() * a.eigenMat().transpose();
    }
  }
  else
  {
    if (transpose)
    {
      res.eigenMat() = a.eigenMat().transpose() * m.eigenMat() * a.eigenMat();
    }
    else
    {
      res.eigenMat() = a.eigenMat() * m.eigenMat() * a.eigenMat().transpose();
    }
  }
}

void MatrixSparse::_prodnormGeneral(MatrixSparse& res,
                                    const MatrixSparse& a,
                                    const VectorDouble& vec,
                                    bool transpose)
{
  Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), vec.size());
  if (transpose)
  {
    res.eigenMat() = a.eigenMat().transpose() * vecm * a.eigenMat();
  }
  else
  {
    res.eigenMat() = a.eigenMat() * vecm * a.eigenMat().transpose();
  }
}

bool MatrixSparse::_areIdenticalGeneral(const MatrixSparse& a, const MatrixSparse& b, bool verbose)
{
  if (a.getNRows() != b.getNRows() || a.getNCols() != b.getNCols())
    return false;

  Eigen::SparseMatrix<double> diff = a.eigenMat() - b.eigenMat();
  bool flag                        = diff.nonZeros() == 0;
  if (!flag && verbose)
    messerr("Matrices are not identical ");
  return flag;
}

}; // namespace gstlrn
