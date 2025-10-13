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
#include "Matrix/MatrixSquare.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Matrix/AMatrix.hpp"
#include "Matrix/MatrixFactory.hpp"

#include <cmath>

namespace gstlrn
{
MatrixSquare::MatrixSquare(Id nrow)
  : MatrixDense(nrow, nrow)
{
}

MatrixSquare::MatrixSquare(const MatrixSquare& r)
  : MatrixDense(r)
{
}

MatrixSquare::MatrixSquare(const AMatrix& m)
  : MatrixDense(m)
{
  if (!m.isSquare())
  {
    messerr("The input matrix should be Square");
    _clear();
    return;
  }
}

MatrixSquare& MatrixSquare::operator=(const MatrixSquare& r)
{
  if (this != &r)
  {
    MatrixDense::operator=(r);
  }
  return *this;
}

MatrixSquare::~MatrixSquare()
{
}

void MatrixSquare::_setNSize(Id nval)
{
  _setNRows(nval);
  _setNCols(nval);
}

void MatrixSquare::resetFromVVD(const VectorVectorDouble& tab, bool byCol)
{
  if (tab.empty()) return;
  Id n1 = static_cast<Id>(tab.size());
  Id n2 = static_cast<Id>(tab[0].size());
  if (n1 != n2)
  {
    messerr("The Matrix should be square");
    messerr("Loading is not performed");
    return;
  }
  AMatrix::resetFromVVD(tab, byCol);
}

double MatrixSquare::trace() const
{
  double res = 0.;
  for (Id irow = 0; irow < getNSize(); irow++)
    res += getValue(irow, irow);
  return res;
}
/**
 * Perform the product: this = t(R1) %*% X %*% R2 + t(R2) %*% X %*% R1
 * @param x: Square matrix
 * @param r1: Left Hand Matrix
 * @param r2: Right Hand Matrix
 * \remarks The number of rows of Y must be equal to the dimension of X
 * \remarks The output matrix is square with dimension equal to the number of columns of Y
 */
void MatrixSquare::innerMatrix(const MatrixSquare& x,
                               const AMatrix& r1,
                               const AMatrix& r2)
{
  auto n = x.getNSize();
  if (n != r1.getNRows())
  {
    my_throw("Incompatible matrix dimensions");
  }
  if (n != r2.getNRows())
  {
    my_throw("Incompatible matrix dimensions");
  }

  for (Id irow = 0; irow < n; irow++)
  {
    for (Id icol = 0; icol < n; icol++)
    {
      double value = 0.;
      for (Id k = 0; k < n; k++)
      {
        for (Id l = 0; l < n; l++)
        {
          value += r1.getValue(k, irow) * x.getValue(k, l) * r2.getValue(l, icol);
          value += r2.getValue(k, irow) * x.getValue(k, l) * r1.getValue(l, icol);
        }
      }
      setValue(irow, icol, value);
    }
  }
}

/*! Multiply the diagonal by a vector */
void MatrixSquare::prodDiagByVector(const VectorDouble& diag)
{
  if (static_cast<Id>(diag.size()) != getNRows())
    my_throw("Argument 'Diag' must match Matrix dimension");

  for (Id i = 0; i < getNRows(); i++)
  {
    setValue(i, i, getValue(i, i) * diag[i]);
  }
}

/*! Divide the diagonal by a vector */
void MatrixSquare::divideDiagByVector(const VectorDouble& diag)
{
  if (static_cast<Id>(diag.size()) != getNRows())
    my_throw("Argument 'Diag' must match Matrix dimension");

  for (Id i = 0; i < getNRows(); i++)
  {
    if (isZero(diag[i]))
      my_throw("Argument 'Diag' may not have too small values");
    setValue(i, i, getValue(i, i) / diag[i]);
  }
}

double MatrixSquare::determinant(void) const
{
  auto neq = getNRows();

  /* Dispatch according to the matrix dimension */

  double deter = 0.;
  switch (neq)
  {
    case 1:
      deter = getValue(0, 0);
      break;

    case 2:
      deter = (getValue(0, 0) * getValue(1, 1) - getValue(1, 0) * getValue(0, 1));
      break;

    case 3:
      deter = (getValue(0, 0) * getValue(1, 1) * getValue(2, 2) + getValue(0, 1) * getValue(1, 2) * getValue(2, 0) + getValue(1, 0) * getValue(2, 1) * getValue(0, 2) - getValue(2, 0) * getValue(1, 1) * getValue(0, 2) - getValue(1, 0) * getValue(0, 1) * getValue(2, 2) - getValue(2, 1) * getValue(1, 2) * getValue(0, 0));
      break;

    default:
      Id neqm1 = neq - 1;
      MatrixSquare c(neqm1); // TODO Merge MatrixSquare and MatrixSquare into MatrixSquare

      for (Id j1 = 0; j1 < neq; j1++)
      {
        for (Id i = 1; i < neq; i++)
        {
          Id j2 = 0;
          for (Id j = 0; j < neq; j++)
          {
            if (j == j1) continue;
            c.setValue(i - 1, j2, getValue(i, j));
            j2++;
          }
        }
        deter += pow(-1.0, j1 + 2.0) * getValue(0, j1) * c.determinant();
      }
      break;
  }

  return (deter);
}

double MatrixSquare::normVec(const VectorDouble& vec)
{
  if (getNRows() != static_cast<Id>(vec.size()))
  {
    messerr("Wrong dimension of 'vec' argument");
    return TEST;
  }
  double value = 0.;
  for (Id irow = 0; irow < getNRows(); irow++)
    for (Id icol = 0; icol < getNCols(); icol++)
      value += vec[irow] * getValue(irow, icol) * vec[icol];

  return value;
}

/*****************************************************************************/
/*!
 **  Performs the 'this' %*% diag(c) where c is a vector
 **
 ** \param[in]  mode  0: c as is; 1: sqrt(c); 2: 1/c; 3: 1/sqrt(c)
 ** \param[in]  c     vector
 **
 *****************************************************************************/
void MatrixSquare::prodByDiagInPlace(Id mode, const VectorDouble& c)
{
  auto neq = getNRows();
  for (Id i1 = 0; i1 < neq; i1++)
    for (Id i2 = 0; i2 < neq; i2++)
    {
      double val = c[i2];
      if (mode == 1)
        val = sqrt(val);
      else if (mode == 2)
        val = 1. / val;
      else if (mode == 3)
        val = 1. / sqrt(val);
      setValue(i1, i2, getValue(i1, i2) * val);
    }
}

bool MatrixSquare::_isNumbersValid(Id nrows, Id ncols) const
{
  AMatrix::_isNumbersValid(nrows, ncols);
  if (nrows != ncols)
  {
    messerr("Arguments 'nrows' and 'ncols' should be equal for Square Matrices");
    return false;
  }
  return true;
}
/**
 * Converts a VectorVectorDouble into a Matrix
 * Note: the input argument is stored by row (if coming from [] specification)
 * @param X Input VectorVectorDouble argument
 * @return The returned square matrix
 *
 * @remark: the matrix is transposed implicitly while reading
 */
MatrixSquare* MatrixSquare::createFromVVD(const VectorVectorDouble& X)
{
  Id nrow = static_cast<Id>(X.size());
  Id ncol = static_cast<Id>(X[0].size());
  if (nrow != ncol)
  {
    messerr("The matrix does not seem to be square");
    return nullptr;
  }

  auto* mat = new MatrixSquare(nrow);
  mat->_fillFromVVD(X);
  return mat;
}

MatrixSquare* MatrixSquare::createFromVD(const VectorDouble& X,
                                         Id nrow,
                                         bool byCol,
                                         bool invertColumnOrder)
{
  Id ncol = nrow;
  if (nrow * ncol != static_cast<Id>(X.size()))
  {
    messerr("Inconsistency between arguments 'nrow'(%d) and 'ncol'(%d)", nrow, ncol);
    messerr("and the dimension of the input Vector (%d)", static_cast<Id>(X.size()));
  }
  auto* mat = new MatrixSquare(nrow);

  Id lec = 0;
  if (byCol)
  {
    for (Id irow = 0; irow < nrow; irow++)
      for (Id icol = 0; icol < ncol; icol++)
      {
        Id jcol = (invertColumnOrder) ? ncol - icol - 1 : icol;
        mat->setValue(irow, jcol, X[lec++]);
      }
  }
  else
  {
    for (Id icol = 0; icol < ncol; icol++)
      for (Id irow = 0; irow < nrow; irow++)
      {
        Id jcol = (invertColumnOrder) ? ncol - icol - 1 : icol;
        mat->setValue(irow, jcol, X[lec++]);
      }
  }
  return mat;
}

/**
 * LU Decomposition of a square matrix (not necessarily symmetric)
 * @param tls  Output square matrix containing lower triangle (stored columnwise)
 * @param tus  Output square matrix containing upper triangle (stored linewise)
 * @param eps  Tolerance
 *
 * @remarks The output matrices 'tus'  and 'tls' must be dimensioned beforehand
 */
Id MatrixSquare::decomposeLU(MatrixSquare& tls,
                             MatrixSquare& tus,
                             double eps)
{
  auto neq = getNRows();
  tls.fill(0.);
  tus.fill(0.);

  for (Id i = 0; i < neq; i++)
    tls.setValue(i, i, 1.);

  for (Id i = 0; i < neq; i++)
  {
    Id ip1 = i + 1;
    Id im1 = i - 1;

    for (Id j = 0; j < neq; j++)
    {
      tus.setValue(i, j, getValue(i, j));
      if (im1 >= 0)
      {
        for (Id k = 0; k <= im1; k++)
        {
          tus.setValue(i, j, tus.getValue(i, j) - tls.getValue(i, k) * tus.getValue(k, j));
        }
      }
    }
    if (ip1 < neq)
    {
      for (Id j = ip1; j < neq; j++)
      {
        tls.setValue(j, i, getValue(j, i));
        if (im1 >= 0)
        {
          for (Id k = 0; k <= im1; k++)
          {
            tls.setValue(j, i, tls.getValue(j, i) - tls.getValue(j, k) * tus.getValue(k, i));
          }
        }
        double pivot = tus.getValue(i, i);
        if (abs(pivot) < eps) return 1;
        tls.setValue(j, i, tls.getValue(j, i) / pivot);
      }
    }
  }
  return 0;
}

Id MatrixSquare::_invertLU()
{
  auto neq = getNRows();

  // Perform the LU decomposition
  MatrixSquare tls(neq);
  MatrixSquare tus(neq);
  MatrixSquare ais(neq);
  ais.fill(0.);
  if (decomposeLU(tls, tus)) return 1;

  VectorDouble b(neq);
  VectorDouble x(neq);
  for (Id i = 0; i < neq; i++)
  {
    // Preparing the right-hand side vector (column of the identity matrix)
    VH::fill(b, 0.);
    b[i] = 1.;

    if (_solveLU(tus, tls, b.data(), x.data())) return 1;

    for (Id j = 0; j < neq; j++)
      ais.setValue(i, j, x[j]);
  }

  // Copy the inverse matrix in the input matrix
  for (Id irow = 0; irow < neq; irow++)
    for (Id icol = 0; icol < neq; icol++)
      setValue(irow, icol, ais.getValue(irow, icol));
  return 0;
}

Id MatrixSquare::_solveLU(const MatrixSquare& tus,
                          const MatrixSquare& tls,
                          const double* b,
                          double* x)
{
  auto neq = getNRows();
  VectorDouble y(neq);
  if (_forwardLU(tls, b, y.data())) return 1;
  if (_backwardLU(tus, y.data(), x)) return 1;
  return 0;
}

/*****************************************************************************/
/*!
 **  Get the solution of a linear system (after LU decomposition)
 **     L * x = b
 **
 ** \return  Error returned code
 **
 ** \param[in]  tls  Square matrix containing lower triangle (stored column-wise)
 ** \param[in]  b    matrix (dimension neq)
 ** \param[in]  eps  Tolerance
 **
 ** \param[out] x    resulting matrix (dimension neq)
 **
 *****************************************************************************/
Id MatrixSquare::_forwardLU(const MatrixSquare& tls, const double* b, double* x, double eps)
{
  auto neq = getNRows();
  for (Id i = 0; i < neq; i++) x[i] = 0.;
  for (Id i = 0; i < neq; i++)
  {
    double tmp = b[i];
    for (Id j = 0; j < i; j++)
      tmp -= tls.getValue(i, j) * x[j];

    double pivot = tls.getValue(i, i);
    if (abs(pivot) < eps) return 1;
    x[i] = tmp / pivot;
  }
  return 0;
}

/*****************************************************************************/
/*!
 **  Get the solution of a linear system (after LU decomposition)
 **     U * x = b
 **
 ** \return  Error returned code
 **
 ** \param[in]  tus  square matrix containing upper triangle (stored line-wise)
 ** \param[in]  b    matrix (dimension neq)
 ** \param[in]  eps  Tolerance
 **
 ** \param[out] x    resulting matrix (dimension neq)
 **
 ** \remark As the Upper matrix is stored line-wise, it is considered
 ** \remark as the transposed version of a lower triangle
 **
 *****************************************************************************/
Id MatrixSquare::_backwardLU(const MatrixSquare& tus, const double* b, double* x, double eps)
{
  auto neq = getNRows();
  for (Id i = neq - 1; i >= 0; i--)
  {
    double tmp = b[i];
    for (Id j = i + 1; j < neq; j++)
      tmp -= tus.getValue(i, j) * x[j];

    double pivot = tus.getValue(i, i);
    if (abs(pivot) < eps) return 1;
    x[i] = tmp / pivot;
  }
  return 0;
}

MatrixSquare* prodNormMatMat(const MatrixDense* a,
                             const MatrixDense* m,
                             bool transpose)
{
  Id nrow   = (transpose) ? a->getNCols() : a->getNRows();
  auto* mat = new MatrixSquare(nrow);
  mat->prodNormMatMatInPlace(a, m, transpose);
  return mat;
}

MatrixSquare* prodNormMat(const MatrixDense& a, bool transpose)
{
  Id nsym   = (transpose) ? a.getNCols() : a.getNRows();
  auto* mat = new MatrixSquare(nsym);
  mat->prodNormMatInPlace(&a, transpose);
  return mat;
}

MatrixSquare* prodNormMatVec(const MatrixDense& a, const VectorDouble& vec, bool transpose)
{
  Id nsym   = (transpose) ? a.getNCols() : a.getNRows();
  auto* mat = new MatrixSquare(nsym);
  mat->prodNormMatVecInPlace(&a, vec, transpose);
  return mat;
}

/**
 * @brief Create a square matrix from three diagonals
 *
 * @param vecdiag Vector of information along the diagonal
 * @param vecinf  Vector of information along the line below the diagonal
 * @param vecsup  Vector of information along the line above the diagonal
 * @return Returned matrix
 */
MatrixSquare* MatrixSquare::createFromTridiagonal(const VectorDouble& vecdiag,
                                                  const VectorDouble& vecinf,
                                                  const VectorDouble& vecsup)
{
  // Preliminary checks
  if (vecdiag.size() != vecinf.size() || vecdiag.size() != vecsup.size())
  {
    messerr("Arguments 'vecdiag'(%d), 'vecinf'(%d) and 'vecsup'(%d) should have same dimension",
            vecdiag.size(), vecinf.size(), vecsup.size());
    return nullptr;
  }

  /* Initializations */
  Id neq    = static_cast<Id>(vecdiag.size());
  auto* res = new MatrixSquare(neq);
  res->fill(0.);

  for (Id i = 0; i < neq; i++)
  {
    res->setValue(i, i, vecdiag[i]);
    if (i > 0)
      res->setValue(i, i - 1, vecinf[i]); // sous-diagonale
    if (i < neq - 1)
      res->setValue(i, i + 1, vecsup[i]); // sur-diagonale
  }
  return res;
}

Id MatrixSquare::computeEigen(bool optionPositive)
{
  return MatrixDense::_computeEigen(optionPositive);
}

} // namespace gstlrn