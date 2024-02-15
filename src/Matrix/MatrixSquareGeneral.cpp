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
#include "geoslib_old_f.h"

#include "Matrix/MatrixSquareGeneral.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/AException.hpp"

MatrixSquareGeneral::MatrixSquareGeneral(int nrow, int opt_eigen)
  : AMatrixSquare(nrow, opt_eigen)
  , _squareMatrix()
{
  _allocate();
}

MatrixSquareGeneral::MatrixSquareGeneral(const MatrixSquareGeneral &r) 
  : AMatrixSquare(r),
    _squareMatrix()
{
  _recopyLocal(r);
}

MatrixSquareGeneral::MatrixSquareGeneral(const AMatrix &m)
  : AMatrixSquare(m),
    _squareMatrix()
{
  const MatrixSquareGeneral* matrixLoc = dynamic_cast<const MatrixSquareGeneral*>(&m);
  if (matrixLoc != nullptr)
    _recopyLocal(*matrixLoc);
  else
  {
    _allocate();
    AMatrix::copyElements(m);
  }
}

MatrixSquareGeneral& MatrixSquareGeneral::operator= (const MatrixSquareGeneral &r)
{
  if (this != &r)
  {
    AMatrixSquare::operator=(r);
    _recopyLocal(r);
  }
  return *this;
}

MatrixSquareGeneral::~MatrixSquareGeneral()
{
  _deallocate();
}

/**
 * Converts a VectorVectorDouble into a Matrix
 * Note: the input argument is stored by row (if coming from [] specification)
 * @param X Input VectorVectorDouble argument
 * @param opt_eigen Option for use of Eigen Library
 * @return The returned square matrix
 *
 * @remark: the matrix is transposed implicitly while reading
 */
MatrixSquareGeneral* MatrixSquareGeneral::createFromVVD(const VectorVectorDouble& X, int opt_eigen)
{
  int nrow = (int) X.size();
  int ncol = (int) X[0].size();
  if (nrow != ncol)
  {
    messerr("The matrix does not seem to be square");
    return nullptr;
  }

  MatrixSquareGeneral* mat = new MatrixSquareGeneral(nrow, opt_eigen);
  mat->_fillFromVVD(X);
  return mat;
}

MatrixSquareGeneral* MatrixSquareGeneral::createFromVD(const VectorDouble &X,
                                                       int nrow,
                                                       bool byCol,
                                                       int opt_eigen,
                                                       bool invertColumnOrder)
{
  int ncol = nrow;
  if (nrow * ncol != (int) X.size())
  {
    messerr("Inconsistency between arguments 'nrow'(%d) and 'ncol'(%d)", nrow, ncol);
    messerr("and the dimension of the input Vector (%d)", (int) X.size());
  }
  MatrixSquareGeneral* mat = new MatrixSquareGeneral(nrow, opt_eigen);

  int lec = 0;
  if (byCol)
  {
    for (int irow = 0; irow < nrow; irow++)
      for (int icol = 0; icol < ncol; icol++)
      {
        int jcol = (invertColumnOrder) ? ncol - icol - 1 : icol;
        mat->setValue(irow, jcol, X[lec++]);
      }
  }
  else
  {
    for (int icol = 0; icol < ncol; icol++)
      for (int irow = 0; irow < nrow; irow++)
      {
        int jcol = (invertColumnOrder) ? ncol - icol - 1 : icol;
        mat->setValue(irow, jcol, X[lec++]);
      }
  }
  return mat;
}

double MatrixSquareGeneral::_getValue(int irow, int icol) const
{
  if (isFlagEigen())
    return AMatrixDense::_getValue(irow, icol);
  else
    return _getValueLocal(irow, icol);
}

double MatrixSquareGeneral::_getValueByRank(int irank) const
{
  if (isFlagEigen())
    return AMatrixDense::_getValueByRank(irank);
  else
    return _getValueLocal(irank);
}

double& MatrixSquareGeneral::_getValueRef(int irow, int icol)
{
  if (isFlagEigen())
    return AMatrixDense::_getValueRef(irow, icol);
  else
    return _getValueRefLocal(irow, icol);
}

void MatrixSquareGeneral::_setValue(int irow, int icol, double value)
{
  if (isFlagEigen())
    AMatrixDense::_setValue(irow, icol, value);
  else
    _setValueLocal(irow, icol, value);
}

void MatrixSquareGeneral::_setValueByRank(int irank, double value)
{
  if (isFlagEigen())
    AMatrixDense::_setValueByRank(irank, value);
  else
    _setValueLocal(irank, value);
}

/**
 * Returns 'y' = 'this' %*% 'x'
 * @param x  Input Vector
 * @param y Output Vector
 * @param transpose True if the matrix 'this' must be transposed
 */
void MatrixSquareGeneral::_prodMatVecInPlacePtr(const double *x, double *y, bool transpose) const
{
  if (isFlagEigen())
    AMatrixDense::_prodMatVecInPlacePtr(x, y, transpose);
  else
    _prodMatVecInPlacePtrLocal(x, y, transpose);
}

/**
 * Returns 'y' = 'this' %*% 'x'
 * @param x  Input Vector
 * @param y Output Vector
 * @param transpose True if the matrix 'this' must be transposed
 */
void MatrixSquareGeneral::_prodVecMatInPlacePtr(const double *x, double *y, bool transpose) const
{
  if (isFlagEigen())
    AMatrixDense::_prodVecMatInPlacePtr(x, y, transpose);
  else
    _prodVecMatInPlacePtrLocal(x, y, transpose);
}

void MatrixSquareGeneral::_transposeInPlace()
{
  if (isFlagEigen())
    AMatrixDense::_transposeInPlace();
  else
    _transposeInPlaceLocal();
}

int MatrixSquareGeneral::_invert()
{
  if (isFlagEigen())
    return AMatrixDense::_invert();
  else
    return _invertLocal();
}

void MatrixSquareGeneral::_deallocate()
{
  if (isFlagEigen())
    AMatrixDense::_deallocate();
  else
  {
    // Specific code for this class should be placed here
  }
}

void MatrixSquareGeneral::_allocate()
{
  if (isFlagEigen())
    AMatrixDense::_allocate();
  else
    _allocateLocal();
}

int MatrixSquareGeneral::_getMatrixPhysicalSize() const
{
  return(getNRows() * getNCols());
}

int MatrixSquareGeneral::_solve(const VectorDouble& /*b*/, VectorDouble& /*x*/) const
{
  my_throw("Invert method is limited to Square Symmetrical Matrices");
  return 0;
}

MatrixSquareGeneral* MatrixSquareGeneral::createReduce(const VectorInt &validRows) const
{
  // Order and shrink the input vectors
  VectorInt localValidRows = VH::filter(validRows, 0, getNRows());
  int newNRows = (int) localValidRows.size();
  if (newNRows <= 0)
  {
    messerr("The new Matrix has no Row left");
    return nullptr;
  }

  MatrixSquareGeneral* res = new MatrixSquareGeneral(newNRows);
  res->copyReduce(this, localValidRows, localValidRows);

  return res;
}

/// ==========================================================================
/// The subsequent methods rely on the specific local storage ('squareMatrix')
/// ==========================================================================

void MatrixSquareGeneral::_allocateLocal()
{
  _squareMatrix.resize(_getMatrixPhysicalSize());
}

void MatrixSquareGeneral::_recopyLocal(const MatrixSquareGeneral &r)
{
  _squareMatrix = r._squareMatrix;
}

double MatrixSquareGeneral::_getValueLocal(int irow, int icol) const
{
  if (!_isIndexValid(irow, icol)) return TEST;
  int rank = _getIndexToRank(irow, icol);
  return _squareMatrix[rank];
}

double MatrixSquareGeneral::_getValueLocal(int irank) const
{
  return _squareMatrix[irank];
}

double& MatrixSquareGeneral::_getValueRefLocal(int irow, int icol)
{
  int rank = _getIndexToRank(irow,icol);
  return _squareMatrix[rank];
}

void MatrixSquareGeneral::_setValueLocal(int irow, int icol, double value)
{
  if (!_isIndexValid(irow, icol)) return;
  int rank = _getIndexToRank(irow, icol);
  _squareMatrix[rank] = value;
}

void MatrixSquareGeneral::_setValueLocal(int irank, double value)
{
  if (!_isRankValid(irank)) return;
  _squareMatrix[irank] = value;
}

void MatrixSquareGeneral::_prodMatVecInPlacePtrLocal(const double *x, double *y, bool transpose) const
{
  if (transpose)
    matrix_product_safe(1, getNRows(), getNCols(), x, _squareMatrix.data(), y);
  else
    matrix_product_safe(getNRows(), getNCols(), 1, _squareMatrix.data(), x, y);
}

void MatrixSquareGeneral::_prodVecMatInPlacePtrLocal(const double *x, double *y, bool transpose) const
{
  if (transpose)
    matrix_product_safe(getNRows(), getNCols(), 1, _squareMatrix.data(), x, y);
  else
    matrix_product_safe(1, getNRows(), getNCols(), x, _squareMatrix.data(), y);
}

void MatrixSquareGeneral::_transposeInPlaceLocal()
{
  int nrow = getNRows();
  int ncol = getNCols();
  VectorDouble old = _squareMatrix;
  matrix_transpose(nrow, ncol, _squareMatrix.data(), old.data());
  _squareMatrix = old;
  _setNCols(nrow);
  _setNRows(ncol);
}

int MatrixSquareGeneral::_invertLocal()
{
  if (getNRows() <= 3)
    return matrix_invreal(_squareMatrix.data(), getNRows());
  else
  {
    int error = matrix_LU_invert(getNRows(), _squareMatrix.data());
    return error;
  }
}

MatrixSquareGeneral* prodNormMatMat(const AMatrixDense &a,
                                    const AMatrixDense &m,
                                    bool transpose)
{
  int nrow = (transpose) ? a.getNCols() : a.getNRows();
  MatrixSquareGeneral *mat = new MatrixSquareGeneral(nrow, a.isFlagEigen());
  mat->prodNormMatMatInPlace(a, m, transpose);
  return mat;
}

MatrixSquareGeneral* prodNormMat(const AMatrixDense &a, const VectorDouble& vec, bool transpose)
{
  int nsym = (transpose) ? a.getNCols() : a.getNRows();
  MatrixSquareGeneral *mat = new MatrixSquareGeneral(nsym, a.isFlagEigen());
  mat->prodNormMatInPlace(a, vec, transpose);
  return mat;
}

