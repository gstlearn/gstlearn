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

double MatrixSquareGeneral::_getValue(int irow, int icol) const
{
  if (_isFlagEigen())
    return AMatrixDense::_getValue(irow, icol);
  else
    return _getValueLocal(irow, icol);
}

double MatrixSquareGeneral::_getValueByRank(int irank) const
{
  if (_isFlagEigen())
    return AMatrixDense::_getValueByRank(irank);
  else
    return _getValueLocal(irank);
}

double& MatrixSquareGeneral::_getValueRef(int irow, int icol)
{
  if (_isFlagEigen())
    return AMatrixDense::_getValueRef(irow, icol);
  else
    return _getValueRefLocal(irow, icol);
}

void MatrixSquareGeneral::_setValue(int irow, int icol, double value)
{
  if (_isFlagEigen())
    AMatrixDense::_setValue(irow, icol, value);
  else
    _setValueLocal(irow, icol, value);
}

void MatrixSquareGeneral::_setValueByRank(int irank, double value)
{
  if (_isFlagEigen())
    AMatrixDense::_setValueByRank(irank, value);
  else
    _setValueLocal(irank, value);
}

/**
 * Right product of this by in gives out
 * @param inv  Input Vector
 * @param outv Output Vector
 */
void MatrixSquareGeneral::_prodVectorInPlace(const double *inv, double *outv) const
{
  if (_isFlagEigen())
    AMatrixDense::_prodVectorInPlace(inv, outv);
  else
    _prodVectorLocal(inv, outv);
}

void MatrixSquareGeneral::_transposeInPlace()
{
  if (_isFlagEigen())
    AMatrixDense::_transposeInPlace();
  else
    _transposeInPlaceLocal();
}

int MatrixSquareGeneral::_invert()
{
  if (_isFlagEigen())
    return AMatrixDense::_invert();
  else
    return _invertLocal();
}

void MatrixSquareGeneral::_deallocate()
{
  if (_isFlagEigen())
    AMatrixDense::_deallocate();
  else
  {
    // Specific code for this class should be placed here
  }
}

void MatrixSquareGeneral::_allocate()
{
  if (_isFlagEigen())
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

/**
 * Converts a VectorVectorDouble into a Matrix
 * Note: the input argument is stored by row (if coming from [] specification)
 * @param X Input VectorVectorDouble argument
 * @return The returned matrix
 *
 * @remark: the matrix is transposed implicitly while reading
 */
MatrixSquareGeneral* MatrixSquareGeneral::createFromVVD(const VectorVectorDouble& X)
{
  int nrow = (int) X.size();
  int ncol = (int) X[0].size();
  if (nrow != ncol)
  {
    messerr("The matrix does not seem to be square");
    return nullptr;
  }

  MatrixSquareGeneral* mat = new MatrixSquareGeneral(nrow);
  mat->_fillFromVVD(X);
  return mat;
}

MatrixSquareGeneral* MatrixSquareGeneral::reduce(const VectorInt &validRows) const
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

/**
 * Right product of this by in gives out
 * @param inv  Input Vector
 * @param outv Output Vector
 */
void MatrixSquareGeneral::_prodVectorLocal(const double *inv, double *outv) const
{
  matrix_product_safe(getNRows(),getNCols(),1,_squareMatrix.data(),inv,outv);
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
