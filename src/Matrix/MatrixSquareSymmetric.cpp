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

#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/AMatrixSquare.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/AException.hpp"

MatrixSquareSymmetric::MatrixSquareSymmetric(int nrow, int opt_eigen)
    : AMatrixSquare(nrow, opt_eigen),
      _squareSymMatrix()
{
  _allocate();
}

MatrixSquareSymmetric::MatrixSquareSymmetric(const MatrixSquareSymmetric &r) 
  : AMatrixSquare(r),
   _squareSymMatrix()
{
  _recopyLocal(r);
}

MatrixSquareSymmetric::MatrixSquareSymmetric(const AMatrix &m)
    : AMatrixSquare(m),
      _squareSymMatrix()
{
  if (!m.isSymmetric())
  {
    messerr("The input matrix should be Symmetric");
    _clear();
    return;
  }
  const MatrixSquareSymmetric* matrixLoc = dynamic_cast<const MatrixSquareSymmetric*>(&m);
  if (matrixLoc != nullptr)
    _recopyLocal(*matrixLoc);
  else
  {
    _allocate();
    AMatrix::copyElements(m);
  }
}

MatrixSquareSymmetric& MatrixSquareSymmetric::operator= (const MatrixSquareSymmetric &r)
{
  if (this != &r)
  {
    AMatrixSquare::operator=(r);
    _recopyLocal(r);
  }
  return *this;
}

MatrixSquareSymmetric::~MatrixSquareSymmetric()
{
  _deallocate();
}

/**
 * Converts a VectorVectorDouble into a Square Symmetric Matrix
 * Note: the input argument is stored by row (if coming from [] specification)
 * @param  X Input VectorVectorDouble argument
 * @param opt_eigen Option for use of Eigen Library
 * @return The returned square symmetric matrix
 *
 * @remark: the matrix is transposed implicitly while reading
 */
MatrixSquareSymmetric* MatrixSquareSymmetric::createFromVVD(const VectorVectorDouble& X, int opt_eigen)
{
  int nrow = (int) X.size();
  int ncol = (int) X[0].size();
  if (nrow != ncol)
  {
    messerr("The matrix does not seem to be square");
    return nullptr;
  }
  MatrixSquareSymmetric* mat = new MatrixSquareSymmetric(nrow, opt_eigen);
  mat->_fillFromVVD(X);
  return mat;
}

MatrixSquareSymmetric* MatrixSquareSymmetric::createFromVD(const VectorDouble &X,
                                                           int nrow,
                                                           int opt_eigen)
{
  int ncol = nrow;
  if (nrow * ncol != (int) X.size())
  {
    messerr("Inconsistency between arguments 'nrow'(%d) and 'ncol'(%d)", nrow, ncol);
    messerr("and the dimension of the input Vector (%d)", (int) X.size());
  }
  // Check symmetry
  MatrixRectangular* mattemp = MatrixRectangular::createFromVD(X, nrow, ncol);
  if (! mattemp->isSymmetric())
  {
    messerr("The input matrix does not seem to be Square and symmetric");
    delete mattemp;
    return nullptr;
  }
  delete mattemp;

  MatrixSquareSymmetric *mat = new MatrixSquareSymmetric(nrow, opt_eigen);

  int lec = 0;
  for (int irow = 0; irow < nrow; irow++)
    for (int icol = 0; icol < ncol; icol++)
      mat->setValue(irow, icol, X[lec++]);
  return mat;
}

double MatrixSquareSymmetric::_getValue(int irow, int icol) const
{
  if (_isFlagEigen())
    return AMatrixDense::_getValue(irow, icol);
  else
    return _getValueLocal(irow, icol);
}

double MatrixSquareSymmetric::_getValueByRank(int irank) const
{
  if (_isFlagEigen())
    return AMatrixDense::_getValueByRank(irank);
  else
    return _getValueLocal(irank);
}

double& MatrixSquareSymmetric::_getValueRef(int irow, int icol)
{
  if (_isFlagEigen())
    return AMatrixDense::_getValueRef(irow, icol);
  else
    return _getValueRef(irow, icol);
}

void MatrixSquareSymmetric::_setValue(int irow, int icol, double value)
{
  if (_isFlagEigen())
  {
    // Do not forget to make a symmetrical call (when stored in an Eigen format)
    AMatrixDense::_setValue(irow, icol, value);
    AMatrixDense::_setValue(icol, irow, value);
  }
  else
    _setValueLocal(irow, icol, value);
}

void MatrixSquareSymmetric::_setValueByRank(int irank, double value)
{
  if (_isFlagEigen())
    AMatrixDense::_setValueByRank(irank, value);
  else
    _setValueLocal(irank, value);
}

void MatrixSquareSymmetric::_prodVectorInPlace(const double *inv, double *outv) const
{
  if (_isFlagEigen())
    AMatrixDense::_prodVectorInPlace(inv, outv);
  else
    _prodVectorLocal(inv, outv);
}

/**
 * \warning : values is provided as a square complete matrix
 */
void MatrixSquareSymmetric::_setValues(const double* values, bool byCol)
{
  if (_isFlagEigen())
    AMatrixDense::_setValues(values, byCol);
  else
    _setValuesLocal(values, byCol);
}

int MatrixSquareSymmetric::_invert()
{
  if (_isFlagEigen())
    return AMatrixDense::_invert();
  else
    return _invertLocal();
}

void MatrixSquareSymmetric::_deallocate()
{
  if (_isFlagEigen())
    AMatrixDense::_deallocate();
  else
  {
    delete _eigenVectors;
  }
}

void MatrixSquareSymmetric::_allocate()
{
  if (_isFlagEigen())
    AMatrixDense::_allocate();
  else
    _allocateLocal();
}

int MatrixSquareSymmetric::_getIndexToRank(int irow, int icol) const
{
  if (_isFlagEigen())
    return AMatrixDense::_getIndexToRank(irow, icol);
  else
    return _getIndexToRankLocal(irow, icol);
}

int MatrixSquareSymmetric::_getMatrixPhysicalSize() const
{
  if (_isFlagEigen())
    return AMatrixDense::_getMatrixPhysicalSize();
  else
    return _getMatrixPhysicalSizeLocal();
}

int MatrixSquareSymmetric::_solve(const VectorDouble& b, VectorDouble& x) const
{
  if (_isFlagEigen())
    return AMatrixDense::_solve(b, x);
  else
    return _solveLocal(b, x);
}

bool MatrixSquareSymmetric::_isPhysicallyPresent(int irow, int icol) const
{
  if (icol >  irow) return false;
  return true;
}

/**
 * Perform the product: this = t(X) %*% X
 * @param x: Matrix [nrow, ncol] where ncol = this->getNSize()
 */
void MatrixSquareSymmetric::normSingleMatrix(const AMatrix& x)
{
  if (getNSize() != x.getNCols())
  {
    my_throw("Incompatible matrix dimensions");
  }

  int nout = getNSize();
  int n = x.getNRows();
  for (int irow = 0; irow < nout; irow++)
  {
    for (int icol = 0; icol <= irow; icol++)
    {
      double value = 0.;
      for (int k = 0; k < n; k++)
      {
        value += x.getValue(k,irow) * x.getValue(k,icol);
      }
      setValue(irow,icol,value);
    }
  }
}

/**
 * Perform the product: this = X %*% t(X)
 * @param x: Matrix [nrow, ncol] where nrow = this->getNSize()
 */
void MatrixSquareSymmetric::normTSingleMatrix(const AMatrix& x)
{
  if (getNSize() != x.getNRows())
  {
    my_throw("Incompatible matrix dimensions");
  }

  int nout = getNSize();
  int n = x.getNCols();
  for (int irow = 0; irow < nout; irow++)
  {
    for (int icol = 0; icol <= irow; icol++)
    {
      double value = 0.;
      for (int k = 0; k < n; k++)
      {
        value += x.getValue(irow,k) * x.getValue(icol,k);
      }
      setValue(irow,icol,value);
    }
  }
}

MatrixSquareSymmetric* MatrixSquareSymmetric::reduce(const VectorInt &validRows) const
{
  // Order and shrink the input vectors
  VectorInt localValidRows = VH::filter(validRows, 0, getNRows());
  int newNRows = (int) localValidRows.size();
  if (newNRows <= 0)
  {
    messerr("The new Matrix has no Row left");
    return nullptr;
  }

  MatrixSquareSymmetric* res = new MatrixSquareSymmetric(newNRows);
  res->copyReduce(this, localValidRows, localValidRows);

  return res;
}

int MatrixSquareSymmetric::computeEigen()
{
  if (_isFlagEigen())
    return AMatrixDense::_computeEigen();
  else
    return _computeEigenLocal();
}

int MatrixSquareSymmetric::computeGeneralizedEigen(const MatrixSquareSymmetric& b)
{
  if (_isFlagEigen())
    return AMatrixDense::_computeGeneralizedEigen(b);
  else
    return _computeGeneralizedEigenLocal(b);
}

/// =============================================================================
/// The subsequent methods rely on the specific local storage ('squareSymMatrix')
/// =============================================================================

void MatrixSquareSymmetric::_recopyLocal(const MatrixSquareSymmetric& r)
{
  _squareSymMatrix = r._squareSymMatrix;
  _flagEigenDecompose = r._flagEigenDecompose;
  _eigenValues = r._eigenValues;
  if (_eigenVectors != nullptr) _eigenVectors = r._eigenVectors->clone();
}

double MatrixSquareSymmetric::_getValueLocal(int irow, int icol) const
{
  if (! _isIndexValid(irow,icol)) return TEST;
  int rank = _getIndexToRank(irow,icol);
  return _squareSymMatrix[rank];
}

double MatrixSquareSymmetric::_getValueLocal(int irank) const
{
  return _squareSymMatrix[irank];
}

double& MatrixSquareSymmetric::_getValueRefLocal(int irow, int icol)
{
  int rank = _getIndexToRank(irow, icol);
  return _squareSymMatrix[rank];
}

void MatrixSquareSymmetric::_setValueLocal(int irow, int icol, double value)
{
  if (! _isIndexValid(irow, icol)) return;
  int irank = _getIndexToRank(irow, icol);
  _squareSymMatrix[irank] = value;
}

void MatrixSquareSymmetric::_setValueLocal(int irank, double value)
{
  if (! _isRankValid(irank)) return;
  _squareSymMatrix[irank] = value;
}

void MatrixSquareSymmetric::_prodVectorLocal(const double *inv, double *outv) const
{
  matrix_triangular_product(getNRows(),2,_squareSymMatrix.data(),inv,outv);
}

/**
 * \warning : values is provided as a square complete matrix
 */
void MatrixSquareSymmetric::_setValuesLocal(const double *values, bool byCol)
{
  // Check that the input argument corresponds to a square symmetric matrix
  for (int icol = 0; icol < getNCols(); icol++)
    for (int irow = 0; irow < getNRows(); irow++)
    {
      double val1 = values[icol * getNRows() + irow];
      double val2 = values[irow * getNCols() + icol];
      if (ABS(val1 - val2) > EPSILON10)
      {
        messerr(
            "Argument 'values' must correspond to a Square Symmetric Matrix");
        messerr("- Element[%d,%d] = %lf", icol, irow, val1);
        messerr("- Element(%d,%d) = %lf", irow, icol, val2);
        messerr("Operation is aborted");
        return;
      }
    }

  int ecr = 0;
  for (int icol = 0; icol < getNCols(); icol++)
    for (int irow = 0; irow < getNRows(); irow++, ecr++)
    {
      setValue(irow, icol, values[ecr]);
    }
}

int MatrixSquareSymmetric::_invertLocal()
{
  return matrix_invert_triangle(getNRows(),_squareSymMatrix.data(), -1);
}

void MatrixSquareSymmetric::_allocateLocal()
{
  _squareSymMatrix.resize(_getMatrixPhysicalSize());
}

int MatrixSquareSymmetric::_getIndexToRankLocal(int irow, int icol) const
{
  int rank;
  int n = getNRows();
  if (irow >= icol)
    rank = icol * n + irow - icol * (icol + 1) / 2;
  else
    rank = irow * n + icol - irow * (irow + 1) / 2;
  return rank;
}

int MatrixSquareSymmetric::_getMatrixPhysicalSizeLocal() const
{
  int n = getNRows();
  int size = n * (n + 1) / 2;
  return size;
}

int MatrixSquareSymmetric::_solveLocal(const VectorDouble& b, VectorDouble& x) const
{
  int pivot;
  return matrix_solve(1,_squareSymMatrix.data(),b.data(),x.data(),
                      static_cast<int> (b.size()),1,&pivot);
}

int MatrixSquareSymmetric::_computeEigenLocal()
{
  int nrows = getNRows();
  _flagEigenDecompose = true;
  _eigenValues = VectorDouble(nrows, 0.);
  VectorDouble eigenVectors = VectorDouble(nrows * nrows, 0);

  int err = matrix_eigen(this->getValues().data(), nrows, _eigenValues.data(), eigenVectors.data());
  if (err == 0)
  {
    _eigenVectors = MatrixSquareGeneral::createFromVD(eigenVectors, nrows, false, 0, false);
  }

  return err;
}

int MatrixSquareSymmetric::_computeGeneralizedEigenLocal(const MatrixSquareSymmetric& b)
{
  int nrows = getNRows();
  _flagEigenDecompose = true;
  _eigenValues = VectorDouble(nrows, 0.);
  VectorDouble eigenVectors = VectorDouble(nrows * nrows, 0);

  int err = matrix_geigen(this->getValues().data(),b.getValues().data(), nrows, _eigenValues.data(), eigenVectors.data());
  if (err == 0)
  {
    std::reverse(_eigenValues.begin(), _eigenValues.end());
    _eigenVectors = MatrixSquareGeneral::createFromVD(eigenVectors, nrows, false, 0, true);
  }
  return err;
}
