/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_old_f.h"

#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/AMatrix.hpp"
#include "Basic/AException.hpp"
#include "csparse_d.h"

MatrixRectangular::MatrixRectangular(int nrows, int ncols, bool sparse)
    : AMatrix(nrows, ncols, sparse),
      _rectMatrix()
{
  _allocate();
}

MatrixRectangular::MatrixRectangular(const MatrixRectangular &r)
  : AMatrix(r),
    _rectMatrix(r._rectMatrix)
{
}

MatrixRectangular& MatrixRectangular::operator= (const MatrixRectangular &r)
{
  if (this != &r)
  {
    AMatrix::operator=(r);
    _rectMatrix = r._rectMatrix;
  }
  return *this;
}

MatrixRectangular::~MatrixRectangular()
{
  _deallocate();
}

/**
 * Converts a VectorVectorDouble into a Matrix
 * Note: the input argument is stored by row (if coming from [] specification)
 * @param X Input VectorVectorDouble argument
 * @param sparse True for a Sparse matrix
 * @return The returned rectangular matrix
 *
 * @remark: the matrix is transposed implicitly while reading
 */
MatrixRectangular* MatrixRectangular::createFromVVD(const VectorVectorDouble& X, bool sparse)
{
  int nrow = (int) X.size();
  int ncol = (int) X[0].size();

  MatrixRectangular* mat = new MatrixRectangular(nrow, ncol, sparse);
  for (int irow = 0; irow < nrow; irow++)
    for (int icol = 0; icol < ncol; icol++)
      mat->setValue(irow, icol, X[irow][icol]);

  return mat;
}

MatrixRectangular* MatrixRectangular::createFromVD(const VectorDouble &X,
                                                   int nrow,
                                                   int ncol,
                                                   bool byCol,
                                                   bool sparse)
{
  if (nrow * ncol != (int) X.size())
  {
    messerr("Inconsistency between arguments 'nrow'(%d) and 'ncol'(%d)", nrow, ncol);
    messerr("and the dimension of the input Vector (%d)", (int) X.size());
  }
  MatrixRectangular* mat = new MatrixRectangular(nrow, ncol, sparse);

  int lec = 0;
  if (byCol)
  {
    for (int irow = 0; irow < nrow; irow++)
      for (int icol = 0; icol < ncol; icol++)
        mat->setValue(irow, icol, X[lec++]);
  }
  else
  {
    for (int icol = 0; icol < ncol; icol++)
      for (int irow = 0; irow < nrow; irow++)
        mat->setValue(irow, icol, X[lec++]);
  }
  return mat;
}

double MatrixRectangular::_getValue(int irow, int icol) const
{
  _isIndexValid(irow,icol);
  int rank = _getIndexToRank(irow,icol);
  return _rectMatrix[rank];
}

double MatrixRectangular::_getValue(int irank) const
{
  _isRankValid(irank);
  return _rectMatrix[irank];
}

void MatrixRectangular::_setValue(int irank, double value)
{
  _isRankValid(irank);
  _rectMatrix[irank] = value;
}

void MatrixRectangular::_setValue(int irow, int icol, double value)
{
  _isIndexValid(irow, icol);
  int rank = _getIndexToRank(irow, icol);
  _rectMatrix[rank] = value;
}

void MatrixRectangular::_prodVector(const double *inv, double *outv) const
{
  matrix_product_safe(getNRows(), getNCols(), 1, _rectMatrix.data(), inv, outv);
}

void MatrixRectangular::_transposeInPlace()
{
  if (isSparse())
  {
    AMatrix::transposeInPlace();
  }
  else
  {
    VectorDouble old;
    old.resize(getNRows() * getNCols());
    matrix_transpose(getNRows(), getNCols(), _rectMatrix.data(), old.data());
    _rectMatrix = old;
    int temp = getNCols();
    _setNCols(getNRows());
    _setNRows(temp);
  }
}

void MatrixRectangular::_setValues(const double* values, bool byCol)
{
  if (byCol)
  {
    int ecr = 0;
    for (int icol = 0; icol < getNCols(); icol++)
      for (int irow = 0; irow < getNRows(); irow++, ecr++)
      {
        setValue(irow, icol, values[ecr]);
      }
  }
  else
  {
    int ecr = 0;
    for (int irow = 0; irow < getNRows(); irow++)
      for (int icol = 0; icol < getNCols(); icol++, ecr++)
      {
        setValue(irow, icol, values[ecr]);
      }
  }
}

/*! Gets a reference to the value at row 'irow' and column 'icol' */
double& MatrixRectangular::_getValueRef(int irow, int icol)
{
  _isIndexValid(irow,icol);
  int rank = _getIndexToRank(irow,icol);
  return _rectMatrix[rank];
}

void MatrixRectangular::_deallocate()
{

}

void MatrixRectangular::_allocate()
{
  _rectMatrix.resize(_getMatrixSize(),0.);
}

int MatrixRectangular::_getIndexToRank(int irow, int icol) const
{
  int rank = icol * getNRows() + irow;
  return rank;
}

int MatrixRectangular::_getMatrixSize() const
{
  return (getNRows() * getNCols());
}

int MatrixRectangular::_invert()
{
  my_throw("Invert method is limited to Square matrices");
  return 0;
}

int MatrixRectangular::_solve(const VectorDouble& /*b*/, VectorDouble& /*x*/) const
{
  my_throw("Invert method is limited to Square Symmetrical Matrices");
  return 0;
}

double MatrixRectangular::_determinant() const
{
  my_throw("Determinant method is limited to Square Matrices");
  return TEST;
}

