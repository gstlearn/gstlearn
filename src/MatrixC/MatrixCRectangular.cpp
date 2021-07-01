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
#include "MatrixC/MatrixCRectangular.hpp"
#include "MatrixC/AMatrixC.hpp"
#include "Basic/AException.hpp"
#include "geoslib_f.h"

MatrixCRectangular::MatrixCRectangular(int nrows, int ncols, bool sparse)
    : AMatrixC(nrows, ncols, sparse),
      _rectMatrix()
{
  _allocate();
}

MatrixCRectangular::MatrixCRectangular(const MatrixCRectangular &r)
  : AMatrixC(r)
{
  _recopy(r);
}

MatrixCRectangular& MatrixCRectangular::operator= (const MatrixCRectangular &r)
{
  if (this != &r)
  {
    _deallocate();
    AMatrixC::operator=(r);
    _recopy(r);
  }
  return *this;
}

MatrixCRectangular::~MatrixCRectangular()
{
  _deallocate();
}

IClonable* MatrixCRectangular::clone() const
{
  return new MatrixCRectangular(*this);
}

double MatrixCRectangular::_getValue(int irow, int icol) const
{
  _isIndexValid(irow,icol);
  int rank = _getIndexToRank(irow,icol);
  return _rectMatrix[rank];
}

double MatrixCRectangular::_getValue(int irank) const
{
  _isRankValid(irank);
  return _rectMatrix[irank];
}

void MatrixCRectangular::_setValue(int irank, double value)
{
  _isRankValid(irank);
  _rectMatrix[irank] = value;
}

void MatrixCRectangular::_setValue(int irow, int icol, double value)
{
  _isIndexValid(irow, icol);
  int rank = _getIndexToRank(irow, icol);
  _rectMatrix[rank] = value;
}

void MatrixCRectangular::_prodVector(const double *in, double *out) const
{
  matrix_product_safe(getNRows(), getNCols(), 1, _rectMatrix.data(), in, out);
}

void MatrixCRectangular::_transposeInPlace()
{
  if (isSparse())
  {
    AMatrixC::transposeInPlace();
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

void MatrixCRectangular::_setValues(const double* values, bool byCol)
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
double& MatrixCRectangular::_getValueRef(int irow, int icol)
{
  _isIndexValid(irow,icol);
  int rank = _getIndexToRank(irow,icol);
  return _rectMatrix[rank];
}

void MatrixCRectangular::_deallocate()
{

}

void MatrixCRectangular::_recopy(const MatrixCRectangular &r)
{
  _rectMatrix = r._rectMatrix;
}

void MatrixCRectangular::_allocate()
{
  _rectMatrix.resize(_getMatrixSize(),0.);
}

int MatrixCRectangular::_getIndexToRank(int irow,
                                        int icol) const
{
  int rank = icol * getNRows() + irow;
  return rank;
}

int MatrixCRectangular::_getMatrixSize() const
{
  return (getNRows() * getNCols());
}

int MatrixCRectangular::_invert()
{
  my_throw("Invert method is limited to Square matrices");
  return 0;
}

int MatrixCRectangular::_solve(const VectorDouble& b, VectorDouble& x) const
{
  my_throw("Invert method is limited to Square Symmetrical Matrices");
  return 0;
}

double MatrixCRectangular::_determinant() const
{
  my_throw("Determinant method is limited to Square Matrices");
  return TEST;
}

