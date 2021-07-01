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
#include <MatrixC/MatrixCSGeneral.hpp>
#include "geoslib_f.h"

#include "Basic/AException.hpp"

MatrixCSGeneral::MatrixCSGeneral(int nrow, bool sparse)
  : AMatrixCSquare(nrow, sparse)
  , _squareMatrix()
{
  _allocate();
}

MatrixCSGeneral::MatrixCSGeneral(const MatrixCSGeneral &r) 
  : AMatrixCSquare(r)
{
  _recopy(r);
}

MatrixCSGeneral::MatrixCSGeneral(const AMatrixC &m)
{
  if (m.isEmpty())
    my_throw("The input matrix should be non-empty");
  if (!m.isSquare())
    my_throw("The input matrix should be Square");

  _setNRows(m.getNRows());
  _setNCols(m.getNCols());
  _allocate();
  for (int icol = 0; icol < m.getNCols(); icol++)
    for (int irow = 0; irow < m.getNRows(); irow++)
    {
      setValue(irow, icol, m.getValue(irow, icol));
    }
}

MatrixCSGeneral& MatrixCSGeneral::operator= (const MatrixCSGeneral &r)
{
  if (this != &r)
  {
    _deallocate();
    AMatrixCSquare::operator=(r);
    _recopy(r);
  }
  return *this;
}

MatrixCSGeneral::~MatrixCSGeneral()
{
  _deallocate();
}

IClonable* MatrixCSGeneral::clone() const
{
  return new MatrixCSGeneral(*this);
}

double MatrixCSGeneral::_getValue(int irow, int icol) const
{
  _isIndexValid(irow,icol);
  int rank = _getIndexToRank(irow,icol);
  return _squareMatrix[rank];
}

double MatrixCSGeneral::_getValue(int irank) const
{
  return _squareMatrix[irank];
}

void MatrixCSGeneral::_setValue(int irow, int icol, double value)
{
  _isIndexValid(irow,icol);
  int rank = _getIndexToRank(irow,icol);
  _squareMatrix[rank] = value;
}

void MatrixCSGeneral::_setValue(int irank, double value)
{
  _isRankValid(irank);
  _squareMatrix[irank] = value;
}
/**
 * Right product of this by in gives out
 * @param in  Input Vector
 * @param out Output Vector
 */
void MatrixCSGeneral::_prodVector(const double *in, double *out) const
{
  int nrow = getNRows();
  int ncol = getNCols();
  matrix_product(nrow,ncol,1,_squareMatrix.data(),in,out);
}

void MatrixCSGeneral::_transposeInPlace()
{
  if (isSparse())
  {
    AMatrixC::transposeInPlace();
  }
  else
  {
    int nrow = getNRows();
    int ncol = getNCols();
    VectorDouble old = _squareMatrix;
    matrix_transpose(nrow, ncol, _squareMatrix.data(), old.data());
    _squareMatrix = old;
    _setNCols(nrow);
    _setNRows(ncol);
  }
}

void MatrixCSGeneral::_setValues(const double* values, bool byCol)
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

double MatrixCSGeneral::_determinant() const
{
  return matrix_determinant(getNRows(),_squareMatrix.data());
}

int MatrixCSGeneral::_invert()
{
  return matrix_invreal(_squareMatrix.data(),getNRows());
}

double& MatrixCSGeneral::_getValueRef(int irow, int icol)
{
  _isIndexValid(irow,icol);
  int rank = _getIndexToRank(irow,icol);
  return _squareMatrix[rank];
}

void MatrixCSGeneral::_deallocate()
{
}

void MatrixCSGeneral::_recopy(const MatrixCSGeneral &r)
{
  _squareMatrix = r._squareMatrix;
}

void MatrixCSGeneral::_allocate()
{
  _squareMatrix.resize(_getMatrixSize());
}

int MatrixCSGeneral::_getIndexToRank(int irow, int icol) const
{
  int rank = irow * getNCols() + icol;
  return rank;
}

int MatrixCSGeneral::_getMatrixSize() const
{
  return(getNRows() * getNCols());
}

int MatrixCSGeneral::_solve(const VectorDouble& b, VectorDouble& x) const
{
  my_throw("Invert method is limited to Square Symmetrical Matrices");
  return 0;
}
