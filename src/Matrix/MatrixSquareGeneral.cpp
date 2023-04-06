/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_old_f.h"

#include "Matrix/MatrixSquareGeneral.hpp"
#include "Basic/AException.hpp"

MatrixSquareGeneral::MatrixSquareGeneral(int nrow, bool sparse)
  : AMatrixSquare(nrow, sparse)
  , _squareMatrix()
{
  _allocate();
}

MatrixSquareGeneral::MatrixSquareGeneral(const MatrixSquareGeneral &r) 
  : AMatrixSquare(r)
{
  _recopy(r);
}

MatrixSquareGeneral::MatrixSquareGeneral(const AMatrix &m)
{
  if (m.isEmpty())
  {
    messerr("The input matrix should be non-empty");
    _clear();
    return;
  }
  if (!m.isSquare())
  {
    messerr("The input matrix should be Square");
    _clear();
    return;
  }
  _setNRows(m.getNRows());
  _setNCols(m.getNCols());
  _allocate();
  for (int icol = 0; icol < m.getNCols(); icol++)
    for (int irow = 0; irow < m.getNRows(); irow++)
    {
      setValue(irow, icol, m.getValue(irow, icol));
    }
}

MatrixSquareGeneral& MatrixSquareGeneral::operator= (const MatrixSquareGeneral &r)
{
  if (this != &r)
  {
    _deallocate();
    AMatrixSquare::operator=(r);
    _recopy(r);
  }
  return *this;
}

MatrixSquareGeneral::~MatrixSquareGeneral()
{
  _deallocate();
}

double MatrixSquareGeneral::_getValue(int irow, int icol) const
{
  _isIndexValid(irow,icol);
  int rank = _getIndexToRank(irow,icol);
  return _squareMatrix[rank];
}

double MatrixSquareGeneral::_getValue(int irank) const
{
  return _squareMatrix[irank];
}

void MatrixSquareGeneral::_setValue(int irow, int icol, double value)
{
  _isIndexValid(irow,icol);
  int rank = _getIndexToRank(irow,icol);
  _squareMatrix[rank] = value;
}

void MatrixSquareGeneral::_setValue(int irank, double value)
{
  _isRankValid(irank);
  _squareMatrix[irank] = value;
}
/**
 * Right product of this by in gives out
 * @param inv  Input Vector
 * @param outv Output Vector
 */
void MatrixSquareGeneral::_prodVector(const double *inv, double *outv) const
{
  int nrow = getNRows();
  int ncol = getNCols();
  matrix_product(nrow,ncol,1,_squareMatrix.data(),inv,outv);
}

void MatrixSquareGeneral::_transposeInPlace()
{
  if (isSparse())
  {
    AMatrix::transposeInPlace();
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

void MatrixSquareGeneral::_setValues(const double* values, bool byCol)
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

int MatrixSquareGeneral::_invert()
{
  return matrix_invreal(_squareMatrix.data(),getNRows());
}

double& MatrixSquareGeneral::_getValueRef(int irow, int icol)
{
  _isIndexValid(irow,icol);
  int rank = _getIndexToRank(irow,icol);
  return _squareMatrix[rank];
}

void MatrixSquareGeneral::_deallocate()
{
}

void MatrixSquareGeneral::_recopy(const MatrixSquareGeneral &r)
{
  _squareMatrix = r._squareMatrix;
}

void MatrixSquareGeneral::_allocate()
{
  _squareMatrix.resize(_getMatrixSize());
}

int MatrixSquareGeneral::_getIndexToRank(int irow, int icol) const
{
  // TODO We must check the impact of this modification
  int rank = irow * getNCols() + icol;
//  int rank = icol * getNRows() + irow;
  return rank;
}

int MatrixSquareGeneral::_getMatrixSize() const
{
  return(getNRows() * getNCols());
}

int MatrixSquareGeneral::_solve(const VectorDouble& /*b*/, VectorDouble& /*x*/) const
{
  my_throw("Invert method is limited to Square Symmetrical Matrices");
  return 0;
}

