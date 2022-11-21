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

#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/AMatrixSquare.hpp"
#include "Basic/AException.hpp"

MatrixSquareSymmetric::MatrixSquareSymmetric(int nrow, bool sparse)
: AMatrixSquare(nrow, sparse)
, _squareSymMatrix()
{
  _allocate();
}

MatrixSquareSymmetric::MatrixSquareSymmetric(const MatrixSquareSymmetric &r) 
  : AMatrixSquare(r)
{
  _recopy(r);
}

MatrixSquareSymmetric::MatrixSquareSymmetric(const AMatrix &m)
    : AMatrixSquare(m.getNRows(), m.isSparse()),
      _squareSymMatrix()
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
  if (!m.isSymmetric())
  {
    messerr("The input matrix should be Symmetric");
    _clear();
    return;
  }

  _setNRows(m.getNRows());
  _setNCols(m.getNCols());
  _allocate();
  for (int icol = 0; icol < m.getNCols(); icol++)
    for (int irow = 0; irow < m.getNRows(); irow++)
    {
      setValue(irow,icol,m.getValue(irow,icol));
    }
}

MatrixSquareSymmetric& MatrixSquareSymmetric::operator= (const MatrixSquareSymmetric &r)
{
  if (this != &r)
  {
    _deallocate();
    AMatrixSquare::operator=(r);
    _recopy(r);
  }
  return *this;
}

MatrixSquareSymmetric::~MatrixSquareSymmetric()
{
  _deallocate();
}

double MatrixSquareSymmetric::_getValue(int irow, int icol) const
{
  _isIndexValid(irow,icol);
  int rank = _getIndexToRank(irow,icol);
  return _squareSymMatrix[rank];
}

double MatrixSquareSymmetric::_getValue(int irank) const
{
  return _squareSymMatrix[irank];
}

void MatrixSquareSymmetric::_setValue(int irow, int icol, double value)
{
  _isIndexValid(irow, icol);
  int irank = _getIndexToRank(irow, icol);
  _squareSymMatrix[irank] = value;
}

void MatrixSquareSymmetric::_setValue(int irank, double value)
{
  _isRankValid(irank);
  _squareSymMatrix[irank] = value;
}

/**
 * Loading values from an Upper Triangular matrix
 * @param nsize
 * @param tab
 */
void MatrixSquareSymmetric::initMatTri(int     nsize,
                             double* tab)
{
  _isNumberValid(nsize,nsize);
  _setNSize(nsize);
  _allocate();
  for (int i=0; i<_getMatrixSize(); i++) _squareSymMatrix[i] = tab[i];
}

void MatrixSquareSymmetric::_prodVector(const double *inv, double *outv) const
{
  matrix_triangular_product(getNRows(),2,_squareSymMatrix.data(),inv,outv);
}

/**
 * \warning : values is provided as a square complete matrix
 */
void MatrixSquareSymmetric::_setValues(const double* values, bool /*byCol*/)
{
  // Check that the input argument corresponds to a square symmetric matrix
  for (int icol = 0; icol < getNCols(); icol++)
    for (int irow = 0; irow < getNRows(); irow++)
    {
      double val1 = values[icol * getNRows() + irow];
      double val2 = values[irow * getNCols() + icol];
      if (ABS(val1 - val2) > EPSILON10)
      {
        messerr("Argument 'values' must correspond to a Square Symmetric Matrix");
        messerr("- Element[%d,%d] = %lf",icol,irow,val1);
        messerr("- Element(%d,%d) = %lf",irow,icol,val2);
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

int MatrixSquareSymmetric::_invert()
{
  return matrix_invert_triangle(getNRows(),_squareSymMatrix.data(), -1);
}

double& MatrixSquareSymmetric::_getValueRef(int irow, int icol)
{
  _isIndexValid(irow,icol);
  int rank = _getIndexToRank(irow,icol);
  return _squareSymMatrix[rank];
}

void MatrixSquareSymmetric::_deallocate()
{
}

void MatrixSquareSymmetric::_recopy(const MatrixSquareSymmetric &r)
{
  _squareSymMatrix = r._squareSymMatrix;
}

void MatrixSquareSymmetric::_allocate()
{
  _squareSymMatrix.resize(_getMatrixSize());
}

int MatrixSquareSymmetric::_getIndexToRank(int irow, int icol) const
{
  int rank;

  int n = getNRows();
  if (irow >= icol)
    rank = icol * n + irow - icol * (icol + 1) / 2;
  else
    rank = irow * n + icol - irow * (irow + 1) / 2;
  return rank;
}

int MatrixSquareSymmetric::_getMatrixSize() const
{
  int n = getNRows();
  int size = n * (n + 1) / 2;
  return size;
}

int MatrixSquareSymmetric::_solve(const VectorDouble& b, VectorDouble& x) const
{
  int pivot;
  return matrix_solve(1,_squareSymMatrix.data(),b.data(),x.data(),
                      static_cast<int> (b.size()),1,&pivot);
}

String MatrixSquareSymmetric::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

   if (isSparse())
   {
     sstr << AMatrix::toString(strfmt);
   }
   else
   {
     sstr << "- Number of rows    = " <<  getNRows() << std::endl;
     sstr << "- Number of columns = " <<  getNCols() << std::endl;
     sstr << toMatrixSymmetric(String(), VectorString(), VectorString(),
                               true, getNCols(), getValues());
   }
  return sstr.str();
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

