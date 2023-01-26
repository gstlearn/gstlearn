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
#include <Matrix/AMatrixSquare.hpp>
#include <Matrix/MatrixFactory.hpp>
#include <Matrix/AMatrix.hpp>
#include "Basic/AException.hpp"
#include <math.h>

AMatrixSquare::AMatrixSquare(int nrow, bool sparse)
  : AMatrix(nrow, nrow, sparse)
{
}

AMatrixSquare::AMatrixSquare(const AMatrixSquare &r)
  : AMatrix(r)
{
}

AMatrixSquare& AMatrixSquare::operator= (const AMatrixSquare &r)
{
  if (this != &r)
  {
    AMatrix::operator=(r);
  }
  return *this;
}

AMatrixSquare::~AMatrixSquare()
{
}

void AMatrixSquare::_setNSize(int nval)
{
  _setNRows(nval);
  _setNCols(nval);
}

/**
 * Perform the product: this = t(Y) %*% X %*% Y
 * @param x: Square matrix
 * @param y: Matrix (possibly rectangular)
 * \remarks The number of rows of Y must be equal to the dimension of X
 * \remarks The output matrix is square with dimension equal to the number of columns of Y
 */
void AMatrixSquare::normMatrix(const AMatrixSquare& x, const AMatrix& y)
{
  if (x.getNSize() != y.getNRows())
  {
    my_throw("Incompatible matrix dimensions");
  }

  int n    = x.getNSize();
  int nout = y.getNCols();
  for (int irow = 0; irow < nout; irow++)
  {
    for (int icol = 0; icol < nout; icol++)
    {
      double value = 0.;
      for (int k = 0; k < n; k++)
      {
        for (int l = 0; l < n; l++)
        {
          value += y.getValue(k,irow) * x.getValue(k,l) * y.getValue(l,icol);
        }
      }
      setValue(irow,icol,value);
    }
  }
}


double AMatrixSquare::trace() const
{
  double res = 0.;
  for (int irow = 0; irow < getNSize(); irow++)
    res += getValue(irow,irow);
  return res;

}
/**
 * Perform the product: this = Y %*% X %*% t(Y)
 * @param x: Square matrix
 * @param y: Matrix (possibly rectangular)
 * \remarks The number of columns of Y must be equal to the dimension of X
 * \remarks The output matrix is square with dimension equal to the number of rows of Y
 */
void AMatrixSquare::normTMatrix(const AMatrixSquare& x, const AMatrix& y)
{
  if (x.getNSize() != y.getNCols())
  {
    my_throw("Incompatible matrix dimensions");
  }

  int n    = x.getNSize();
  int nout = y.getNRows();
  for (int irow = 0; irow < nout; irow++)
  {
    for (int icol = 0; icol < nout; icol++)
    {
      double value = 0.;
      for (int k = 0; k < n; k++)
      {
        for (int l = 0; l < n; l++)
        {
          value += y.getValue(irow,k) * x.getValue(k,l) * y.getValue(icol,l);
        }
      }
      setValue(irow,icol,value);
    }
  }
}

/**
 * Perform the product: this = t(R1) %*% X %*% R2 + t(R2) %*% X %*% R1
 * @param x: Square matrix
 * @param r1: Left Hand Matrix
 * @param r2: Right Hand Matrix
 * \remarks The number of rows of Y must be equal to the dimension of X
 * \remarks The output matrix is square with dimension equal to the number of columns of Y
 */
void AMatrixSquare::innerMatrix(const AMatrixSquare& x,
                                 const AMatrix& r1,
                                 const AMatrix& r2)
{
  int n = x.getNSize();
  if (n != r1.getNRows())
  {
    my_throw("Incompatible matrix dimensions");
  }
  if (n != r2.getNRows())
  {
    my_throw("Incompatible matrix dimensions");
  }


  for (int irow = 0; irow < n; irow++)
  {
    for (int icol = 0; icol < n; icol++)
    {
      double value = 0.;
      for (int k = 0; k < n; k++)
      {
        for (int l = 0; l < n; l++)
        {
          value += r1.getValue(k,irow) * x.getValue(k,l) * r2.getValue(l,icol);
          value += r2.getValue(k,irow) * x.getValue(k,l) * r1.getValue(l,icol);
        }
      }
      setValue(irow,icol,value);
    }
  }
}

bool AMatrixSquare::_isNumberValid(int nrows, int ncols) const
{
  AMatrix::_isNumbersValid(nrows,ncols);
  if (nrows != ncols)
  {
    messerr("Arguments 'nrows' and 'ncols' should be equal for Square Matrices");
    return false;
  }
  return true;
}

/*! Multiply the diagonal by a vector */
void AMatrixSquare::prodDiagByVector(const VectorDouble& diag)
{
  if ((int) diag.size() != getNRows())
    my_throw("Argument 'Diag' must match Matrix dimension");

  for (int i = 0; i < getNRows(); i++)
  {
    setValue(i, i, getValue(i, i) * diag[i]);
  }
}

/*! Divide the diagonal by a vector */
void AMatrixSquare::divideDiagByVector(const VectorDouble& diag)
{
  if ((int) diag.size() != getNRows())
    my_throw("Argument 'Diag' must match Matrix dimension");

  for (int i = 0; i < getNRows(); i++)
  {
    if (ABS(diag[i]) < EPSILON10)
      my_throw("Argument 'Diag' may not have too small values");
    setValue(i, i, getValue(i, i) / diag[i]);
  }
}

double AMatrixSquare::_determinant(void) const
{
  int neq = getNRows();

  /* Dispatch according to the matrix dimension */

  double deter = 0.;
  switch (neq)
  {
    case 1:
      deter = getValue(0, 0);
      break;

    case 2:
      deter = (getValue(0, 0) * getValue(1, 1)
               - getValue(1, 0) * getValue(0, 1));
      break;

    case 3:
      deter = (getValue(0, 0) * getValue(1, 1) * getValue(2, 2)
               + getValue(0, 1) * getValue(1, 2) * getValue(2, 0)
               + getValue(1, 0) * getValue(2, 1) * getValue(0, 2)
               - getValue(2, 0) * getValue(1, 1) * getValue(0, 2)
               - getValue(1, 0) * getValue(0, 1) * getValue(2, 2)
               - getValue(2, 1) * getValue(1, 2) * getValue(0, 0));
      break;

    default:

      int neqm1 = neq - 1;
      AMatrixSquare* c = MatrixFactory::createMatrixSquare(this,neqm1);

      for (int j1 = 0; j1 < neq; j1++)
      {
        for (int i = 1; i < neq; i++)
        {
          int j2 = 0;
          for (int j = 0; j < neq; j++)
          {
            if (j == j1) continue;
            c->setValue(i - 1, j2, getValue(i, j));
            j2++;
          }
        }
        deter += pow(-1.0, j1 + 2.0) * getValue(0, j1)
                 * c->_determinant();
      }
      break;
  }

  return(deter);
}
