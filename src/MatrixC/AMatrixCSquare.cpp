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
#include <MatrixC/AMatrixCSquare.hpp>
#include <MatrixC/AMatrixC.hpp>
#include "Basic/AException.hpp"

AMatrixCSquare::AMatrixCSquare(int nrow, bool sparse)
  : AMatrixC(nrow, nrow, sparse)
{
}

AMatrixCSquare::AMatrixCSquare(const AMatrixCSquare &r)
  : AMatrixC(r)
{
}

AMatrixCSquare& AMatrixCSquare::operator= (const AMatrixCSquare &r)
{
  if (this != &r)
  {
    AMatrixC::operator=(r);
  }
  return *this;
}

AMatrixCSquare::~AMatrixCSquare()
{
}

void AMatrixCSquare::_setNSize(int nval)
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
void AMatrixCSquare::normMatrix(const AMatrixCSquare& x, const AMatrixC& y)
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

/**
 * Perform the product: this = t(R1) %*% X %*% R2 + t(R2) %*% X %*% R1
 * @param x: Square matrix
 * @param r1: Left Hand Matrix
 * @param r2: Right Hand Matrix
 * \remarks The number of rows of Y must be equal to the dimension of X
 * \remarks The output matrix is square with dimension equal to the number of columns of Y
 */
void AMatrixCSquare::innerMatrix(const AMatrixCSquare& x,
                                 const AMatrixC& r1,
                                 const AMatrixC& r2)
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

bool AMatrixCSquare::_isNumberValid(int nrows, int ncols) const
{
  AMatrixC::_isNumbersValid(nrows,ncols);
  if (nrows != ncols)
  {
    messerr("Arguments 'nrows' and 'ncols' should be equal for Square Matrices");
    return false;
  }
  return true;
}

/*! Multiply the diagonal by a vector */
void AMatrixCSquare::prodDiagByVector(const VectorDouble& diag)
{
  if ((int) diag.size() != getNRows())
    my_throw("Argument 'Diag' must match Matrix dimension");

  for (int i = 0; i < getNRows(); i++)
  {
    setValue(i, i, getValue(i, i) * diag[i]);
  }
}

/*! Divide the diagonal by a vector */
void AMatrixCSquare::divideDiagByVector(const VectorDouble& diag)
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
