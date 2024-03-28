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
#include <Matrix/AMatrixSquare.hpp>
#include <Matrix/MatrixSquareGeneral.hpp>
#include <Matrix/MatrixFactory.hpp>
#include <Matrix/AMatrix.hpp>
#include "Basic/AException.hpp"
#include <math.h>

AMatrixSquare::AMatrixSquare(int nrow, int opt_eigen)
  : MatrixRectangular(nrow, nrow, opt_eigen)
{
}

AMatrixSquare::AMatrixSquare(const AMatrixSquare &r)
  : MatrixRectangular(r)
{
}

AMatrixSquare::AMatrixSquare(const AMatrix &m)
  : MatrixRectangular(m)
{
  if (!m.isSquare())
  {
    messerr("The input matrix should be Square");
    _clear();
    return;
  }
}

AMatrixSquare& AMatrixSquare::operator= (const AMatrixSquare &r)
{
  if (this != &r)
  {
    MatrixRectangular::operator=(r);
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

double AMatrixSquare::trace() const
{
  double res = 0.;
  for (int irow = 0; irow < getNSize(); irow++)
    res += getValue(irow,irow);
  return res;

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

double AMatrixSquare::determinant(void) const
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
      deter = (getValue(0, 0) * getValue(1, 1) - getValue(1, 0) * getValue(0, 1));
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
      MatrixSquareGeneral c(neqm1); // TODO Merge AMtrixSquare and MatrixSquareGeneral into MatrixSquare

      for (int j1 = 0; j1 < neq; j1++)
      {
        for (int i = 1; i < neq; i++)
        {
          int j2 = 0;
          for (int j = 0; j < neq; j++)
          {
            if (j == j1) continue;
            c.setValue(i - 1, j2, getValue(i, j));
            j2++;
          }
        }
        deter += pow(-1.0, j1 + 2.0) * getValue(0, j1) * c.determinant();
      }
      break;
  }

  return(deter);
}

double AMatrixSquare::normVec(const VectorDouble& vec)
{
  if (getNRows() != (int) vec.size())
  {
    messerr("Wrong dimension of 'vec' argument");
    return TEST;
  }
  double value = 0.;
  for (int irow = 0; irow < getNRows(); irow++)
    for (int icol = 0; icol < getNCols(); icol++)
      value += vec[irow] * getValue(irow, icol) * vec[icol];

  return value;
}

/*****************************************************************************/
/*!
 **  Performs the 'this' %*% diag(c) where c is a vector
 **
 ** \param[in]  mode  0: c as is; 1: sqrt(c); 2: 1/c; 3: 1/sqrt(c)
 ** \param[in]  c     vector
 **
 *****************************************************************************/
void AMatrixSquare::prodByDiagInPlace(int mode, const VectorDouble& c)
{
  int neq = getNRows();
  for (int i1 = 0; i1 < neq; i1++)
    for (int i2 = 0; i2 < neq; i2++)
    {
      double val = c[i2];
      if (mode == 1)
        val = sqrt(val);
      else if (mode == 2)
        val = 1. / val;
      else if (mode == 3) val = 1. / sqrt(val);
      setValue(i1, i2, getValue(i1, i2) * val);
    }
  return;
}

