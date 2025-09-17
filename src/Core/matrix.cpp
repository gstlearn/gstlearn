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
#include "Basic/Utilities.hpp"
#include "Basic/VectorNumT.hpp"
#include "geoslib_define.h"
#include "geoslib_old_f.h"
#include <cmath>

/*! \cond */
#define TRI(i)        (((i) * ((i) + 1)) / 2)
#define SQ(i, j, neq) ((j) * neq + (i))
#define A(i, j)       a[SQ(i, j, neq)]
#define B(i, j)       b[SQ(i, j, neq)]
#define C(i, j)       c[SQ(i, j, neqm1)]
#define X(i, j)       x[SQ(i, j, neq)]
#define AS(i, j)      a[SQ(i, j, neq)]
#define BS(i, j)      b[SQ(i, j, neq)] // Proposition a valider: c'etait j,i
#define XS(i, j)      x[SQ(i, j, neq)] // Proposition a valider: c'etait j,i
#define V1(i, j)      v1[SQ(i, j, n1)]
#define V2(i, j)      v2[SQ(i, j, n2)]
#define V3(i, j)      v3[SQ(i, j, n1)]           /* Warning not to change the last argument: major bug */
#define TU(i, j)      tu[TRI(i) + (j)]           /* for i >= j */
#define TL(i, j)      tl[SQ(i, j, neq) - TRI(j)] /* for i >= j */
#define U(i, j)       u[SQ(j, i, neq)]
#define V(i, j)       v[SQ(j, i, neq)]
/*! \endcond */

namespace gstlrn
{
static double _getTolInvert()
{
  return EPSILON25;
}

/*****************************************************************************/
/*!
 **  Performs the product of two matrices
 **
 ** \param[in]  n1 matrix dimension
 ** \param[in]  n2 matrix dimension
 ** \param[in]  n3 matrix dimension
 ** \param[in]  v1 rectangular matrix (n1,n2)
 ** \param[in]  v2 rectangular matrix (n2,n3)
 **
 ** \param[out] v3 rectangular matrix (n1,n3)
 **
 ** \remark  The matrix v3[] may NOT coincide with one of the two initial ones
 **
 *****************************************************************************/
void matrix_product_safe(Id n1, Id n2, Id n3, const double* v1, const double* v2, double* v3)
{
  Id i1, i2, i3, i4;

  if (v1 == v3 || v2 == v3)
    messageAbort("Violated protection in matrix_product_safe");

  for (i4 = 0; i4 < n1 * n3; i4++) v3[i4] = 0.;

  for (i3 = 0; i3 < n3; i3++)
    for (i1 = 0; i1 < n1; i1++)
      for (i2 = 0; i2 < n2; i2++) V3(i1, i3) += V1(i1, i2) * V2(i2, i3);
}

/*****************************************************************************/
/*!
 **  Performs the product t(G) %*% A %*% G or G %*% A %*% t(G)
 **
 ** \return Error return code
 **
 ** \param[in]  transpose transposition mode
 **                      -1 : transpose the first term
 **                      +1 : transpose the last term
 ** \param[in]  n1        matrix dimension
 ** \param[in]  n2        matrix dimension
 ** \param[in]  v1        rectangular matrix (n1,n2)
 ** \param[in]  a         square matrix (optional)
 **
 ** \param[out] w         square matrix
 **
 ** \remarks According to the value of 'transpose':
 ** \remarks -1: the output array has dimension (n2,n2)
 ** \remarks +1: the output array has dimension (n1,n1)
 ** \remarks According to the value of 'transpose':
 ** \remarks -1: the optional array A has dimension (n1,n1)
 ** \remarks +1: the optional array A has dimension (n2,n2)
 **
 *****************************************************************************/
Id matrix_prod_norme(Id transpose, Id n1, Id n2, const double* v1, const double* a, double* w)
{
  Id i1, j1, i2, j2, ecr, neq;
  double value, vala, vi;

  ecr = 0;
  switch (transpose)
  {
    case -1:
      neq = n1;
      for (i2 = 0; i2 < n2; i2++)
        for (j2 = 0; j2 < n2; j2++)
        {
          value = 0.;

          for (i1 = 0; i1 < n1; i1++)
          {
            vi = V1(i1, i2);
            if (!isZero(vi))
              for (j1 = 0; j1 < n1; j1++)
              {
                if (a != nullptr)
                  vala = AS(i1, j1);
                else
                  vala = (i1 == j1);
                value += vi * vala * V1(j1, j2);
              }
          }
          w[ecr++] = value;
        }
      break;

    case +1:
      neq = n2;
      for (i1 = 0; i1 < n1; i1++)
        for (j1 = 0; j1 < n1; j1++)
        {
          value = 0.;
          for (i2 = 0; i2 < n2; i2++)
          {
            vi = V1(i1, i2);
            if (!isZero(vi))
              for (j2 = 0; j2 < n2; j2++)
              {
                if (a != nullptr)
                  vala = AS(i2, j2);
                else
                  vala = (i2 == j2);
                value += vi * vala * V1(j1, j2);
              }
          }
          w[ecr++] = value;
        }
      break;

    default:
      break;
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Transpose a (square or rectangular) matrix
 **
 ** \param[in]  n1 matrix dimension
 ** \param[in]  n2 matrix dimension
 ** \param[in]  v1 rectangular matrix (n1,n2)
 **
 ** \param[out] w1 rectangular matrix (n2,n1)
 **
 ** \remark  The matrix w1[] may NOT coincide with v1[]
 **
 *****************************************************************************/
void matrix_transpose(Id n1, Id n2, VectorDouble& v1, VectorDouble& w1)
{
  Id ecr = 0;
  for (Id i1 = 0; i1 < n1; i1++)
    for (Id i2 = 0; i2 < n2; i2++)
      w1[ecr++] = V1(i1, i2);
}

/*****************************************************************************/
/*!
 **  Invert a symmetric square matrix
 **  Pivots are assumed to be located on the diagonal
 **
 ** \return  Return code: 0 no error; k if the k-th pivot is zero
 **
 ** \param[in,out] a input matrix, destroyed in computation and replaced by
 **                  resultant inverse
 ** \param[in]  neq  number of equations in the matrix 'a'
 ** \param[in]  rank Type of message when inversion problem is encountered
 **                  >=0: message involves 'rank+1'
 **                  -1:  neutral message
 **                  -2:  no message
 **
 ** \remark  It is unnecessary to edit a message if inversion problem occurs
 **
 *****************************************************************************/
Id matrix_invert(double* a, Id neq, Id rank)
{
  for (Id k = 0; k < neq; k++)
  {
    double biga = A(k, k);
    if (ABS(biga) < _getTolInvert())
    {
      if (rank >= 0)
        messerr("Error in matrix inversion (rank=%d) : Pivot #%d is null",
                rank + 1, k + 1);
      else if (rank == -1)
        messerr("Error in matrix inversion : Pivot #%d is null", k + 1);
      return (k + 1);
    }

    for (Id i = 0; i < neq; i++)
      if (i != k) A(i, k) = -A(i, k) / biga;

    for (Id i = 0; i < neq; i++)
    {
      double hold = A(i, k);
      if (i != k)
        for (Id j = 0; j < neq; j++)
          if (j != k) A(i, j) += hold * A(k, j);
    }

    for (Id j = 0; j < neq; j++)
      if (j != k) A(k, j) /= biga;

    A(k, k) = 1. / biga;
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Calculate the determinant of the square matrix (full storage)
 **
 ** \return  Value of the determinant
 **
 ** \param[in]  neq    Size of the matrix
 ** \param[in]  b      Square matrix to be checked
 **
 *****************************************************************************/
double matrix_determinant(Id neq, const VectorDouble& b)
{
  switch (neq)
  {
    case 1:
      return B(0, 0);

    case 2:
      return (B(0, 0) * B(1, 1) - B(1, 0) * B(0, 1));

    case 3:
      return ((B(0, 0) * B(1, 1) * B(2, 2) + B(0, 1) * B(1, 2) * B(2, 0) + B(1, 0) * B(2, 1) * B(0, 2) - B(2, 0) * B(1, 1) * B(0, 2) - B(1, 0) * B(0, 1) * B(2, 2) - B(2, 1) * B(1, 2) * B(0, 0)));

    default:
      /* Core allocation */
      double deter = 0.;
      Id neqm1     = neq - 1;
      VectorDouble c(neqm1 * neqm1, 0.);

      for (Id j1 = 0; j1 < neq; j1++)
      {
        for (Id i = 1; i < neq; i++)
        {
          Id j2 = 0;
          for (Id j = 0; j < neq; j++)
          {
            if (j == j1) continue;
            C(i - 1, j2) = B(i, j);
            j2++;
          }
        }
        deter += pow(-1.0, j1 + 2.0) * B(0, j1) * matrix_determinant(neqm1, c);
      }
      return deter;
  }
}

/*****************************************************************************/
/*!
 **  Performs the Cholesky triangular decomposition of a definite
 **  positive symmetric matrix
 **         A = t(TL) * TL
 **
 ** \return  Return code: >0 rank of zero pivot or 0 if no error
 **
 ** \param[in]  a   symmetric matrix
 ** \param[in]  neq number of equations in the system
 **
 ** \param[out] tl  Lower triangular matrix defined by column
 **
 ** \remark  the matrix a[] is destroyed during the calculations
 **
 *****************************************************************************/
Id matrix_cholesky_decompose(const double* a, double* tl, Id neq)
{
  double prod;
  Id ip, jp, kp;

  for (ip = 0; ip < neq; ip++)
    for (jp = 0; jp <= ip; jp++)
      TL(ip, jp) = AS(ip, jp);

  for (ip = 0; ip < neq; ip++)
  {
    prod = TL(ip, ip);
    for (kp = 0; kp < ip; kp++)
      prod -= TL(ip, kp) * TL(ip, kp);
    if (prod < 0.) return (ip + 1);
    TL(ip, ip) = sqrt(prod);

    for (jp = ip + 1; jp < neq; jp++)
    {
      prod = TL(jp, ip);
      for (kp = 0; kp < ip; kp++)
        prod -= TL(ip, kp) * TL(jp, kp);
      if (TL(ip, ip) <= 0.) return (ip + 1);
      TL(jp, ip) = prod / TL(ip, ip);
    }
  }
  return (0);
}

} // namespace gstlrn
