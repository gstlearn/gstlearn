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
#include "geoslib_enum.h"

#include "Enum/ECst.hpp"

#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/OptCst.hpp"

#include <math.h>
#include <cmath>

/*! \cond */
#define TRI(i)        (((i) * ((i) + 1)) / 2)
#define SQ(i,j,neq)   ((j) * neq + (i))
#define VECTOR(i,j)    vector[SQ(i,j,neq)]
#define A(i,j)         a[SQ(i,j,neq)]
#define ILS(i,j)       ils[SQ(i,j,neq)]
#define B(i,j)         b[SQ(i,j,neq)]
#define C(i,j)         c[SQ(i,j,neqm1)]
#define X(i,j)         x[SQ(i,j,neq)]
#define AMT(i,j)       amat[SQ(i,j,neq)]
#define HMT(i,j)       hmat[SQ(i,j,neq)]
#define AIMT(i,j)      aimat[SQ(i,j,neq)]
#define HA(i,j)        ha[SQ(i,j,neq)]
#define AI(i,j)        ai[SQ(i,j,neq)]
#define TEMP(i,j)      temp[SQ(i,j,na)]
#define AS(i,j)        a[SQ(i,j,neq)]
#define AT(i,j)        at[TRI(j)+(i)] /* for j >= i */
#define BS(i,j)        b[SQ(i,j,neq)] // Proposition a valider: c'etait j,i
#define XS(i,j)        x[SQ(i,j,neq)] // Proposition a valider: c'etait j,i
#define V1(i,j)        v1[SQ(i,j,n1)]
#define V2(i,j)        v2[SQ(i,j,n2)]
#define V3(i,j)        v3[SQ(i,j,n1)] /* Warning not to change the last argument: major bug */
#define V4(i,j)        v4[SQ(i,j,n1)]
#define TU(i,j)        tu[TRI(i) + (j)]       /* for i >= j */
#define TL(i,j)        tl[SQ(i,j,neq)-TRI(j)] /* for i >= j */
#define XL(i,j)        xl[SQ(i,j,neq)-TRI(j)] /* for i >= j */
#define AL(i,j)        al[SQ(i,j,neq)-TRI(j)] /* for i >= j */
#define TABEMT(i,j)    tabemat[SQ(i,j,neq)]
#define TABIMT(i,j)    tabimat[SQ(i,j,neq)]
#define TABOUT(i,j)    tabout[SQ(i,j,neq)]
#define EIGVEC(i,j)    eigvec[SQ(j,i,neq)]
#define A2(i,j)        a2[SQ(i,j,2 * neq)]
#define U(i,j)         u[SQ(j,i,neq)]
#define V(i,j)         v[SQ(j,i,neq)]
/*! \endcond */

static double _getTolInvert()
{
  return 1.e-20;
}
static double _getTolInvGen()
{
  return 1.e-20;
}
static double _getEpsMatrix()
{
  return 2.3e-16;
}

/*****************************************************************************/
/*!
 **  Calculates the eigen values and eigen vectors of a symmetric square matrix
 **             A X = X l
 **
 ** \return  Return code:
 ** \return   0 no error
 ** \return   1 convergence problem
 **
 ** \param[in]  a_in  square symmetric matrix (dimension = neq*neq)
 ** \param[in]  neq   matrix dimension
 **
 ** \param[out] value  matrix of the eigen values (dimension: neq)
 ** \param[out] vector matrix of the eigen vectors (dimension: neq*neq)
 **
 ** \remark  a_in is protected
 **
 *****************************************************************************/
int matrix_eigen(const double *a_in, int neq, double *value, double *vector)

{
  double a11, a12, a13, a21, a22, a23, a33, a34;
  double bb, cc, co, hold, s, si, v1, v2, bigk, qj, pp, temp;
  double *a, *work[4], *tmp;
  int *ind, i, j, k, ji, ki, kj, n1, n2, i1, i2, iter, error;

  /* Initializations */

  error = 1;
  ind   = nullptr;
  a = tmp = nullptr;
  for (i = 0; i < 4; i++)
    work[i] = nullptr;

  a34 = 0.0;
  if (neq == 1)
  {
    value[0] = a_in[0];
    VECTOR(0,0)= 1;
    return (0);
  }

  a = (double*) mem_alloc(sizeof(double) * neq * neq, 1);
  for (i = 0; i < neq * neq; i++)
    a[i] = a_in[i];

  for (i = 0; i < 4; i++)
  {
    work[i] = (double*) mem_alloc(sizeof(double) * neq, 1);
    for (j = 0; j < neq; j++)
      work[i][j] = 0.;
  }

  work[0][0] = AS(0, 0);
  if (neq <= 2)
  {
    work[0][1] = AS(1, 1);
    work[1][1] = AS(1, 0);
  }
  else
  {
    for (j = 1; j < neq; j++)
    {
      work[0][j] = AS(j, j);
      for (i = 0; i < j; i++)
        AS(i,j)= AS(j,i);
      }

      for (i=0; i<neq-2; i++)
      {
        pp = 0.;
        i1 = i+1;
        for (j=i1; j<neq; j++) pp += AS(i,j)*AS(i,j);
        work[1][i1] = SIGN (AS(i,i1),-sqrt(pp));
        if (pp <= 0) continue;
        hold = pp - work[1][i1]*AS(i,i1);
        AS(i,i1) -= work[1][i1];
        for (ki=i1; ki<neq; ki++)
        {
          qj = 0.;
          for (kj=i1; kj<=ki; kj++)
          qj += AS(kj,ki) * AS(i,kj);
          for (kj=ki+1; kj<neq; kj++)
          qj += AS(ki,kj) * AS(i,kj);
          work[2][ki] = qj/hold;
        }
        bigk = 0.;
        for (kj=i1; kj<neq; kj++)
        bigk += AS(i,kj) * work[2][kj];
        bigk /= 2.0 * hold;
        for (kj=i1; kj<neq; kj++)
        work[2][kj] -= bigk * AS(i,kj);
        for (ki=i1; ki<neq; ki++)
        for (kj=ki; kj<neq; kj++)
        AS(ki,kj) -=
        work[2][ki] * AS(i,kj) + work[2][kj] * AS(i,ki);
      }

      for (i=1; i<neq; i++)
      {
        hold = work[0][i];
        work[0][i] = AS(i,i);
        AS(i,i) = hold;
      }
      work[1][neq-1] = AS(neq-2,neq-1);
    }

  work[3][0] = work[0][0];
  for (i = 1; i < neq; i++)
    for (j = 0; j < 2; j++)
      work[3 - j][i] = work[j][i];

  iter = 0;
  do
  {
    n1 = n2 = 0;
    for (i2 = neq - 1; i2 > 0 && n2 == 0; i2--)
    {
      n1 = 0;
      for (i1 = i2; i1 > 0 && n1 == 0; i1--)
        if (ABS(work[2][i1]) <= _getEpsMatrix() * neq
                                * (ABS(work[3][i1-1]) + ABS(work[3][i1])))
          n1 = i1;
      if (n1 != i2) n2 = i2;
    }
    if (n2 < 1) break;

    bb = (work[3][n2] - work[3][n2 - 1]) / 2.;
    cc = work[2][n2] * work[2][n2];
    a22 = work[3][n1];
    a12 = a22 - work[3][n2];
    if (! isZero(bb) || ! isZero(cc)) a12 -= cc / (bb + SIGN(bb, sqrt(bb * bb + cc)));
    a23 = work[2][n1 + 1];
    a13 = a23;
    for (i = n1; i < n2; i++)
    {
      a33 = work[3][i + 1];
      if (i != n2 - 1) a34 = work[2][i + 2];
      s = sqrt(a12 * a12 + a13 * a13);
      si = a13 / s;
      co = a12 / s;
      if (i != n1) work[2][i] = s;
      a11 = co * a22 + si * a23;
      a12 = co * a23 + si * a33;
      a13 = si * a34;
      a21 = co * a23 - si * a22;
      a22 = co * a33 - si * a23;
      a23 = co * a34;
      work[3][i] = a11 * co + a12 * si;
      a12 = -a11 * si + a12 * co;
      work[2][i + 1] = a12;
      a22 = a22 * co - a21 * si;
    }
    work[3][n2] = a22;
    iter++;
  }
  while (iter < 10 * neq && n2 != -1);

  for (i = 0; i < neq; i++)
  {
    value[i] = work[0][i];
    work[2][i] = work[1][i];
    for (j = 0; j < neq; j++)
      VECTOR(j,i)= 0.;
    VECTOR(i,i)= 1.;
  }

  k = 0;
  for (n2 = neq - 1; n2 >= 1; n2--)
  {
    hold = work[3][n2];
    iter = 0;
    do
    {
      bb = (value[n2] - value[n2 - 1]) / 2.;
      cc = work[2][n2] * work[2][n2];
      a22 = value[n2];
      if (! isZero(bb) || ! isZero(cc)) a22 += cc / (bb + SIGN(bb, sqrt(bb * bb + cc)));
      for (i = 0; i < n2; i++)
        if (ABS(hold-a22) > ABS(work[3][i] - a22))
        {
          hold = work[3][i];
          work[3][i] = work[3][n2];
          work[3][n2] = hold;
        }

      n1 = 0;
      for (i1 = n2; i1 > 0 && n1 == 0; i1--)
        if (ABS(work[2][i1]) <= _getEpsMatrix() * neq
                                * (ABS(value[i1-1]) + ABS(value[i1]))) n1 = i1;
      if (n2 == n1) break;

      if (iter >= 3) hold = a22;
      k++;
      a22 = value[n1];
      a12 = a22 - hold;
      a23 = work[2][n1 + 1];
      a13 = a23;
      for (i = n1; i < n2; i++)
      {
        a33 = value[i + 1];
        if (i != n2 - 1) a34 = work[2][i + 2];
        s = SIGN(a12, sqrt(a12 * a12 + a13 * a13));
        si = a13 / s;
        co = a12 / s;
        for (ji = 0; ji <= MIN(neq - 1, i + k); ji++)
        {
          v1 = VECTOR(ji, i);
          v2 = VECTOR(ji, i + 1);
          VECTOR(ji,i)= v1*co+v2*si;
          VECTOR(ji,i+1)= v2*co-v1*si;
        }
        if (i != n1) work[2][i] = s;
        a11 = co * a22 + si * a23;
        a12 = co * a23 + si * a33;
        a13 = si * a34;
        a21 = co * a23 - si * a22;
        a22 = co * a33 - si * a23;
        a23 = co * a34;
        value[i] = a11 * co + a12 * si;
        a12 = -a11 * si + a12 * co;
        work[2][i + 1] = a12;
        a22 = a22 * co - a21 * si;
      }
      value[n2] = a22;
      iter++;
    }
    while (iter < 20 && n2 != n1);
    if (iter == 20) goto label_end;
  }

  for (j = 0; j < neq; j++)
  {
    v2 = VECTOR(0,j)*VECTOR(0,j);
    v1 = v2 * work[0][0];
    for (i=1; i<neq; i++)
    {
      v2 += VECTOR(i,j)*VECTOR(i,j);
      v1 += VECTOR(i,j)*
      (2.*work[1][i]*VECTOR(i-1,j)+work[0][i]*VECTOR(i,j));
    }
    value[j] = v1/v2;
  }

  if (neq > 2) for (j = 0; j < neq; j++)
    for (i = neq - 2; i > 0; i--)
      if (! isZero(work[1][i]))
      {
        pp = 0.;
        for (ki = i; ki < neq; ki++)
          pp += AS(i-1,ki)* VECTOR(ki,j);
        pp /= (AS(i-1,i)*work[1][i]);
        for (ki = i; ki < neq; ki++)
          VECTOR(ki,j)+= pp * AS(i-1,ki);
        }

        /* Sort the eigen values and the corresponding vectors */

  ind = (int*) mem_alloc(sizeof(int) * neq, 1);
  tmp = (double*) mem_alloc(sizeof(double) * neq * neq, 1);
  for (i = 0; i < neq; i++)
    ind[i] = i;
  ut_sort_double(0, neq, ind, value);
  for (i = 0; i < neq; i++)
    for (j = 0; j < neq; j++)
      tmp[i + neq * j] = VECTOR(i, ind[j]);
  for (i = 0; i < neq * neq; i++)
    vector[i] = tmp[i];

  /* Sorting in decreasing order */

  for (i = 0; i < neq / 2; i++)
  {
    k = neq - i - 1;
    temp = value[i];
    value[i] = value[k];
    value[k] = temp;
    for (j = 0; j < neq; j++)
    {
      temp = VECTOR(j, i);
      VECTOR(j,i)= VECTOR(j,k);
      VECTOR(j,k)= temp;
    }
  }

  /* Set the error returned code */

  error = 0;

  label_end: a = (double*) mem_free((char* ) a);
  ind = (int*) mem_free((char* ) ind);
  tmp = (double*) mem_free((char* ) tmp);
  for (i = 0; i < 4; i++)
    work[i] = (double*) mem_free((char* ) work[i]);

  if (error)
    print_matrix("Eigen matrix", 0, 1, neq, neq, NULL, a_in);

  return (error);
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
void matrix_product_safe(int n1,
                         int n2,
                         int n3,
                         const double *v1,
                         const double *v2,
                         double *v3)
{
  int i1, i2, i3, i4;

  if (v1 == v3 || v2 == v3)
    messageAbort("Violated protection in matrix_product_safe");

  for (i4 = 0; i4 < n1 * n3; i4++)
    v3[i4] = 0.;

  for (i3 = 0; i3 < n3; i3++)
    for (i1 = 0; i1 < n1; i1++)
      for (i2 = 0; i2 < n2; i2++)
        V3(i1,i3) += V1(i1,i2) * V2(i2,i3);
  return;
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
int matrix_prod_norme(int transpose,
                      int n1,
                      int n2,
                      const double *v1,
                      const double *a,
                      double *w)
{
  int i1, j1, i2, j2, ecr, neq;
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
            if (! isZero(vi)) for (j1 = 0; j1 < n1; j1++)
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
            if (! isZero(vi)) for (j2 = 0; j2 < n2; j2++)
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
void matrix_transpose(int n1, int n2, VectorDouble& v1, VectorDouble& w1)
{
  int ecr = 0;
  for (int i1 = 0; i1 < n1; i1++)
    for (int i2 = 0; i2 < n2; i2++)
      w1[ecr++] = V1(i1, i2);
  return;
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
int matrix_invert(double *a, int neq, int rank)
{
  for (int k = 0; k < neq; k++)
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

    for (int i = 0; i < neq; i++)
      if (i != k) A(i,k)= -A(i,k) / biga;

    for (int i = 0; i < neq; i++)
    {
      double hold = A(i, k);
      if (i != k)
        for (int j = 0; j < neq; j++)
          if (j != k) A(i,j)+= hold * A(k,j);
    }

    for (int j = 0; j < neq; j++)
      if (j != k) A(k,j)/= biga;

    A(k,k)= 1. / biga;
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
double matrix_determinant(int neq, const VectorDouble& b)
{
  switch (neq)
  {
    case 1:
      return B(0,0);

    case 2:
      return (B(0,0)* B(1,1) - B(1,0) * B(0,1));

      case 3:
      return ((B(0,0) * B(1,1) * B(2,2)
          + B(0,1) * B(1,2) * B(2,0)
          + B(1,0) * B(2,1) * B(0,2)
          - B(2,0) * B(1,1) * B(0,2)
          - B(1,0) * B(0,1) * B(2,2)
          - B(2,1) * B(1,2) * B(0,0)));

    default:

      /* Core allocation */
      double deter = 0.;
      int neqm1 = neq - 1;
      VectorDouble c(neqm1 * neqm1,0.);

      for (int j1=0; j1<neq; j1++)
      {
        for (int i=1; i<neq; i++)
        {
          int j2 = 0;
          for (int j=0; j<neq; j++)
          {
            if (j == j1) continue;
            C(i-1,j2) = B(i,j);
            j2++;
          }
        }
        deter += pow(-1.0,j1+2.0) * B(0,j1) * matrix_determinant(neqm1,c);
      }
      return deter;
    }
  return TEST;
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
int matrix_cholesky_decompose(const double *a, double *tl, int neq)
{
  double prod;
  int ip, jp, kp;

  for (ip = 0; ip < neq; ip++)
    for (jp = 0; jp <= ip; jp++)
      TL(ip,jp) = AS(ip,jp);

  for (ip = 0; ip < neq; ip++)
  {
    prod = TL(ip, ip);
    for (kp = 0; kp < ip; kp++)
      prod -= TL(ip,kp) * TL(ip,kp);
    if (prod < 0.) return (ip + 1);
    TL(ip,ip) = sqrt(prod);

    for (jp = ip + 1; jp < neq; jp++)
    {
      prod = TL(jp, ip);
      for (kp = 0; kp < ip; kp++)
        prod -= TL(ip,kp) * TL(jp,kp);
      if (TL(ip,ip)<= 0.) return(ip+1);
      TL(jp,ip) = prod / TL(ip,ip);
    }
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Performs the product between a triangular and a square matrix
 **  TL is the lower triangular matrix and X is a square matrix
 **
 ** \param[in]  mode Type of calculations:
 **             0 : X=TU%*%A
 **             1 : X=TL%*%A
 **             2 : X=A%*%TU
 **             3 : X=A%*%TL
 **             4 : X=t(A)%*%TU
 **             5 : X=t(A)%*%TL
 ** \param[in]  neq  number of equations in the system
 ** \param[in]  nrhs number of columns in x
 ** \param[in]  tl   Triangular matrix defined by column (dimension: neq * neq)
 ** \param[in]  a    matrix (dimension neq * nrhs)
 **
 ** \param[out] x    resulting matrix (dimension neq * nrhs)
 **
 *****************************************************************************/
void matrix_cholesky_product(int mode,
                             int neq,
                             int nrhs,
                             const double *tl,
                             const double *a,
                             double *x)
{
  int irhs, i, j, n1, n2;
  double val, *v2;
  const double *v1;

  if (mode == 0)
  {
    for (irhs = 0; irhs < nrhs; irhs++)
      for (i = 0; i < neq; i++)
      {
        val = 0.;
        for (j = i; j < neq; j++)
          val += TL(j,i) * AS(j,irhs);
        XS(i,irhs)= val;
      }
    }
    else if (mode == 1)
    {
      for (irhs=0; irhs<nrhs; irhs++)
      for (i=0; i<neq; i++)
      {
        val = 0.;
        for (j=0; j<=i; j++)
        val += TL(i,j) * AS(j,irhs);
        XS(i,irhs) = val;
      }
    }
    else if (mode == 2)
    {
      v2 = x;
      n2 = nrhs;
      for (irhs=0; irhs<nrhs; irhs++)
      for (i=0; i<neq; i++)
      {
        val = 0.;
        for (j=0; j<=i; j++)
        val += AS(irhs,j) * TL(i,j);
        V2(irhs,i) = val;
      }
    }
    else if (mode == 3)
    {
      v2 = x;
      n2 = nrhs;
      for (irhs=0; irhs<nrhs; irhs++)
      for (i=0; i<neq; i++)
      {
        val = 0.;
        for (j=i; j<neq; j++)
        val += AS(irhs,j) * TL(j,i);
        V2(irhs,i) = val;
      }
    }
    else if (mode == 4)
    {
      v1 = a;
      n1 = nrhs;
      v2 = x;
      n2 = nrhs;
      for (irhs=0; irhs<nrhs; irhs++)
      for (i=0; i<neq; i++)
      {
        val = 0.;
        for (j=0; j<=i; j++)
        val += V1(irhs,j) * TL(i,j);
        V2(irhs,i) = val;
      }
    }
    else if (mode == 5)
    {
      v1 = a;
      n1 = nrhs;
      v2 = x;
      n2 = nrhs;
      for (irhs=0; irhs<nrhs; irhs++)
      for (i=0; i<neq; i++)
      {
        val = 0.;
        for (j=i; j<neq; j++)
        val += V1(irhs,j) * TL(j,i);
        V2(irhs,i) = val;
      }
    }
  return;
}

/*****************************************************************************/
/*!
 **  Invert the Cholesky matrix
 **
 ** \param[in]  neq  number of equations in the system
 ** \param[in]  tl   lower triangular matrix defined by column
 **
 ** \param[out] xl   lower triangular inverted matrix defined by column
 **
 *****************************************************************************/
void matrix_cholesky_invert(int neq, const double *tl, double *xl)
{
  for (int i = 0; i < neq; i++)
  {
    for (int j = 0; j < i; j++)
    {
      double sum = 0.;
      for (int l = j; l < i; l++)
        sum += TL(i,l) * XL(l,j);
      XL(i,j) = - sum / TL(i,i);
    }
    XL(i,i) = 1. / TL(i,i);
  }
}

/****************************************************************************/
/*!
 **  Calculate the generalized inverse of a square symmetric matrix
 **
 ** \return  Error returned code
 **
 ** \param[in]  a         Symmetric matrix to be inverted
 ** \param[in]  neq       Number of equations
 **
 ** \param[out] tabout    Inverted matrix
 ** \param[out] cond      Condition number (MAX(ABS(eigval))/MIN(ABS(eigval)))
 **
 ** \remark The input and output matrices can matchprint
 **
 *****************************************************************************/
int matrix_invgen(double *a, int neq, double *tabout, double *cond)
{
  double *eigvec, *eigval, value, valcond;
  int i, j, k, error;

  /* Initializations */

  error = 1;
  eigvec = eigval = nullptr;

  /* Core allocation */

  eigval = (double*) mem_alloc(sizeof(double) * neq, 0);
  if (eigval == nullptr) goto label_end;
  eigvec = (double*) mem_alloc(sizeof(double) * neq * neq, 0);
  if (eigvec == nullptr) goto label_end;

  /* Calculate the eigen vectors */

  if (matrix_eigen(a, neq, eigval, eigvec)) goto label_end;
  valcond = MAX(ABS(eigval[0]), ABS(eigval[neq-1]));
  if (cond != nullptr) *cond = valcond;

  /* Calculate the generalized inverse */

  for (i = 0; i < neq; i++)
    for (j = 0; j < neq; j++)
    {
      value = 0.;
      for (k = 0; k < neq; k++)
      {
        if (ABS(eigval[k]) > valcond * _getTolInvGen()) value += EIGVEC(k,i)* EIGVEC(k,j) / eigval[k];
      }
      TABOUT(i,j) = value;
    }

    /* Set the error returned code */

  error = 0;

  label_end: eigval = (double*) mem_free((char* ) eigval);
  eigvec = (double*) mem_free((char* ) eigvec);
  return (error);
}

/*****************************************************************************/
/*!
 **  Perform a linear combinaison of matrices or vectors
 **            [C] = coeffa * [A] + coeffb * [B]
 **
 ** \param[in]  nval    Number of elements of the matrices or vectors
 ** \param[in]  coeffa  Coefficient applied to the first matrix or vector
 ** \param[in]  a       First matrix or vector (not used if NULL)
 ** \param[in]  coeffb  Coefficient applied to the second matrix or vector
 ** \param[in]  b       Second matrix or vector (not used if NULL)
 **
 ** \param[out] c       Resulting matrix or vector
 **
 ** \remark  No test is performed to check that the three matrices or vectors
 ** \remark  have the same dimensions
 ** \remark  Matrix c[] can coincide with matrices a[] or b[]
 **
 *****************************************************************************/
void matrix_combine(int nval,
                    double coeffa,
                    double *a,
                    double coeffb,
                    double *b,
                    double *c)
{
  int i;
  double value;

  for (i = 0; i < nval; i++)
  {
    value = 0.;
    if (a != nullptr) value += coeffa * a[i];
    if (b != nullptr) value += coeffb * b[i];
    c[i] = value;
  }
}

/*****************************************************************************/
/*!
 **  Transform a tridiagonal non-symmetric matrix into a symmetric one
 **
 ** \return  Error return code (see remarks)
 **
 ** \param[in]  vecdiag Vector for the main diagonal
 ** \param[in]  vecinf  Vector for the subdiagonal (in the last neq-1 positions)
 ** \param[in]  vecsup  Vector for the superdiagonal (in the first neq positions)
 ** \param[in]  neq    matrix dimension
 **
 ** \param[out] eigvec square symmetric matrix (dimension: neq * neq)
 ** \param[out] eigval vector (dimensionL neq)
 **
 ** \remark Given the nonsymmetric tridiagonal matrix, we must have the products
 ** \remark of corresponding pairs of off-diagonal elements are all non-negative,
 ** \remark and zero only when both factors are zero
 **
 *****************************************************************************/
int matrix_eigen_tridiagonal(const double *vecdiag,
                             const double *vecinf,
                             const double *vecsup,
                             int neq,
                             double *eigvec,
                             double *eigval)
{
  double *b, *e, h;
  int i, j;

  /* Initializations */

  e = (double*) mem_alloc(sizeof(double) * neq, 1);
  b = (double*) mem_alloc(sizeof(double) * neq * neq, 1);

  for (i = 1; i < neq; i++)
  {
    h = vecinf[i] * vecsup[i - 1];
    if (h < 0) return (1);
    if (isZero(h))
    {
      if (! isZero(vecinf[i]) || ! isZero(vecsup[i - 1])) return (2);
      e[i] = 0.;
    }
    else
    {
      e[i] = sqrt(h);
    }
  }

  for (i = 0; i < neq * neq; i++)
    b[i] = 0.;
  for (i = 0; i < neq; i++)
  {
    B(i,i)= vecdiag[i];
    if (i > 0) B(i,i-1) = B(i-1,i) = e[i];
  }

  /* Compute the eigen eigval and eigen eigvec */

  matrix_eigen(b, neq, eigval, eigvec);

  e[0] = 1.;
  for (i = 1; i < neq; i++)
  {
    if (! isZero(e[i]))
      e[i] *= e[i - 1] / vecsup[i - 1];
    else
      e[i] = 1.;
  }
  for (i = 0; i < neq; i++)
    for (j = 1; j < neq; j++)
      EIGVEC(i,j)*= e[j];

      /* Core deallocation */

  b = (double*) mem_free((char* ) b);
  e = (double*) mem_free((char* ) e);
  return (0);
}
