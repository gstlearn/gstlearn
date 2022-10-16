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
#include "geoslib_enum.h"

#include "Enum/ECst.hpp"

#include "Basic/Utilities.hpp"
#include "Basic/OptCst.hpp"

#include <math.h>
#include <cmath>

/*! \cond */
#define TRI(i)        (((i) * ((i) + 1)) / 2)
#define SQ(i,j,neq)   ((j) * neq + (i))
#define VECTOR(i,j)    vector[SQ(i,j,neq)]
#define A(i,j)         a[SQ(i,j,neq)]
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

static double Epsilon = 1.0e-06;

static int LNG_LHS = 0;
static int LNG_RHS = 0;
static double *LHS_TAB = NULL;
static double *RHS_TAB = NULL;

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
static double _getEpsSVD()
{
  return 1.0e-5;
}

/*****************************************************************************/
/*!
 **  Solve a system of linear equations with symmetric coefficient
 **  matrix upper triangular part of which is stored columnwise. Use
 **  is restricted to matrices whose leading principal minors have
 **  non-zero determinants.
 **
 ** \return  Return code:  0 no error
 ** \return  -1 if neq<1 +k when zero pivot encountered at k-th iteration
 **
 ** \param[in]  neq  number of equations in the system
 ** \param[in]  nrhs number of right-hand side vectors
 ** \param[in]  at   upper triangular matrix by row (dimension = neq*(neq+1)/2)
 ** \param[in]  b    right-hand side matrix (dimension = neq*nrhs)
 **
 ** \param[out] x: matrix of solutions (dimension = neq*nrhs)
 **
 ** \remark  1 - The algorithm is gauss elimination. pivots are taken along
 ** \remark      main diagonal. There are no interchanges and no search for
 ** \remark      maximal element. The equations remain in the same order as on
 ** \remark      input and the right-hand sides are kept.
 ** \remark      It is therefore possible to solve the system with the last
 ** REMAKRS:      equation removed by simply applying back substitution on the
 ** \remark      triangularized matrix 'a'.
 ** \remark  2 - Return code=k at the rank k indicates that the determinant of
 ** \remark      g(k), the leading principal minor of order k, is zero.
 ** \remark      Generally, after triangularization the diagonal term of the
 ** \remark      i-th row is : at*(i*(i+1)/2) = det(g(i)) / det(g(i-1))
 ** \remark      therefore use of this function is restricted to matrices such
 ** \remark      that det(g(i)) is never zero.
 ** \remark   3- The arrays at and b are modified by this function The
 ** \remark      arrays b and x may not coincide
 **
 *****************************************************************************/
static int st_matrix_solve(double *at, double *b, double *x, int neq, int nrhs)

{
  double pivot, ratio;
  int i, j, k, l;

  for (k = 0; k < neq - 1; k++)
  {
    pivot = AT(k, k);
    if (ABS(pivot) < _getTolInvert()) return (k + 1);
    for (i = k + 1; i < neq; i++)
    {
      ratio = AT(k,i)/ pivot;
      for (j=i; j<neq; j++) AT(i,j) -= AT(k,j) * ratio;
      for (l=0; l<nrhs; l++) BS(i,l) -= BS(k,l) * ratio;
    }
  }

  pivot = AT(neq - 1, neq - 1);
  if (ABS(pivot) < _getTolInvert()) return (neq);

  for (l = 0; l < nrhs; l++)
    XS(neq-1,l)= BS(neq-1,l) / pivot;

  for (l = 0; l < nrhs; l++)
  {
    for (k = neq - 2; k >= 0; k--)
    {
      ratio = BS(k, l);
      for (j = k + 1; j < neq; j++)
        ratio -= AT(k,j)* XS(j,l);
      XS(k,l)= ratio / AT(k,k);
    }
  }

  return (0);
}

/*****************************************************************************/
/*!
 **  Calculates the eigen values and eigen vectors of a symmetric square matrix
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
    if (bb != 0. || cc != 0.) a12 -= cc / (bb + SIGN(bb, sqrt(bb * bb + cc)));
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
      if (bb != 0. || cc != 0.) a22 += cc / (bb + SIGN(bb, sqrt(bb * bb + cc)));
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
      if (work[1][i] != 0)
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

  if (error) print_matrix("Eigen matrix", 0, 1, neq, neq, NULL, a_in);

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
 ** \remark  The matrix v3[] may coincide with one of the two initial ones
 **
 *****************************************************************************/
void matrix_product(int n1,
                    int n2,
                    int n3,
                    const double *v1,
                    const double *v2,
                    double *v3)
{
  double *v4;
  int i1, i2, i3, i4;

  v4 = (double*) mem_alloc(sizeof(double) * n1 * n3, 1);
  for (i4 = 0; i4 < n1 * n3; i4++)
    v4[i4] = 0.;

  for (i3 = 0; i3 < n3; i3++)
    for (i1 = 0; i1 < n1; i1++)
      for (i2 = 0; i2 < n2; i2++)
        V4(i1,i3)+= V1(i1,i2) * V2(i2,i3);

  for (i4 = 0; i4 < n1 * n3; i4++)
    v3[i4] = v4[i4];
  v4 = (double*) mem_free((char* ) v4);

  return;
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

  for (i4 = 0; i4 < n1 * n3; i4++)
    v3[i4] = 0.;

  for (i3 = 0; i3 < n3; i3++)
    for (i1 = 0; i1 < n1; i1++)
      for (i2 = 0; i2 < n2; i2++)
        V3(i1,i3)+= V1(i1,i2) * V2(i2,i3);
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
            if (vi != 0.) for (j1 = 0; j1 < n1; j1++)
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
            if (vi != 0.) for (j2 = 0; j2 < n2; j2++)
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
void matrix_transpose(int n1, int n2, double *v1, double *w1)
{
  int i1, i2, ecr;

  ecr = 0;
  for (i1 = 0; i1 < n1; i1++)
    for (i2 = 0; i2 < n2; i2++)
      w1[ecr++] = V1(i1, i2);

  return;
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
void matrix_int_transpose(int n1, int n2, int *v1, int *w1)
{
  int i1, i2, ecr;

  ecr = 0;
  for (i1 = 0; i1 < n1; i1++)
    for (i2 = 0; i2 < n2; i2++)
      w1[ecr++] = V1(i1, i2);

  return;
}


/*****************************************************************************/
/*!
 **  Transpose a matrix in place
 **
 ** \param[in]  n1 matrix dimension
 ** \param[in]  n2 matrix dimension
 ** \param[in,out]  v1 rectangular matrix (n1,n2)
 **
 *****************************************************************************/
void matrix_transpose_in_place(int n1, int n2, double *v1)
{
  int i1, i2, ecr;
  double *w1;

  w1 = (double*) mem_alloc(sizeof(double) * n1 * n2, 1);

  ecr = 0;
  for (i1 = 0; i1 < n1; i1++)
    for (i2 = 0; i2 < n2; i2++)
      w1[ecr++] = V1(i1, i2);

  for (int i = 0; i < n1 * n2; i++)
    v1[i] = w1[i];

  w1 = (double*) mem_free((char* ) w1);

  return;
}

/*****************************************************************************/
/*!
 **  Transpose an integer matrix in place
 **
 ** \param[in]  n1 matrix dimension
 ** \param[in]  n2 matrix dimension
 ** \param[in,out]  v1 rectangular matrix (n1,n2)
 **
 *****************************************************************************/
void matrix_int_transpose_in_place(int n1, int n2, int *v1)
{
  int i1, i2, ecr;
  int *w1;

  w1 = (int*) mem_alloc(sizeof(int) * n1 * n2, 1);

  ecr = 0;
  for (i1 = 0; i1 < n1; i1++)
    for (i2 = 0; i2 < n2; i2++)
      w1[ecr++] = V1(i1, i2);

  for (int i = 0; i < n1 * n2; i++)
    v1[i] = w1[i];

  w1 = (int*) mem_free((char* ) w1);

  return;
}

/*****************************************************************************/
/*!
 **  Invert a symmetric square matrix
 **  Pivots are assumed to be located on the diagonal
 **
 ** \return  Return code: 0 no error; k if the k-th pivot is zero
 **
 ** \param[in,out] a    input matrix, destroyed in computation and replaced by
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
  int i, j, k;
  double biga, hold;

  for (k = 0; k < neq; k++)
  {
    biga = A(k, k);
    if (ABS(biga) < _getTolInvert())
    {
      if (rank >= 0)
        messerr("Error in matrix inversion (rank=%d) : Pivot #%d is null",
                rank + 1, k + 1);
      else if (rank == -1)
        messerr("Error in matrix inversion : Pivot #%d is null", k + 1);
      return (k + 1);
    }

    for (i = 0; i < neq; i++)
      if (i != k) A(i,k)= -A(i,k) / biga;

    for (i = 0; i < neq; i++)
    {
      hold = A(i, k);
      if (i != k) for (j = 0; j < neq; j++)
        if (j != k) A(i,j)+= hold * A(k,j);
    }

    for (j = 0; j < neq; j++)
      if (j != k) A(k,j)/= biga;

    A(k,k)= 1. / biga;
  }

  return (0);
}

/*****************************************************************************/
/*!
 **  Invert a symmetric square matrix
 **  Pivots are assumed to located on the diagonal
 **
 ** \return  Return code: 0 no error; k if the k-th pivot is zero
 **
 ** \param[in]  a    input matrix
 ** \param[in]  neq  number of equations in the matrix 'a'
 **
 ** \param[out] b    output matrix
 **
 ** \remark  The difference with matrix_invert() is that the output
 ** \remark  matrix is different from input matrix
 **
 *****************************************************************************/
int matrix_invert_copy(const double *a, int neq, double *b)
{
  int i, error;

  /* Copy the input matrix into the output matrix */

  for (i = 0; i < neq * neq; i++)
    b[i] = a[i];

  /* Invert the matrix */

  error = matrix_invert(b, neq, 0);

  return (error);
}

/*****************************************************************************/
/*!
 **  Solve a system of linear equations with symmetric coefficient
 **  matrix upper triangular part of which is stored columnwise. Use
 **  is restricted to matrices whose leading principal minors have
 **  non-zero determinants.
 **
 ** \return  Error return code : 1 for core problem; 0 otherwise
 **
 ** \param[in]  mode Storage: 0 for upper triangular, 1 for lower triangular
 ** \param[in]  a    triangular matrix (dimension = neq*(neq+1)/2)
 ** \param[in]  b    right-hand side matrix (dimension = neq*nrhs)
 ** \param[in]  neq  number of equations in the system
 ** \param[in]  nrhs number of right-hand side vectors
 **
 ** \param[out] x     matrix of solutions (dimension = neq*nrhs)
 ** \param[out] pivot rank of the pivoting error (0 none)
 **
 *****************************************************************************/
int matrix_solve(int mode,
                 const double *a,
                 const double *b,
                 double *x,
                 int neq,
                 int nrhs,
                 int *pivot)
{
  int i, loc_lhs, loc_rhs, error;

  /* Calculate array dimensions */

  error = 1;
  loc_lhs = neq * (neq + 1) / 2;
  loc_rhs = neq * nrhs;

  /* Core management for the L.H.S.*/

  if (loc_lhs > LNG_LHS)
  {
    if (LNG_LHS > 0 && LHS_TAB != nullptr)
    {
      LHS_TAB = (double*) mem_free((char* ) LHS_TAB);
      LNG_LHS = 0;
    }
    LHS_TAB = (double*) mem_alloc(sizeof(double) * loc_lhs, 0);
    if (LHS_TAB == nullptr) goto label_end;
    LNG_LHS = loc_lhs;
  }

  /* Core management for R.H.S. */

  if (loc_rhs > LNG_RHS)
  {
    if (LNG_RHS > 0 && RHS_TAB != nullptr)
    {
      RHS_TAB = (double*) mem_free((char* ) RHS_TAB);
      LNG_RHS = 0;
    }
    RHS_TAB = (double*) mem_alloc(sizeof(double) * loc_rhs, 0);
    if (RHS_TAB == nullptr) goto label_end;
    LNG_RHS = loc_rhs;
  }

  /* Backup the arrays */

  if (mode == 0)
    for (i = 0; i < loc_lhs; i++)
      LHS_TAB[i] = a[i];
  else
    matrix_tl2tu(neq, a, LHS_TAB);
  for (i = 0; i < loc_rhs; i++)
    RHS_TAB[i] = b[i];

  *pivot = st_matrix_solve(LHS_TAB, RHS_TAB, x, neq, nrhs);

  error = 0;

  label_end: return (error);
}

/****************************************************************************/
/*!
 **  Check if a matrix is symmetric
 **
 ** \return  1 if the matrix is symmetric; 0 otherwise
 **
 ** \param[in]  neq      Size of the matrix
 ** \param[in]  a        Symmetric square matrix to be checked
 ** \param[in]  verbose  1 for the verbose option
 **
 *****************************************************************************/
int is_matrix_symmetric(int neq, const double *a, int verbose)
{
  int i, j;
  double ratio;

  for (i = 0; i < neq; i++)
    for (j = 0; j < neq; j++)
    {
      ratio = ABS(A(i,j) + A(j,i));
      if (ratio <= _getEpsMatrix()) ratio = 1.;
      if (ABS(A(i,j) - A(j,i)) / ratio > Epsilon)
      {
        if (verbose)
          messerr("The A[%d,%d]=%lf != A[%d,%d]=%lf", i, j, A(i, j), j, i,
                  A(j, i));
        return (0);
      }
    }
  return (1);
}

/****************************************************************************/
/*!
 **  Check if a matrix is definite positive
 **
 ** \return  1 if the matrix is definite positive; 0 otherwise
 **
 ** \param[in]  neq      Size of the matrix
 ** \param[in]  a        Symmetric square matrix to be checked
 ** \param[in]  verbose  1 for the verbose option
 **
 ** \param[out]  valpro Array of eigen values  (Dimension: neq)
 ** \param[out]  vecpro Array of eigen vectors (Dimension: neq * neq)
 **
 *****************************************************************************/
int is_matrix_definite_positive(int neq,
                                const double *a,
                                double *valpro,
                                double *vecpro,
                                int verbose)
{
  int i, code;

  code = 1;
  if (!is_matrix_symmetric(neq, a, verbose))
  {
    if (verbose) messerr("Matrix is not symmetric");
    return (0);
  }

  /* Verbose option: Print the matrix */

  if (verbose) print_matrix("Matrix to be checked", 0, 1, neq, neq, NULL, a);

  /* Calculate the eigen values and vectors */

  if (matrix_eigen(a, neq, valpro, vecpro)) messageAbort("matrix_eigen");

  /* Verbose option: Print the Eigen values and vectors */

  if (verbose)
  {
    print_matrix("Eigen Values", 0, 1, 1, neq, NULL, valpro);
    print_matrix("Eigen Vectors", 0, 1, neq, neq, NULL, vecpro);
  }

  /* Check if the eigen values are all positive */

  for (i = 0; i < neq; i++)
    if (valpro[i] < 0.)
    {
      if (valpro[i] < -1.0e-10)
      {
        if (verbose)
          messerr("The matrix is not definite positive: Eigen value #%d = %lf",
                  i + 1, valpro[i]);
        code = 0;
        goto label_end;
      }
      else
      {
        valpro[i] = 0.;
      }
    }

  label_end: return (code);
}

/****************************************************************************/
/*!
 **  Check if all the elements of a matrix are non-negative
 **
 ** \return  1 if the matrix is non-negative; 0 otherwise
 **
 ** \param[in]  nrow     Number of rows
 ** \param[in]  ncol     Number of columns
 ** \param[in]  a        Array to be checked
 ** \param[in]  verbose  1 for the verbose option
 **
 *****************************************************************************/
int is_matrix_non_negative(int nrow,
                                           int ncol,
                                           double *a,
                                           int verbose)
{
  int i;

  for (i = 0; i < nrow * ncol; i++)
  {
    if (a[i] < 0.)
    {
      if (verbose) messerr("The matrix is not non-negative");
      return (0);
    }
  }

  return (1);
}

/****************************************************************************/
/*!
 **  Check if all the elements of a matrix are null
 **
 ** \return  1 if the matrix is null; 0 otherwise
 **
 ** \param[in]  nrow     Number of rows
 ** \param[in]  ncol     Number of columns
 ** \param[in]  a        Array to be checked
 ** \param[in]  verbose  1 for the verbose option
 **
 *****************************************************************************/
int is_matrix_null(int nrow, int ncol, const double *a, int verbose)
{
  int i;

  for (i = 0; i < nrow * ncol; i++)
  {
    if (a[i] != 0.)
    {
      if (verbose) messerr("The matrix is not null");
      return (0);
    }
  }

  if (verbose) messerr("The matrix is null");
  return (1);
}

/****************************************************************************/
/*!
 **  Check if a matrix is a correlation matrix
 **
 ** \return  1 if the matrix is a correlation matrix; 0 otherwise
 **
 ** \param[in]  neq    Size of the matrix
 ** \param[in]  a      Symmetric square matrix to be checked
 **
 *****************************************************************************/
int is_matrix_correlation(int neq, double *a)
{
  double *valpro, *vecpro;
  int i, status;

  /* Initializations */

  status = 0;
  valpro = (double*) mem_alloc(sizeof(double) * neq, 1);
  vecpro = (double*) mem_alloc(sizeof(double) * neq * neq, 1);

  /* Check that the matrix is definite positive */

  if (!is_matrix_definite_positive(neq, a, valpro, vecpro, 1)) goto label_end;

  /* Check that the diagonal is 1 */

  for (i = 0; i < neq; i++)
    if (A(i,i)!= 1.) goto label_end;

    /* Set the status */

  status = 1;

  label_end:

  /* Core deallocation */

  valpro = (double*) mem_free((char* ) valpro);
  vecpro = (double*) mem_free((char* ) vecpro);
  return (status);
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
double matrix_determinant(int neq, const double *b)
{
  int i, j, neqm1, j1, j2;
  double *c, deter;

  /* Dispatch according to the matrix dimension */

  deter = 0.;
  switch (neq)
  {
    case 1:
      deter = B(0, 0);
      break;

    case 2:
      deter = B(0,0)* B(1,1) - B(1,0) * B(0,1);
      break;

      case 3:
      deter = ( B(0,0) * B(1,1) * B(2,2)
          + B(0,1) * B(1,2) * B(2,0)
          + B(1,0) * B(2,1) * B(0,2)
          - B(2,0) * B(1,1) * B(0,2)
          - B(1,0) * B(0,1) * B(2,2)
          - B(2,1) * B(1,2) * B(0,0));
      break;

      default:

      /* Core allocation */
      neqm1 = neq - 1;
      c = (double *) mem_alloc(sizeof(double) * neqm1 * neqm1,1);

      for (j1=0; j1<neq; j1++)
      {
        for (i=1; i<neq; i++)
        {
          j2 = 0;
          for (j=0; j<neq; j++)
          {
            if (j == j1)
            continue;
            C(i-1,j2) = B(i,j);
            j2++;
          }
        }
        deter += pow(-1.0,j1+2.0) * B(0,j1) * matrix_determinant(neqm1,c);
      }

      /* Core deallocation */
      c = (double *) mem_free((char *) c);
      break;
    }

  return (deter);
}

/****************************************************************************/
/*!
 **  Calculate the cofactor of the matrix
 **
 ** \return  Value of the determinant
 **
 ** \param[in]  neq    Size of the matrix
 ** \param[in]  a      Square matrix to be checked
 **
 ** \param[out] b      Square cofactor
 **
 *****************************************************************************/
int matrix_cofactor(int neq, double *a, double *b)
{
  int i, j, ii, jj, i1, j1, neqm1;
  double *c, det;

  /* Core allocation */

  neqm1 = neq - 1;

  // Process the case when the matrix A is of dimension 1
  if (neqm1 <= 0)
  {
    B(0,0)= 1.;
    return 0;
  }
  c = (double*) mem_alloc(sizeof(double) * neqm1 * neqm1, 0);
  if (c == nullptr) return (1);

  /* Processing */

  for (j = 0; j < neq; j++)
  {
    for (i = 0; i < neq; i++)
    {

      /* Form the adjoint a_ij */

      i1 = 0;
      for (ii = 0; ii < neq; ii++)
      {
        if (ii == i) continue;
        j1 = 0;
        for (jj = 0; jj < neq; jj++)
        {
          if (jj == j) continue;
          C(i1,j1)= A(ii,jj);
          j1++;
        }
        i1++;
      }

      /* Calculate the determinate */
      det = matrix_determinant(neqm1, c);

      /* Fill in the elements of the cofactor */
      B(i,j)= pow(-1.0, i+j+2.0) * det;
    }
  }

  c = (double*) mem_free((char* ) c);
  return (0);
}

/****************************************************************************/
/*!
 **  Calculate the determinant of the triangular matrix after Cholesky
 **
 ** \return  Value of the determinant
 **
 ** \param[in]  neq Size of the matrix
 **
 ** \param[in]  tl  Lower triangular matrix defined by column
 **
 *****************************************************************************/
double matrix_cholesky_determinant(int neq, double *tl)
{
  int i;
  double deter;

  deter = 1.;
  for (i = 0; i < neq; i++)
    deter *= TL(i, i);
  return (deter);
}

/****************************************************************************/
/*!
 **  Check if a matrix is a rotation matrix
 **
 ** \return  1 if the matrix is a rotation matrix; 0 otherwise
 **
 ** \param[in]  neq      Size of the matrix
 ** \param[in]  a        Square matrix to be checked
 ** \param[in]  verbose  1 for the verbose option
 **
 ** \remark  A rotation matrix must be orthogonal with determinant equal to 1
 **
 *****************************************************************************/
int is_matrix_rotation(int neq, const double *a, int verbose)
{
  double deter, comp, prod;
  int i, j, k;

  /* Check product of matrix by its transpose and compare to unity matrix */

  for (i = 0; i < neq; i++)
    for (j = 0; j < neq; j++)
    {
      prod = 0.;
      for (k = 0; k < neq; k++)
        prod += A(i,k)* A(j,k);
      comp = (i == j) ? 1 :
                        0.;
      if (ABS(prod - comp) > Epsilon)
      {
        if (verbose)
          messerr("The element (A*At)[%d,%d] = %lf (should be %lf)", i + 1,
                  j + 1, prod, comp);
        return (0);
      }
    }

  deter = matrix_determinant(neq, a);

  if (ABS(deter - 1.) > Epsilon)
  {
    if (verbose) messerr("The Determinant = %f (should be 1)", deter);
    return (0);
  }

  return (1);
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
int matrix_cholesky_decompose(const double *a,
                                              double *tl,
                                              int neq)
{
  double prod;
  int ip, jp, kp;

  for (ip = 0; ip < neq; ip++)
    for (jp = 0; jp <= ip; jp++)
      TL(ip,jp)= AS(ip,jp);

  for (ip = 0; ip < neq; ip++)
  {
    prod = TL(ip, ip);
    for (kp = 0; kp < ip; kp++)
      prod -= TL(ip,kp)* TL(ip,kp);
    if (prod < 0.) return (ip + 1);
    TL(ip,ip)= sqrt(prod);

    for (jp = ip + 1; jp < neq; jp++)
    {
      prod = TL(jp, ip);
      for (kp = 0; kp < ip; kp++)
        prod -= TL(ip,kp)* TL(jp,kp);
      if (TL(ip,ip)<= 0.) return(ip+1);
      TL(jp,ip)= prod / TL(ip,ip);
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
                             double *tl,
                             double *a,
                             double *x)
{
  int irhs, i, j, n1, n2;
  double val, *v1, *v2;

  if (mode == 0)
  {
    for (irhs = 0; irhs < nrhs; irhs++)
      for (i = 0; i < neq; i++)
      {
        val = 0.;
        for (j = i; j < neq; j++)
          val += TL(j,i)* AS(j,irhs);
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
 **  Get the solution of a linear system (after Cholesky decomposition)
 **
 ** \return  Error returned code
 ** \return   0 : Success
 ** \return  -1 : Core allocation problem
 ** \return  -2 : Singular matrix (zero pivot)
 **
 ** \param[in]  neq  number of equations in the system
 ** \param[in]  tl   lower triangular matrix define by column
 ** \param[in]  b    matrix (dimension neq)
 **
 ** \param[out] x    resulting matrix (dimension neq)
 **
 *****************************************************************************/
int matrix_cholesky_solve(int neq,
                                          double *tl,
                                          double *b,
                                          double *x)
{
  int i, j, error;
  double sum, pivot, *r;

  /* Core allocation */

  error = -2;
  r = (double*) mem_alloc(sizeof(double) * neq, 0);
  if (r == nullptr) return (-1);

  /* Linear system associated to the lower triangular system */

  for (i = 0; i < neq; i++)
  {
    sum = b[i];
    for (j = 0; j < i; j++)
      sum -= r[j] * TL(i, j);
    pivot = TL(i, i);
    if (ABS(pivot) < _getTolInvert()) goto label_end;
    r[i] = sum / pivot;
  }

  /* Linear system associated to the upper triangular system */

  for (i = neq - 1; i >= 0; i--)
  {
    sum = r[i];
    for (j = i + 1; j < neq; j++)
      sum -= x[j] * TL(j, i);
    x[i] = sum / TL(i, i);
  }

  /* Set the error returned code */

  error = 0;

  label_end: r = (double*) mem_free((char* ) r);
  return (error);
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
void matrix_cholesky_invert(int neq, double *tl, double *xl)
{
  double sum;
  int i, j, l;

  for (i = 0; i < neq; i++)
  {
    for (j = 0; j < i; j++)
    {
      sum = 0.;
      for (l = j; l < i; l++)
        sum += TL(i,l)* XL(l,j);
      XL(i,j)= - sum / TL(i,i);
    }
    XL(i,i) = 1. / TL(i,i);
  }
}

/*****************************************************************************/
/*!
 **  Performs the product B = TL * A * TU or TU * A * TL
 **  where TL,TU is a triangular matrix and A a square symmetric matrix
 **
 ** \param[in]  mode  0: TL * A * TU; 1: TU * A * TL
 ** \param[in]  neq  number of equations in the system
 ** \param[in]  tl   Triangular matrix defined by column
 ** \param[in]  a    Square symmetric matrix (optional)
 **
 ** \param[out] b    Square symmetric matrix
 **
 *****************************************************************************/
void matrix_cholesky_norme(int mode,
                                           int neq,
                                           double *tl,
                                           double *a,
                                           double *b)
{
  int i, j, k, l;
  double val, vala;

  for (i = 0; i < neq; i++)
    for (j = 0; j < neq; j++)
    {
      val = 0.;
      if (mode == 0)
      {
        for (l = 0; l <= j; l++)
          for (k = 0; k <= i; k++)
          {
            if (a != nullptr)
              vala = AS(k, l);
            else
              vala = (k == l);
            val += TL(i,k)* vala * TL(j,l);
          }
        }
        else
        {
          for (l=j; l<neq; l++)
          for (k=i; k<neq; k++)
          {
            if (a != nullptr)
            vala = AS(k,l);
            else
            vala = (k == l);
            val += TL(k,i) * vala * TL(l,j);
          }
        }
      BS(i,j)= val;
    }
  return;
}

/*****************************************************************************/
/*!
 **  Invert the matrix (after Cholesky decomposition)
 **
 ** \return  Error returned code
 ** \return   0 : Success
 ** \return  -1 : Core allocation problem
 ** \return  -2 : Singular matrix (zero pivot)
 **
 ** \param[in]  neq  number of equations in the system
 ** \param[in]  tl   lower triangular matrix defined by column
 **
 ** \param[out] xl   invert lower triangular matrix
 **
 *****************************************************************************/
int matrix_cholesky_to_invert(int neq, double *tl, double *xl)
{
  int i, j, k, error;
  double *r, *ek, *res, sum, pivot;

  /* Initializations */

  r = ek = res = nullptr;

  /* Look for a zero pivot */

  for (i = 0; i < neq; i++)
  {
    pivot = TL(i, i);
    if (ABS(pivot) < _getTolInvert()) return (-2);
  }

  /* Core allocation */

  error = -1;
  r = (double*) mem_alloc(sizeof(double) * neq, 0);
  if (r == nullptr) goto label_end;
  ek = (double*) mem_alloc(sizeof(double) * neq, 0);
  if (ek == nullptr) goto label_end;
  res = (double*) mem_alloc(sizeof(double) * neq, 0);
  if (res == nullptr) goto label_end;

  /* Loop on the canonical basis */

  for (k = 0; k < neq; k++)
  {
    for (i = 0; i < neq; i++)
      ek[i] = 0.;
    ek[k] = 1.;

    /* Linear system associated to the lower triangular factor */

    for (i = 0; i < neq; i++)
    {
      sum = ek[i];
      for (j = 0; j < i; j++)
        sum -= r[j] * TL(i, j);
      pivot = TL(i, i);
      r[i] = sum / pivot;
    }

    /* Linear system associated to the upper triangular factor */

    for (i = neq - 1; i >= k; i--)
    {
      sum = r[i];
      for (j = i + 1; j < neq; j++)
        sum -= res[j] * TL(j, i);
      res[i] = sum / TL(i, i);
      XL(i,k)= res[i];
    }
  }

  /* Set the error returned code */

  error = 0;

  label_end: r = (double*) mem_free((char* ) r);
  ek = (double*) mem_free((char* ) ek);
  res = (double*) mem_free((char* ) res);
  return (error);
}

/*****************************************************************************/
/*!
 **  Performs the product of a symmetric matrix by a vector
 **
 ** \param[in]  neq    Dimension of the matrix
 ** \param[in]  mode   1 if the Lower matrix is stored linewise
 **                      (or if the Upper matrix is stored columnwise)
 **                    2 if the Lower matrix is stored columnwise
 **                      (or the Upper matrix is stored linewise)
 ** \param[in]  al     Lower triangular matrix defined by column
 ** \param[in]  b      Vector
 **
 ** \param[out] x      Resulting product vector
 **
 *****************************************************************************/
void matrix_triangular_product(int neq,
                                               int mode,
                                               const double *al,
                                               const double *b,
                                               double *x)
{
  int i, j;
  const double *at;
  double value;

  if (mode == 1)
  {
    at = al;
    for (i = 0; i < neq; i++)
    {
      value = 0.;
      for (j = 0; j <= i; j++)
        value += AT(j,i)* b[j];
      for (j = i + 1; j < neq; j++)
        value += AT(i,j)* b[j];
      x[i] = value;
    }
  }
  else
  {
    for (i = 0; i < neq; i++)
    {
      value = 0.;
      for (j = 0; j <= i; j++)
        value += AL(i,j)* b[j];
      for (j = i + 1; j < neq; j++)
        value += AL(j,i)* b[j];
      x[i] = value;
    }
  }
  return;
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
 ** \remark The input and output matrices can match
 **
 *****************************************************************************/
int matrix_invgen(double *a,
                                  int neq,
                                  double *tabout,
                                  double *cond)
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
      TABOUT(i,j)= value;
    }

    /* Set the error returned code */

  error = 0;

  label_end: eigval = (double*) mem_free((char* ) eigval);
  eigvec = (double*) mem_free((char* ) eigvec);
  return (error);
}

/****************************************************************************/
/*!
 **  Calculate the norm of a vector
 **
 ** \return  Value of the norm
 **
 ** \param[in]  a         Input vector
 ** \param[in]  neq       Number of equations
 **
 *****************************************************************************/
double matrix_norm(double *a, int neq)
{
  double value;
  int i;

  value = 0.;
  for (i = 0; i < neq; i++)
    value += a[i] * a[i];
  return (value);
}

/****************************************************************************/
/*!
 **  Calculate the vector product of two vectors (only in 3-D)
 **
 ** \param[in]  a         First vector
 ** \param[in]  b         Second vector
 **
 ** \param[out] v         Vector product
 **
 *****************************************************************************/
void vector_product(double *a, double *b, double *v)
{
  v[0] = a[1] * b[2] - a[2] * b[1];
  v[1] = a[2] * b[0] - a[0] * b[2];
  v[2] = a[0] * b[1] - a[1] * b[0];
}

/****************************************************************************/
/*!
 ** Translate a point by a vector
 **
 ** \param[in]  ndim      Space dimension
 ** \param[in]  a         Coordinates of the target point (Dimension: ndim)
 ** \param[in]  v         Vector coordinates (Dimension: ndim)
 **
 ** \param[out] b         Coordinates of the resulting point (Dimension: ndim)
 **
 *****************************************************************************/
void vector_translate(int ndim, double *a, double *v, double *b)
{
  for (int i = 0; i < ndim; i++)
    b[i] = a[i] + v[i];
}

/****************************************************************************/
/*!
 **  Calculate the inner product of two vectors
 **
 ** \return  Value of the inner product
 **
 ** \param[in]  a         First vector
 ** \param[in]  b         Second vector
 ** \param[in]  neq       Space dimension
 **
 *****************************************************************************/
double inner_product(const double *a, const double *b, int neq)
{
  double value;
  int i;

  value = 0.;
  for (i = 0; i < neq; i++)
    value += a[i] * b[i];
  return (value);
}

/****************************************************************************/
/*!
 **  Calculate the square norm of a vector issued from the
 **  inner product relative to a matrix A
 **
 ** \return  Value of the norm
 **
 ** \param[in]  b         Vector (Dimension: neq)
 ** \param[in]  a         Matrix (Symmetric, Dimension neq x neq)
 ** \param[in]  neq       Space dimension
 ** \param[in]  subneq    Dimension of the subspace to be considered
 **
 *****************************************************************************/
double matrix_normA(double *b, double *a, int neq, int subneq)
{
  double value;
  int i, j;

  value = 0.;
  for (i = 0; i < (subneq - 1); i++)
    for (j = i + 1; j < subneq; j++)
      value += b[i] * A(i, j) * b[j];
  value *= 2.;

  for (i = 0; i < neq; i++)
    value += b[i] * A(i, i) * b[i];

  return (value);
}

/****************************************************************************/
/*!
 **  Perform the partial-pivoting Gaussian elimination
 **
 ** \param[in]  a         Matrix to be inverted
 ** \param[in]  neq       Matrix dimension
 ** \param[in]  indx      Working array (Dimension: neq)
 ** \param[in]  c         Working array (Dimension: neq)
 **
 *****************************************************************************/
static void st_elgs(double *a, int neq, int *indx, double *c)
{
  int i, j, k, itmp;
  double c1, pi, pi1, pj;

  /* Initialize the index */

  k = 0;
  for (i = 0; i < neq; i++)
    indx[i] = i;

  /* Find the rescaling factors, one from each row */

  for (i = 0; i < neq; i++)
  {
    c1 = 0;
    for (j = 0; j < neq; j++)
    {
      if (ABS(A(i,j)) > c1) c1 = ABS(A(i,j));
    }
    c[i] = c1;
  }

  /* Search the pivoting (largest) element from each column */

  for (j = 0; j < neq - 1; j++)
  {
    pi1 = 0;
    for (i = j; i < neq; i++)
    {
      pi = ABS(A(indx[i],j)) / c[indx[i]];
      if (pi > pi1)
      {
        pi1 = pi;
        k = i;
      }
    }

    /* Interchange the rows via indx[] to record pivoting order */

    itmp = indx[j];
    indx[j] = indx[k];
    indx[k] = itmp;
    for (i = j + 1; i < neq; i++)
    {
      pj = A(indx[i],j)/ A(indx[j],j);

      /* Record pivoting ratios below the diagonal */

      A(indx[i],j) = pj;

      /* Modify other elements accordingly */

      for (k = j+1; k < neq; k++)
      A(indx[i],k) -= pj * A(indx[j],k);
    }
  }
}

/*****************************************************************************/
/*!
 **  Invert a symmetric square matrix
 **
 ** \return  Return code: 0 no error; k if the k-th pivot is zero
 **
 ** \param[in,out] a input matrix, destroyed in computation and replaced by
 **                  resultant inverse
 ** \param[in]  neq  number of equations in the matrix 'a'
 **
 *****************************************************************************/
int matrix_invsym(double *a, int neq)
{
  double *x, *b, *c, ratio;
  int *indx, i, j, k;

  /* Initializations */

  b = c = x = nullptr;
  indx = nullptr;

  /* Core allocation */

  b = (double*) mem_alloc(sizeof(double) * neq * neq, 1);
  x = (double*) mem_alloc(sizeof(double) * neq * neq, 1);
  c = (double*) mem_alloc(sizeof(double) * neq, 1);
  indx = (int*) mem_alloc(sizeof(int) * neq, 1);

  /* Processing */

  for (i = 0; i < neq; i++)
    for (j = 0; j < neq; j++)
      B(i,j)= 0.;

  for (i = 0; i < neq; i++)
    B(i,i)= 1.;

  st_elgs(a, neq, indx, c);

  for (i = 0; i < neq - 1; i++)
    for (j = i + 1; j < neq; j++)
      for (k = 0; k < neq; k++)
        B(indx[j],k)= B(indx[j],k) - A(indx[j],i) * B(indx[i],k);

  for (i = 0; i < neq; i++)
  {
    ratio = A(indx[neq - 1], neq - 1);
    if (ABS(ratio) < _getTolInvert()) return (neq);
    X(neq-1,i)= B(indx[neq-1],i) / ratio;
    for (j = neq - 2; j >= 0; j = j - 1)
    {
      X(j,i)= B(indx[j],i);
      for (k = j+1; k < neq; k++)
      X(j,i) = X(j,i) - A(indx[j],k) * X(k,i);
      ratio = A(indx[j],j);
      if (ABS(ratio) < _getTolInvert()) return(indx[j]+1);
      X(j,i) /= ratio;
    }
  }

  /* Copy the result */

  for (i = 0; i < neq; i++)
    for (j = 0; j < neq; j++)
      A(i,j)= X(i,j);

      /* Core deallocation */

  b = (double*) mem_free((char* ) b);
  c = (double*) mem_free((char* ) c);
  c = (double*) mem_free((char* ) x);
  indx = (int*) mem_free((char* ) indx);
  return (0);
}

/*****************************************************************************/
/*!
 **  Invert the sign of a matrix
 **
 ** \param[in]  neq    number of cells in the matrix
 ** \param[in,out] a   input/output matrix
 **
 *****************************************************************************/
void matrix_invsign(int neq, double *a)
{
  int i;

  for (i = 0; i < neq; i++)
    a[i] = -a[i];
}

/*****************************************************************************/
/*!
 **  Fill a square matrix with a triangular matrix
 **
 ** \param[in]  mode   0: TL (upper); 1: TL (lower)
 ** \param[in]  neq    number of equations in the system
 ** \param[in]  tl     Triangular matrix (lower part)
 **
 ** \param[out] a      Resulting square matrix
 **
 *****************************************************************************/
void matrix_triangle_to_square(int mode,
                                               int neq,
                                               double *tl,
                                               double *a)
{
  for (int i = 0; i < neq * neq; i++)
    a[i] = 0.;

  for (int i = 0; i < neq; i++)
    for (int j = 0; j < neq; j++)
    {
      if (mode == 0)
      {
        if (j <= i) AS(i,j)= TL(i,j);
      }
      else
      {
        if (j >= i) AS(i,j) = TL(j,i);
      }
    }
  }

  /*****************************************************************************/
  /*!
   **  Transform a symmetrical matrix (entered as triangle)
   **  into a square matrix
   **
   ** \param[in]  neq    number of equations in the system
   ** \param[in]  tl     Upper Triangular matrix (columnwise)
   **
   ** \param[out] a      Resulting square matrix
   **
   *****************************************************************************/
void matrix_tri2sq(int neq, double *tl, double *a)
{
  int i, j;

  for (i = 0; i < neq; i++)
    for (j = 0; j < neq; j++)
    {
      AS(i,j)= (j < i) ? TL(i,j) : TL(j,i);
    }
  }

  /*****************************************************************************/
  /*!
   **  Transform a square matrix into a triangular one
   **
   ** \param[in]  mode   0: TL (upper); 1: TL (lower)
   ** \param[in]  neq    number of equations in the system
   ** \param[in]  a      Input square (symmetric) matrix
   **
   ** \param[out] tl     Triangular matrix (lower part)
   **
   ** \remark: No test is performed to check that the input matrix is symmetric
   **
   *****************************************************************************/
void matrix_square_to_triangle(int mode,
                                               int neq,
                                               double *a,
                                               double *tl)
{
  int i, j;

  for (i = 0; i < neq; i++)
    for (j = 0; j < neq; j++)
    {
      if (mode == 0)
      {
        if (j <= i) TL(i,j)= AS(i,j);
      }
      else
      {
        if (j >= i) TL(j,i) = AS(i,j);
      }
    }
  }

  /*****************************************************************************/
  /*!
   **  Calculate the product of 'tl' (lower triangle) by its transpose
   **
   ** \param[in]  neq    number of equations in the system
   ** \param[in]  tl     Lower triangular matrix defined by column
   **
   ** \param[out] a      Resulting square matrix
   **
   *****************************************************************************/
void matrix_produit_lu(int neq, double *tl, double *a)
{
  int i, j, k;

  for (i = 0; i < neq; i++)
    for (j = 0; j < neq; j++)
    {
      A(i,j)= 0.;
      for (k=0; k<neq; k++)
      {
        if (k > i || k > j) continue;
        A(i,j) += TL(i,k) * TL(j,k);
      }
    }
  return;
}

/*****************************************************************************/
/*!
 **  Calculate the product of 'tl' (lower triangle) by its transpose
 **  and store the result in a VectorDouble
 **
 **  \return The resulting VectorDouble
 **
 ** \param[in]  neq    number of equations in the system
 ** \param[in]  tl     Lower triangular matrix defined by column
 **
 *****************************************************************************/
VectorDouble matrix_produit_lu_VD(int neq, double *tl)
{
  VectorDouble a;
  a.resize(neq * neq);

  matrix_produit_lu(neq, tl, a.data());
  return a;
}

/*****************************************************************************/
/*!
 **  Checks if the product of the matrix by its inverse
 **  is close enough to the identity
 **
 ** \return  1 if the product is identity; 0 otherwise
 **
 ** \param[in]  a      Matrix
 ** \param[in]  b      Inverse Matrix
 ** \param[in]  neq    number of equations in the system
 **
 ** \param[out] errmax  Maximum error encountered
 **
 *****************************************************************************/
int is_matrix_product_identity(int neq,
                                               double *a,
                                               double *b,
                                               double *errmax)
{
  double *x, compare, valmax;
  int i, j, error;

  x = (double*) mem_alloc(sizeof(double) * neq * neq, 1);

  matrix_product(neq, neq, neq, a, b, x);

  valmax = 0.;
  for (i = 0; i < neq; i++)
    for (j = 0; j < neq; j++)
    {
      compare = (i == j) ? 1. :
                           0.;
      valmax = ABS(X(i,j) - compare);
    }
  error = (valmax <= Epsilon);

  *errmax = valmax;
  x = (double*) mem_free((char* ) x);
  return (error);
}

/*****************************************************************************/
/*!
 **  Invert a symmetric square matrix (SVD method)
 **
 ** \return Error return code
 **
 ** \param[in]  a2   pre-allocated with 2n rows and n columns. The matrix to be
 **                  decomposed is contained in the first n rows of A.
 **                  On return, the n first rows of A contain the product US
 **                  and the lower n rows contain V (not V').
 ** \param[in]  neq  number of equations in the matrix 'a'
 **
 ** \param[out] s    the square of the singular values.
 **
 *****************************************************************************/
static int st_svd(double *a2, double *s, int neq)
{
  int i, j, k, EstColRank, RotCount, SweepCount, slimit;
  double e2, tol, vt, p, x0, y0, q, r, c0, s0, d1, d2;

  /* Initializations */

  SweepCount = 0;
  EstColRank = neq;
  RotCount = neq;
  slimit = (neq < 120) ? 30 :
                         neq / 4;
  e2 = 10.0 * neq * _getEpsSVD() * _getEpsSVD();
  tol = 0.10 * _getEpsSVD();

  /* Processing */

  for (i = 0; i < neq; i++)
  {
    for (j = 0; j < neq; j++)
      A2(neq+i,j)= 0.0;
      A2(neq+i,i) = 1.0;
    }

  while (RotCount != 0 && SweepCount++ <= slimit)
  {
    RotCount = EstColRank * (EstColRank - 1) / 2;
    for (j = 0; j < EstColRank - 1; j++)
      for (k = j + 1; k < EstColRank; k++)
      {
        p = q = r = 0.0;
        for (i = 0; i < neq; i++)
        {
          x0 = A2(i, j);
          y0 = A2(i, k);
          p += x0 * y0;
          q += x0 * x0;
          r += y0 * y0;
        }
        s[j] = q;
        s[k] = r;
        if (q >= r)
        {
          if (q <= e2 * s[0] || ABS(p) <= tol * q)
            RotCount--;
          else
          {
            p /= q;
            r = 1.0 - r / q;
            vt = sqrt(4. * p * p + r * r);
            c0 = sqrt(0.5 * (1.0 + r / vt));
            s0 = p / (vt * c0);
            for (i = 0; i < 2 * neq; i++)
            {
              d1 = A2(i, j);
              d2 = A2(i, k);
              A2(i,j)= d1*c0+d2*s0;
              A2(i,k)= -d1*s0+d2*c0;
            }
          }
        }
        else
        {
          p /= r;
          q = q/r - 1.0;
          vt = sqrt(4.0*p*p + q*q);
          s0 = sqrt(0.5*(1.0 - q/vt));
          if (p < 0.) s0 = -s0;
          c0 = p / (vt*s0);
          for (i=0; i<2*neq; i++)
          {
            d1 = A2(i,j);
            d2 = A2(i,k);
            A2(i,j) = d1*c0+d2*s0;
            A2(i,k) = -d1*s0+d2*c0;
          }
        }
      }
    while (EstColRank > 2 && s[EstColRank - 1] <= s[0] * tol + tol * tol)
      EstColRank--;
  }

  if (SweepCount > slimit)
  {
    messerr("Warning: Reached maximum number of sweeps (%d) in SVD routine",
            slimit);
    return (1);
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Invert a symmetric square matrix (SVD method)
 **
 ** \param[in]  neq  number of equations in the matrix 'a'
 ** \param[in]  s    input matrix (square)
 **
 ** \param[out] u       first decomposition matrix
 ** \param[out] v       second decomposition matrix
 ** \param[out] tabout  inverse matrix
 **
 ** \remark This routine has been adapted from a Pascal implementation
 ** \remark (c) 1988 J. C. Nash, "Compact numerical methods for computers",
 ** \remark Hilger 1990.
 **
 *****************************************************************************/
void matrix_svd_inverse(int neq,
                                        double *s,
                                        double *u,
                                        double *v,
                                        double *tabout)
{
  int i, j, k, lec, number;
  double maxval, thresh, value;

  /* Unloading the results */

  number = 0;
  maxval = 0.;
  for (i = 0; i < neq; i++)
    if (s[i] > maxval) maxval = s[i];
  thresh = _getEpsSVD() * neq * maxval;
  thresh = 1.e-06;
  for (i = 0; i < neq; i++)
    if (s[i] < thresh) number++;

  /* Compute the inverse */

  for (i = lec = 0; i < neq; i++)
    for (j = 0; j < neq; j++, lec++)
    {
      value = 0.;
      for (k = 0; k < neq; k++)
        if (s[k] > thresh) value += V(k,i)* (1./s[k]) * U(k,j);
      tabout[lec] = value;
    }
}

/*****************************************************************************/
/*!
 **  Invert a symmetric square matrix (SVD method)
 **
 ** \return  Error return code
 **
 ** \param[in]  mat  input matrix, destroyed in computation and replaced by
 **                  resultant inverse
 ** \param[in]  neq  number of equations in the matrix 'a'
 **
 ** \param[out] rank instability rank (if > 0)
 **
 ** \remark This routine has been adapted from a Pascal implementation
 ** \remark (c) 1988 J. C. Nash, "Compact numerical methods for computers",
 ** \remark Hilger 1990.
 **
 *****************************************************************************/
int matrix_invsvdsym(double *mat, int neq, int rank)
{
  int i, j, lec, error;
  double *a2, *s, *u, *v;

  /* Core allocation */

  error = 1;
  a2 = (double*) mem_alloc(2 * neq * neq * sizeof(double), 1);
  u = (double*) mem_alloc(neq * neq * sizeof(double), 1);
  v = (double*) mem_alloc(neq * neq * sizeof(double), 1);
  s = (double*) mem_alloc(neq * sizeof(double), 1);

  /* Loading */

  for (i = lec = 0; i < neq; i++)
    for (j = 0; j < neq; j++, lec++)
      A2(i,j)= mat[lec];

      /* Calling the SVD routine */

  if (st_svd(a2, s, (int) neq)) goto label_end;

  /* Unloading the results */

  for (i = 0; i < neq; i++)
    s[i] = sqrt(s[i]);

  for (i = 0; i < neq; i++)
    for (j = 0; j < neq; j++)
      U(j,i)= (s[j] > 0.) ? A2(i,j) / s[j] : 0.;

  for (i = 0; i < neq; i++)
    for (j = 0; j < neq; j++)
      V(j,i)= A2(i+neq,j);

      /* Compute the inverse */

  matrix_svd_inverse(neq, s, u, v, mat);

  /* Set the error return code */

  error = 0;

  /* Core deallocation */

  label_end: if (error)
  {
    if (rank >= 0)
      messerr("Error in SVD decomposition (rank=%d)", rank + 1);
    else
      messerr("Error in SVD decomposition");
  }
  a2 = (double*) mem_free((char* ) a2);
  u = (double*) mem_free((char* ) u);
  v = (double*) mem_free((char* ) v);
  s = (double*) mem_free((char* ) s);

  return (error);
}

/*****************************************************************************/
/*!
 **  Check if the index match one element of a list
 **
 ** \return  1 if the target index is present in the list; -1 otherwise
 **
 ** \param[in]  index   Target index
 ** \param[in]  nitem   Number of items in the list
 ** \param[in]  items   List of items
 **
 *****************************************************************************/
static int st_match_index(int index, int nitem, int *items)
{
  int i;

  if (nitem <= 0) return (-1);
  for (i = 0; i < nitem; i++)
  {
    if (items[i] == index) return (1);
  }
  return (-1);
}

/*****************************************************************************/
/*!
 **  Manage a matrix and derive a submatrix
 **
 ** \param[in]  nrows   Number of rows of the input matrix
 ** \param[in]  ncols   Number of columns of the input matrix
 ** \param[in]  nr      Number of rows of interest
 ** \li                  >0 : for extraction
 ** \li                   0 : no action
 ** \li                  <0 : for suppression
 ** \param[in]  nc      Number of columns of interest
 ** \li                  >0 : for extraction
 ** \li                   0 : no action
 ** \li                  <0 : for suppression
 ** \param[in]  rowsel  Array of rows indices of interest (starting from 0)
 **                     (Dimension: ABS(nr)
 ** \param[in]  colsel  Array of columns indices of interest (starting from 0)
 **                     (Dimension: ABS(nr)
 ** \param[in]  v1      Input rectangular matrix (Dimension: nrows * ncols)
 **
 ** \param[out] v2      Output rectangular matrix
 **
 *****************************************************************************/
void matrix_manage(int nrows,
                   int ncols,
                   int nr,
                   int nc,
                   int *rowsel,
                   int *colsel,
                   double *v1,
                   double *v2)
{
  int irow, icol, ecr, lec, flag_col, flag_row;

  lec = ecr = 0;
  for (icol = 0; icol < ncols; icol++)
  {
    flag_col = st_match_index(icol, ABS(nc), colsel);
    for (irow = 0; irow < nrows; irow++, lec++)
    {
      flag_row = st_match_index(irow, ABS(nr), rowsel);
      if (nr * flag_row < 0 || nc * flag_col < 0) continue;
      v2[ecr++] = v1[lec];
    }
  }
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
 **  Fill the matrix to take the symmetry into account
 **
 ** \param[in]  neq   matrix dimension
 ** \param[in,out]  a square symmetric matrix (dimension = neq*neq)
 **
 ** \remark  We assume that, in the input matrix, the elements A[i,j] where
 ** \remark  i >= j have already been filled
 **
 *****************************************************************************/
void matrix_fill_symmetry(int neq, double *a)
{
  int i, j;

  for (i = 0; i < neq; i++)
    for (j = i; j < neq; j++)
      A(j,i)= A(i,j);
    }

    /*****************************************************************************/
    /*!
     **  Calculate the maximum of the absolute values of a vector
     **
     ** \param[in]  nval  vector dimension
     ** \param[in]  tab   vector
     **
     *****************************************************************************/
double matrix_norminf(int nval, double *tab)
{
  double value, retval;
  int i;

  retval = 0.;
  for (i = 0; i < nval; i++)
  {
    value = ABS(tab[i]);
    if (value > retval) retval = value;
  }
  return (retval);
}

/*****************************************************************************/
/*!
 **  Calculate the sum of the absolute values of a vector
 **
 ** \param[in]  nval  vector dimension
 ** \param[in]  tab   vector
 **
 *****************************************************************************/
double matrix_norml1(int nval, double *tab)
{
  double retval;
  int i;

  retval = 0.;
  for (i = 0; i < nval; i++)
    retval += ABS(tab[i]);
  return (retval);
}

/*****************************************************************************/
/*!
 **  Performs the t(A) %*% A where A is a square matrix
 **
 ** \param[in]  neq matrix dimension for A
 ** \param[in]  a   square matrix
 **
 ** \param[out] b   square matrix
 **
 ** \remark Matrices a() and b() may not coincide
 **
 *****************************************************************************/
void matrix_square(int neq, double *a, double *b)
{
  int i1, i2, i3;
  double value;

  for (i1 = 0; i1 < neq; i1++)
    for (i3 = 0; i3 < neq; i3++)
    {
      value = 0.;
      for (i2 = 0; i2 < neq; i2++)
        value += A(i1,i2)* A(i3,i2);
      B(i1,i3)= value;
    }

  return;
}

/*****************************************************************************/
/*!
 **  Performs the B = t(A) %*% A where A is a square matrix
 **
 ** \return the returned VectorDouble
 **
 ** \param[in]  neq matrix dimension for A
 ** \param[in]  a   square matrix (VectorDouble)
 **
 *****************************************************************************/
VectorDouble matrix_square_VD(int neq, const VectorDouble &a)
{
  VectorDouble b;
  int i1, i2, i3;
  double value;

  b.resize(neq * neq);
  for (i1 = 0; i1 < neq; i1++)
    for (i3 = 0; i3 < neq; i3++)
    {
      value = 0.;
      for (i2 = 0; i2 < neq; i2++)
        value += A(i1,i2)* A(i3,i2);
      B(i1,i3)= value;
    }

  return b;
}

/*****************************************************************************/
/*!
 **  Performs the A %*% diag(c) where A is a square matrix and c a vector
 **
 ** \param[in]  mode  0: c as is; 1: sqrt(c); 2: 1/c; 3: 1/sqrt(c)
 ** \param[in]  neq   matrix dimension for A
 ** \param[in]  a     square matrix
 ** \param[in]  c     vector
 **
 ** \param[out] b   square matrix
 **
 ** \remark Matrices a() and b() may coincide
 **
 *****************************************************************************/
void matrix_product_by_diag(int mode,
                                            int neq,
                                            double *a,
                                            double *c,
                                            double *b)
{
  double val;
  int i1, i2;

  for (i1 = 0; i1 < neq; i1++)
    for (i2 = 0; i2 < neq; i2++)
    {
      val = c[i2];
      if (mode == 1)
        val = sqrt(val);
      else if (mode == 2)
        val = 1. / val;
      else if (mode == 3) val = 1. / sqrt(val);
      B(i1,i2)= A(i1,i2) * val;
    }

  return;
}

/*****************************************************************************/
/*!
 **  Performs the A <- A %*% diag(c) where A is a square matrix and c a vector
 **
 ** \param[in]  mode  0: c as is; 1: sqrt(c); 2: 1/c; 3: 1/sqrt(c)
 ** \param[in]  neq   matrix dimension for A
 ** \param[in]  a     square matrix
 ** \param[in]  c     vector
 **
 *****************************************************************************/
void matrix_product_by_diag_VD(int mode,
                                               int neq,
                                               VectorDouble a,
                                               const VectorDouble &c)
{
  double val;
  int i1, i2;

  for (i1 = 0; i1 < neq; i1++)
    for (i2 = 0; i2 < neq; i2++)
    {
      val = c[i2];
      if (mode == 1)
        val = sqrt(val);
      else if (mode == 2)
        val = 1. / val;
      else if (mode == 3) val = 1. / sqrt(val);
      A(i1,i2)*= val;
    }

  return;
}

/*****************************************************************************/
/*!
 **  Performs the a1 * A + b1 * B
 **
 ** \param[in]  neq matrix dimension for A
 ** \param[in]  a1  Coefficient applied to matrix a()
 ** \param[in]  a   square matrix
 ** \param[in]  b1  Coefficient applied to matrix b()
 ** \param[in]  b   square matrix
 **
 ** \param[out] x   square matrix
 **
 ** \remark Matrix x() can coincide with a() or b()
 **
 *****************************************************************************/
void matrix_linear(int neq,
                                   double a1,
                                   double *a,
                                   double b1,
                                   double *b,
                                   double *x)
{
  int i1, i2;

  for (i1 = 0; i1 < neq; i1++)
    for (i2 = 0; i2 < neq; i2++)
      X(i1,i2)= a1 * A(i1,i2) + b1 * B(i1,i2);
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
    if (h == 0)
    {
      if (vecinf[i] != 0. || vecsup[i - 1] != 0.) return (2);
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
    if (e[i] != 0.0)
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

/*****************************************************************************/
/*!
 **  Solve a linear system: H %*% g = x
 **
 ** \return  Error return code
 **
 ** \param[in]  neq     number of equations in the system
 ** \param[in]  hmat    symmetric square matrix (Dimension: neq * neq)
 ** \param[in]  gmat    right-hand side vector (Dimension: neq)
 **
 ** \param[out] xmat    solution vector (Dimension: neq)
 **
 ** \remark In output, hmat contains the inverse matrix
 **
 *****************************************************************************/
int matrix_qo(int neq, double *hmat, double *gmat, double *xmat)
{
  int error;
  double cond;
  error = matrix_invgen(hmat, neq, hmat, &cond);
  if (error || cond > 1e13) return (1);
  matrix_product(neq, neq, 1, hmat, gmat, xmat);
  return (0);
}

/*****************************************************************************/
/*!
 **  Minimize 1/2 t(x) %*% H %*% x + t(g) %*% x under the constraints
 **  t(A) %*% x = b
 **
 ** \return  Error return code
 **
 ** \param[in]  flag_invert Tells if the inverse has already been calculated
 ** \param[in]  neq    number of equations in the system
 ** \param[in]  hmat   symmetric square matrix (Dimension: neq * neq)
 ** \param[in]  gmat   right-hand side vector (Dimension: neq)
 ** \param[in]  na     Number of equalities
 ** \param[in]  amat   matrix for inequalities (Dimension: neq * na)
 ** \param[in]  bmat   inequality vector (Dimension: na)
 ** \param[in]  xmat   solution of the linear system with no constraint.
 **                    On return, solution with constraints (Dimension: neq)
 **
 ** \param[out] lambda working vector (Dimension: na)
 **
 ** \remark In input:
 ** \remark If flag_invert== 1, H is provided as the generalized inverse
 ** \remark and x contains the solution of the linear system with no constraint
 ** \remark If flag_invert==0, H is the primal matrix
 **
 ** \remark In output, H contains the inverse matrix
 **
 *****************************************************************************/
int matrix_qoc(int flag_invert,
                               int neq,
                               double *hmat,
                               double *gmat,
                               int na,
                               double *amat,
                               double *bmat,
                               double *xmat,
                               double *lambda)
{
  int i, j, k, error, error_int;
  double *ha, *temp, *evec, value, cond;

  /* Initializations */

  error = 1;
  temp = evec = ha = nullptr;

  /* Core allocation */

  ha = (double*) mem_alloc(sizeof(double) * neq * na, 1);
  temp = (double*) mem_alloc(sizeof(double) * na * na, 1);
  evec = (double*) mem_alloc(sizeof(double) * na, 1);

  /* Preliminary solution of the linear system with no constraint */

  if (!flag_invert)
  {
    if (matrix_qo(neq, hmat, gmat, xmat)) goto label_end;
  }

  /* Product HA = H %*% A */

  for (i = 0; i < neq; i++)
    for (j = 0; j < na; j++)
    {
      value = 0.;
      for (k = 0; k < neq; k++)
        value += HMT(i,k)* AMT(k,j);
      HA(i,j)= value;
    }

    /* Product TEMP = t(A) %*% H %*% A */

  for (i = 0; i < na; i++)
    for (j = 0; j < na; j++)
    {
      value = 0.;
      for (k = 0; k < neq; k++)
        value += AMT(k,i)* HA(k,j);
      TEMP(i,j)= value;
    }

    /* Generalized inverse of TEMP */

  error_int = matrix_invgen(temp, na, temp, &cond);
  if (error_int || cond > 1e13) goto label_end;

  /* Evaluate evec = t(A) %*% x - b */

  for (i = 0; i < na; i++)
  {
    value = 0.;
    for (j = 0; j < neq; j++)
      value += AMT(j,i)* xmat[j];
    evec[i] = value - bmat[i];
  }

  /* Evaluate lambda = TEMP %*% evec */

  for (i = 0; i < na; i++)
  {
    value = 0.;
    for (j = 0; j < na; j++)
      value += TEMP(i,j)* evec[j];
    lambda[i] = value;
  }

  /* Evaluate x = x - H %*% A %*% lambda */

  for (i = 0; i < neq; i++)
  {
    value = 0.;
    for (j = 0; j < na; j++)
      value += HA(i,j)* lambda[j];
    xmat[i] -= value;
  }

  /* Set the error return code */

  error = 0;

  label_end: ha = (double*) mem_free((char* ) ha);
  temp = (double*) mem_free((char* ) temp);
  evec = (double*) mem_free((char* ) evec);
  return (error);
}

/*****************************************************************************/
/*!
 **  Count the number of active constraints
 **
 ** \return  Number of active constraints
 **
 ** \param[in]  nai    Number of constraints
 ** \param[in]  active Array of constraint status
 **
 *****************************************************************************/
static int st_count_active(int nai, int *active)
{
  int number, i;

  number = 0;
  for (i = 0; i < nai; i++)
    if (active[i]) number++;
  return (number);
}

/*****************************************************************************/
/*!
 **  Concatenate the equality and the active inequality material
 **
 **  \return The total number of constraints
 **
 ** \param[in]  nae      Number of equalities
 ** \param[in]  nai      Number of inequalities
 ** \param[in]  active   Array of active/non active inequalities
 ** \param[in]  neq      First dimension of the array
 ** \param[in]  tabemat  Equality material (Dimension: neq * nai)
 ** \param[in]  tabimat  Inequality material
 **
 ** \param[out] tabout   Output array
 **
 *****************************************************************************/
static int st_copy_active(int nae,
                          int nai,
                          int *active,
                          int neq,
                          double *tabemat,
                          double *tabimat,
                          double *tabout)
{
  int i, j, number;

  /* Copy the equalities */

  number = 0;
  for (i = 0; i < nae; i++)
  {
    for (j = 0; j < neq; j++)
      TABOUT(j,number)= TABEMT(j,i);
      number++;
    }

    /* Copy the active inequalities */

  for (i = 0; i < nai; i++)
  {
    if (active[i] == 0) continue;
    for (j = 0; j < neq; j++)
      TABOUT(j,number)= TABIMT(j,i);
    number++;
  }
  return (number);
}

/*****************************************************************************/
/*!
 **  Calculate how constraints are fulfilled
 **
 **  \return Count of the constraints not fulfilled
 **
 ** \param[in]  neq      First dimension of the array
 ** \param[in]  nai      Number of inequalities
 ** \param[in]  active   Array of active/non active inequalities (optional)
 ** \param[in]  aimat    Inequality material (Dimension: neq * nai)
 ** \param[in]  bimat    right-hand side for inequalities (Dimension: nai)
 ** \param[out] xmat     solution of the linear system with no constraint (neq)
 **
 ** \param[out] vmat     matrix of errors (if not NULL)
 ** \param[out] flag     array specifying if constraint is active (if not NULL)
 **
 *****************************************************************************/
static int st_calcul_error(int neq,
                           int nai,
                           int *active,
                           double *aimat,
                           double *bimat,
                           double *xmat,
                           double *vmat,
                           int *flag)
{
  int i, j, ecr, number, flag_active;
  double value, ecart;
  double eps = 1.e-10;

  number = ecr = 0;
  for (i = 0; i < nai; i++)
  {
    if (active != nullptr && active[i]) continue;

    /* Calculate: T(a) %*% x */

    value = 0.;
    for (j = 0; j < neq; j++)
      value += AIMT(j,i)* xmat[j];

      /* Calculate: T(a) %*% x - b */

    ecart = value - bimat[i];

    /* Store the results */

    if (vmat != nullptr) vmat[ecr] = ecart;
    flag_active = (ecart < -eps);
    if (flag != nullptr) flag[ecr] = flag_active;
    if (flag_active) number++;
    ecr++;
  }
  return (number);
}

/*****************************************************************************/
/*!
 **  Minimize 1/2 t(x) %*% H %*% x + t(g) %*% x under the constraints
 **  t(Ae) %*% x = be and t(Ai) %*% x = bi
 **
 ** \return  Error return code
 **
 ** \param[in]     neq    number of equations in the system
 ** \param[in]     hmat   symmetric square matrix (Dimension: neq * neq)
 ** \param[in]     gmat   right-hand side vector (Dimension: neq)
 ** \param[in]     nae    Number of equalities
 ** \param[in]     aemat  matrix for equalities (Dimension: neq * nae)
 ** \param[in]     bemat  right-hand side for equalities (Dimension: nae)
 ** \param[in]     nai    Number of inequalities
 ** \param[in]     aimat  matrix for inequalities (Dimension: neq * nai)
 ** \param[in]     bimat  right-hand side for inequalities (Dimension: nai)
 **
 ** \param[in,out] xmat solution of the linear system with constraints (neq)
 **
 ** REMARKS:    The initial xmat has to satisfied all the constraints.
 **
 *****************************************************************************/
int matrix_qoci(int neq,
                double *hmat,
                double *gmat,
                int nae,
                double *aemat,
                double *bemat,
                int nai,
                double *aimat,
                double *bimat,
                double *xmat)
{
  int *active, error, sortie, namax, ncur, i, j, first, lec;
  double *beimat, *aeimat, *vmat, *lambda, *xcand, omega, omin, value;

  /* Initializations */

  error = 1;
  active = nullptr;
  lambda = beimat = aeimat = vmat = xcand = nullptr;
  namax = nae + nai;

  /* Case when there is no equality nor inequality constraints */

  if (namax <= 0)
  {
    if (matrix_qo(neq, hmat, gmat, xmat)) goto label_end;
    error = 0;
    goto label_end;
  }

  /* Core allocation */

  active = (int*) mem_alloc(sizeof(int) * nai, 1);
  xcand = (double*) mem_alloc(sizeof(double) * neq, 1);
  lambda = (double*) mem_alloc(sizeof(double) * namax, 1);
  vmat = (double*) mem_alloc(sizeof(double) * namax, 1);
  beimat = (double*) mem_alloc(sizeof(double) * namax, 1);
  aeimat = (double*) mem_alloc(sizeof(double) * neq * namax, 1);

  /* We first perform the optimization with equality constraints only */

  if (matrix_qoc(0, neq, hmat, gmat, nae, aemat, bemat, xcand, lambda))
    goto label_end;
  if (nai <= 0)
  {
    for (i = 0; i < neq; i++)
      xmat[i] = xcand[i];
    error = 0;
    return (0);
  }

  /* Evaluate the array active */

  if (st_calcul_error(neq, nai, NULL, aimat, bimat, xcand, NULL, active) == 0)
  {
    for (i = 0; i < neq; i++)
      xmat[i] = xcand[i];
    error = 0;
    goto label_end;
  }

  /* Implicit loop */

  sortie = 0;
  while (!sortie)
  {

    /* Construct the inequality matrices reduced to the active constraints */

    ncur = st_copy_active(nae, nai, active, neq, aemat, aimat, aeimat);
    ncur = st_copy_active(nae, nai, active, 1, bemat, bimat, beimat);
    if (matrix_qoc(1, neq, hmat, gmat, ncur, aeimat, beimat, xcand, lambda))
      goto label_end;

    if (st_calcul_error(neq, nai, active, aimat, bimat, xcand, vmat, NULL) == 0)
    {
      for (i = 0; i < neq; i++)
        xmat[i] = xcand[i];

      /* Look for the constraint that should not be used */

      first = -1;
      lec = nae;
      for (i = 0; i < nai; i++)
      {
        if (active[i] == 0) continue;
        active[i] = lambda[lec] >= 0;
        if (active[i]) first = i;
        lec++;
      }

      if (st_count_active(nai, active) == 0)
      {
        /* If no constraint has been used, end of the implicit loop */
        sortie = 1;
      }
      else
      {
        /* Otherwise, relax the first active constraint */

        active[first] = 0;
      }
    }
    else
    {

      /* Find an admissible solution between previous and new candidates */

      first = -1;
      omin = 1.e30;
      for (i = 0; i < nai; i++)
      {
        if (active[i]) continue;
        value = 0.;
        for (j = 0; j < neq; j++)
          value += AIMT(j,i)* (xcand[j] - xmat[j]);
        omega = vmat[i] / value;
        if (omega > omin) continue;
        first = i;
        omin = omega;
      }

      for (i = 0; i < neq; i++)
        xmat[i] += omin * (xcand[i] - xmat[i]);
      active[first] = 1;
    }
  }

  /* Set the error return code */

  error = 0;

  label_end: active = (int*) mem_free((char* ) active);
  xcand = (double*) mem_free((char* ) xcand);
  lambda = (double*) mem_free((char* ) lambda);
  vmat = (double*) mem_free((char* ) vmat);
  aeimat = (double*) mem_free((char* ) aeimat);
  beimat = (double*) mem_free((char* ) beimat);
  return (error);
}

/*****************************************************************************/
/*!
 **  Calculates the range of values within a rectangular matrix
 **
 ** \param[in]  n1 matrix dimension
 ** \param[in]  n2 matrix dimension
 ** \param[in]  v1 rectangular matrix (n1,n2)
 **
 ** \param[out] mini   Minimum value
 ** \param[out] maxi   Maximum value
 ** \param[out] norme1 Norm L1 of the matrix
 ** \param[out] norme2 Norm L2 of the matrix
 **
 *****************************************************************************/
void matrix_range(int n1,
                                  int n2,
                                  double *v1,
                                  double *mini,
                                  double *maxi,
                                  double *norme1,
                                  double *norme2)
{
  double value;
  int i1, i2;

  (*mini) = 1.e30;
  (*maxi) = -1.e30;
  (*norme1) = 0.;
  (*norme2) = 0.;

  for (i1 = 0; i1 < n1; i1++)
    for (i2 = 0; i2 < n2; i2++)
    {
      value = V1(i1, i2);
      if (value < (*mini)) (*mini) = value;
      if (value > (*maxi)) (*maxi) = value;
      (*norme1) += ABS(value);
      (*norme2) += value * value;
    }
  (*norme2) = sqrt(*norme2);
  return;
}

/*****************************************************************************/
/*!
 **  Concatenate two matrices
 **
 ** \return Pointer on the newly created concatenated matrix (or NULL)
 **
 ** \param[in]  mode Concatenation type:
 **                  1 : by row    (n31=n11+n21; n32=n12=n22)
 **                  2 : by column (n31=n11=n21; n32=n12+n22)
 ** \param[in]  n11  First dimension of the first matrix
 ** \param[in]  n12  Second dimension of the first matrix
 ** \param[in]  a1   Pointer to the first matrix
 ** \param[in]  n21  First dimension of the second matrix
 ** \param[in]  n22  Second dimension of the second matrix
 ** \param[in]  a2   Pointer to the second matrix
 **
 ** \param[out] n31  First dimension of the output matrix
 ** \param[out] n32  Second dimension of the output matrix
 **
 *****************************************************************************/
double* matrix_bind(int mode,
                                    int n11,
                                    int n12,
                                    double *a1,
                                    int n21,
                                    int n22,
                                    double *a2,
                                    int *n31,
                                    int *n32)
{
  double *a, *v1, *v2;
  int error, n1, n2, n3, i, j, neq;

  /* Initializations */

  error = 1;
  (*n31) = (*n32) = 0;
  a = nullptr;

  /* Preliminary tests */

  if (mode == 1)
  {
    if (n12 != n22)
    {
      messerr("Binding by row: Input matrices must share same column number");
      goto label_end;
    }
    (*n31) = n1 = n11 + n21;
    (*n32) = n2 = n12;
  }
  else if (mode == 2)
  {
    if (n11 != n21)
    {
      messerr("Binding by column: Input matrices must share same row number");
      goto label_end;
    }
    (*n31) = n1 = n11;
    (*n32) = n2 = n12 + n22;
  }
  else
  {
    messerr("The concatenation mode must be 1 or 2");
    goto label_end;
  }

  /* Core allocation */

  a = (double*) mem_alloc(sizeof(double) * (*n31) * (*n32), 0);
  if (a == nullptr) goto label_end;

  /* Copy the first matrix */

  v1 = a1;
  n1 = n11;
  n2 = n12;
  neq = (*n31);
  for (i = 0; i < n1; i++)
    for (j = 0; j < n2; j++)
      A(i,j)= V1(i,j);

      /* Concatenate the second matrix */

  v2 = a2;
  n2 = n21;
  n3 = n22;
  if (mode == 1)
  {
    // By row

    for (i = 0; i < n2; i++)
      for (j = 0; j < n3; j++)
        A(n11+i,j)= V2(i,j);
      }
      else
      {
        // By column

        for (i=0; i<n2; i++)
        for (j=0; j<n3; j++)
        A(i,n12+j) = V2(i,j);
      }

      /* Error return code */

  error = 0;

  label_end: if (error) a = (double*) mem_free((char* ) a);
  return (a);
}

/*****************************************************************************/
/*!
 **  Find the minimum or maximum value within an array
 **
 ** \return Rank of the minimum or maximum value
 **
 ** \param[in]  mode 1 for minimum and 2 for maximum
 ** \param[in]  ntab Number of samples in the array
 ** \param[in]  tab  Array of values
 **
 *****************************************************************************/
int matrix_get_extreme(int mode, int ntab, double *tab)
{
  int i, ibest;
  double vbest;

  /* Dispatch */

  if (mode == 1)
  {

    // Find the minimum

    ibest = -1;
    vbest = 1.e30;
    for (i = 0; i < ntab; i++)
    {
      if (FFFF(tab[i])) continue;
      if (tab[i] > vbest) continue;
      vbest = tab[i];
      ibest = i;
    }
  }
  else
  {

    // Find the maximum

    ibest = -1;
    vbest = -1.e30;
    for (i = 0; i < ntab; i++)
    {
      if (FFFF(tab[i])) continue;
      if (tab[i] < vbest) continue;
      vbest = tab[i];
      ibest = i;
    }
  }
  return (ibest);
}

/*****************************************************************************/
/*!
 **  Invert a symmetric square matrix
 **  Pivots are assumed to located on the diagonal
 **
 ** \return  Return code: 0 no error; k if the k-th pivot is zero
 **
 ** \param[in]  a    input matrix
 ** \param[in]  neq  number of equations in the matrix 'a'
 **
 ** \param[out] b    output matrix
 **
 ** \remark  The difference with matrix_invert() is that the output
 ** \remark  matrix is different from input matrix
 **
 *****************************************************************************/
int matrix_invreal_copy(const double *a, int neq, double *b)
{
  int i, error;

  /* Copy the input matrix into the output matrix */

  for (i = 0; i < neq * neq; i++)
    b[i] = a[i];

  /* Invert the matrix */

  error = matrix_invreal(b, neq);

  return (error);
}

/*****************************************************************************/
/*!
 **  Invert a square real full matrix
 **
 ** \return  Error return code
 **
 ** \param[in]  mat  input matrix, destroyed in computation and replaced by
 **                  resultant inverse
 ** \param[in]  neq  number of equations in the matrix 'a'
 **
 *****************************************************************************/
int matrix_invreal(double *mat, int neq)
{
  double *cofac, det;
  int error;

  /* Initialization */

  error = 1;
  cofac = nullptr;

  /* Calculate the determinant */

  det = matrix_determinant(neq, mat);
  //if (ABS(det) < _getEpsSVD()) return (1);
  if (ABS(det) < 1.e-12) return (1);
  if (std::isnan(det))
  {
    print_matrix("Mat", 0, 1, neq, neq, NULL, mat);
    messageAbort("Values NAN found in matrix");
  }

  if (neq > 1)
  {

    /* Core allocation */

    cofac = (double*) mem_alloc(sizeof(double) * neq * neq, 0);
    if (cofac == nullptr) goto label_end;

    /* Calculate the cofactor */

    if (matrix_cofactor(neq, mat, cofac)) goto label_end;

    /* Transpose the cofactor to obtain the adjoint matrix */

    matrix_transpose(neq, neq, cofac, mat);
  }

  /* Final normation */

  for (int i = 0; i < neq * neq; i++)
    mat[i] /= det;

  /* Set the error return code */

  error = 0;

  label_end: cofac = (double*) mem_free((char* ) cofac);
  return (error);
}

/*****************************************************************************/
/*!
 **  Set a square matrix to Identity
 **
 ** \param[in]  neq  dimension of the matrix
 ** \param[in]  a    Square matrix
 **
 *****************************************************************************/
void matrix_set_identity(int neq, double *a)
{
  for (int i = 0; i < neq; i++)
    for (int j = 0; j < neq; j++)
      A(i,j)= (i == j) ? 1 : 0;
    }

    /*****************************************************************************/
    /*!
     **  Invert a symmetric square matrix (stored as triangular)
     **
     ** \return  Error returned code
     **
     ** \param[in,out] tl input matrix, destroyed in computation and replaced by
     **                   resultant inverse
     ** \param[in]  neq  number of equations in the matrix 'a'
     ** \param[in]  rank Type of message when inversion problem is encountered
     **                  >=0: message involves 'rank+1'
     **                  -1:  neutral message
     **                  -2:  no message
     **
     ** \remark  It is unnecessary to edit a message if inversion problem occurs
     **
     *****************************************************************************/
int matrix_invert_triangle(int neq, double *tl, int rank)
{
  double *a;
  int error;

  a = (double*) mem_alloc(sizeof(double) * neq * neq, 1);

  matrix_tri2sq(neq, tl, a);
  error = matrix_invert(a, neq, rank);
  matrix_square_to_triangle(0, neq, a, tl);

  a = (double*) mem_free((char* ) a);
  return error;
}

/*****************************************************************************/
/*!
 **  Change the storage of triangular matrices (stored columnwise)
 **
 ** \param[in]  neq  Matrix dimension
 ** \param[in]  tl   Lower Triangular Matrix
 ** \param[in]  tu   Upper Triangular Matrix
 **
 *****************************************************************************/
void matrix_tl2tu(int neq, const double *tl, double *tu)
{
  for (int i = 0; i < neq; i++)
    for (int j = 0; j <= i; j++)
      TU(i,j)= TL(i,j);
    }
