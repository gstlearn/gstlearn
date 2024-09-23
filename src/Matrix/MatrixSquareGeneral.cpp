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
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Basic/VectorHelper.hpp"

#define SQ(i,j,neq)   ((j) * neq + (i))
#define A(i,j)         a[SQ(i,j,neq)]
#define B(i,j)         b[SQ(i,j,neq)]
#define C(i,j)         c[SQ(i,j,neqm1)]

MatrixSquareGeneral::MatrixSquareGeneral(int nrow)
  : AMatrixSquare(nrow)
{
}

MatrixSquareGeneral::MatrixSquareGeneral(const MatrixSquareGeneral &r) 
  : AMatrixSquare(r)
{
}

MatrixSquareGeneral::MatrixSquareGeneral(const AMatrix &m)
  : AMatrixSquare(m)
{
  if (!m.isSquare())
  {
    messerr("The input matrix should be Square");
    _clear();
    return;
  }
}

MatrixSquareGeneral& MatrixSquareGeneral::operator= (const MatrixSquareGeneral &r)
{
  if (this != &r)
  {
    AMatrixSquare::operator=(r);
  }
  return *this;
}

MatrixSquareGeneral::~MatrixSquareGeneral()
{
}

/**
 * Converts a VectorVectorDouble into a Matrix
 * Note: the input argument is stored by row (if coming from [] specification)
 * @param X Input VectorVectorDouble argument
 * @return The returned square matrix
 *
 * @remark: the matrix is transposed implicitly while reading
 */
MatrixSquareGeneral* MatrixSquareGeneral::createFromVVD(const VectorVectorDouble& X)
{
  int nrow = (int) X.size();
  int ncol = (int) X[0].size();
  if (nrow != ncol)
  {
    messerr("The matrix does not seem to be square");
    return nullptr;
  }

  MatrixSquareGeneral* mat = new MatrixSquareGeneral(nrow);
  mat->_fillFromVVD(X);
  return mat;
}

MatrixSquareGeneral* MatrixSquareGeneral::createFromVD(const VectorDouble &X,
                                                       int nrow,
                                                       bool byCol,
                                                       bool invertColumnOrder)
{
  int ncol = nrow;
  if (nrow * ncol != (int) X.size())
  {
    messerr("Inconsistency between arguments 'nrow'(%d) and 'ncol'(%d)", nrow, ncol);
    messerr("and the dimension of the input Vector (%d)", (int) X.size());
  }
  MatrixSquareGeneral* mat = new MatrixSquareGeneral(nrow);

  int lec = 0;
  if (byCol)
  {
    for (int irow = 0; irow < nrow; irow++)
      for (int icol = 0; icol < ncol; icol++)
      {
        int jcol = (invertColumnOrder) ? ncol - icol - 1 : icol;
        mat->setValue(irow, jcol, X[lec++]);
      }
  }
  else
  {
    for (int icol = 0; icol < ncol; icol++)
      for (int irow = 0; irow < nrow; irow++)
      {
        int jcol = (invertColumnOrder) ? ncol - icol - 1 : icol;
        mat->setValue(irow, jcol, X[lec++]);
      }
  }
  return mat;
}

MatrixSquareGeneral* prodNormMatMat(const AMatrixDense* a,
                                    const AMatrixDense* m,
                                    bool transpose)
{
  int nrow = (transpose) ? a->getNCols() : a->getNRows();
  MatrixSquareGeneral *mat = new MatrixSquareGeneral(nrow);
  mat->prodNormMatMatInPlace(a, m, transpose);
  return mat;
}

MatrixSquareGeneral* prodNormMat(const AMatrixDense &a, const VectorDouble& vec, bool transpose)
{
  int nsym = (transpose) ? a.getNCols() : a.getNRows();
  MatrixSquareGeneral *mat = new MatrixSquareGeneral(nsym);
  mat->prodNormMatInPlace(a, vec, transpose);
  return mat;
}

int MatrixSquareGeneral::_invertLU()
{
  int neq = getNRows();

  // Perform the LU decomposition
  MatrixSquareGeneral tls(neq);
  MatrixSquareGeneral tus(neq);
  MatrixSquareGeneral ais(neq);
  ais.fill(0.);
  if (decomposeLU(tls, tus)) return 1;

  VectorDouble b(neq);
  VectorDouble x(neq);
  for (int i = 0; i < neq; i++)
  {
    // Preparing the right-hand side vector (column of the identity matrix)
    VH::fill(b, 0.);
    b[i] = 1.;

    if (_solveLU(tus, tls, b.data(), x.data())) return 1;

    for (int j = 0; j < neq; j++)
      ais.setValue(i, j, x[j]);
  }

  // Copy the inverse matrix in the input matrix
  for (int irow = 0; irow < neq; irow++)
    for (int icol = 0; icol < neq; icol++)
      setValue(irow, icol, ais.getValue(irow, icol));
  return 0;
}

int MatrixSquareGeneral::_solveLU(const MatrixSquareGeneral& tus,
                                  const MatrixSquareGeneral& tls,
                                  const double *b,
                                  double *x)
{
  int neq = getNRows();
  VectorDouble y(neq);
  if (_forwardLU(tls, b, y.data())) return 1;
  if (_backwardLU(tus, y.data(), x)) return 1;
  return 0;
}

/*****************************************************************************/
/*!
 **  Get the solution of a linear system (after LU decomposition)
 **     L * x = b
 **
 ** \return  Error returned code
 **
 ** \param[in]  tls  Square matrix containing lower triangle (stored column-wise)
 ** \param[in]  b    matrix (dimension neq)
 ** \param[in]  eps  Tolerance
 **
 ** \param[out] x    resulting matrix (dimension neq)
 **
 *****************************************************************************/
int MatrixSquareGeneral::_forwardLU(const MatrixSquareGeneral& tls, const double *b, double *x, double eps)
{
  int neq = getNRows();
  for (int i = 0; i < neq; i++) x[i] = 0.;
  for (int i = 0; i < neq; i++)
  {
    double tmp = b[i];
    for (int j = 0; j < i; j++)
      tmp -= tls.getValue(i,j) * x[j];

    double pivot = tls.getValue(i,i);
    if (abs(pivot) < eps) return 1;
    x[i] = tmp / pivot;
  }
  return 0;
}

/*****************************************************************************/
/*!
 **  Get the solution of a linear system (after LU decomposition)
 **     U * x = b
 **
 ** \return  Error returned code
 **
 ** \param[in]  tus  square matrix containing upper triangle (stored line-wise)
 ** \param[in]  b    matrix (dimension neq)
 ** \param[in]  eps  Tolerance
 **
 ** \param[out] x    resulting matrix (dimension neq)
 **
 ** \remark As the Upper matrix is stored line-wise, it is considered
 ** \remark as the transposed version of a lower triangle
 **
 *****************************************************************************/
int MatrixSquareGeneral::_backwardLU(const MatrixSquareGeneral& tus, const double *b, double *x, double eps)
{
  int neq = getNRows();
  for (int i = neq-1; i >= 0; i--)
  {
    double tmp = b[i];
    for (int j = i+1; j < neq; j++)
      tmp -= tus.getValue(i,j) * x[j];

    double pivot = tus.getValue(i,i);
    if (abs(pivot) < eps) return 1;
    x[i] = tmp / pivot;
  }
  return 0;
}

/**
 * LU Decomposition of a square matrix (not necessarily symmetric)
 * @param tls  Output square matrix containing lower triangle (stored columnwise)
 * @param tus  Output square matrix containing upper triangle (stored linewise)
 * @param eps  Tolerance
 *
 * @remarks The output matrices 'tus'  and 'tls' must be dimensioned beforehand
 */
int MatrixSquareGeneral::decomposeLU(MatrixSquareGeneral& tls,
                                     MatrixSquareGeneral& tus,
                                     double eps)
{
  int neq = getNRows();
  tls.fill(0.);
  tus.fill(0.);

  for (int i = 0; i < neq; i++)
    tls.setValue(i, i, 1.);

  for (int i = 0; i < neq; i++)
  {
    int ip1 = i + 1;
    int im1 = i - 1;

    for (int j = 0; j < neq; j++)
    {
      tus.setValue(i,j,getValue(i,j));
      if (im1 >= 0)
      {
        for (int k = 0; k <= im1; k++)
        {
          tus.setValue(i,j, tus.getValue(i,j) - tls.getValue(i,k) * tus.getValue(k,j));
        }
      }
    }
    if (ip1 < neq)
    {
      for (int j = ip1; j < neq; j++)
      {
        tls.setValue(j,i,getValue(j,i));
        if (im1 >= 0)
        {
          for (int k = 0; k <= im1; k++)
          {
            tls.setValue(j,i, tls.getValue(j,i) - tls.getValue(j,k) * tus.getValue(k,i));
          }
        }
        double pivot = tus.getValue(i,i);
        if (abs(pivot) < eps) return 1;
        tls.setValue(j,i, tls.getValue(j,i) / pivot);
      }
    }
  }
  return 0;
}

