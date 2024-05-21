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

#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/AMatrixSquare.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/AException.hpp"

#define TRI(i)        (((i) * ((i) + 1)) / 2)
#define SQ(i,j,neq)   ((j) * neq + (i))
#define AT(i,j)        at[TRI(j)+(i)] /* for j >= i */
#define AL(i,j)        al[SQ(i,j,neq)-TRI(j)] /* for i >= j */
#define BS(i,j)        b[SQ(i,j,neq)] // Proposition a valider: c'etait j,i
#define XS(i,j)        x[SQ(i,j,neq)] // Proposition a valider: c'etait j,i
#define AS(i,j)        a[SQ(i,j,neq)]
#define TL(i,j)        tl[SQ(i,j,neq)-TRI(j)] /* for i >= j */

#define _TL(i,j)       _tl[SQ(i,j,neq)-TRI(j)] /* for i >= j */
#define _XL(i,j)       _xl[SQ(i,j,neq)-TRI(j)] /* for i >= j */

#define HA(i,j)        ha[SQ(i,j,neq)]

MatrixSquareSymmetric::MatrixSquareSymmetric(int nrow, int opt_eigen)
    : AMatrixSquare(nrow, opt_eigen),
      _squareSymMatrix(),
      _flagCholeskyDecompose(false),
      _flagCholeskyInverse(false),
      _tl(),
      _xl(),
      _factor()
{
  _allocate_();
}

MatrixSquareSymmetric::MatrixSquareSymmetric(const MatrixSquareSymmetric &r) 
  : AMatrixSquare(r),
   _squareSymMatrix(),
   _flagCholeskyDecompose(r._flagCholeskyDecompose),
   _flagCholeskyInverse(r._flagCholeskyInverse),
   _tl(),
   _xl(),
   _factor()
{
  _recopy(r);
}

MatrixSquareSymmetric::MatrixSquareSymmetric(const AMatrix &m)
    : AMatrixSquare(m),
      _squareSymMatrix(),
      _flagCholeskyDecompose(false),
      _flagCholeskyInverse(false),
      _tl(),
      _xl(),
      _factor()
{
  // TODO: clean this code or move it upwards
  if (!m.isSymmetric())
  {
    messerr("The input matrix should be Symmetric");
    _clear();
    return;
  }
  const MatrixSquareSymmetric* matrixLoc = dynamic_cast<const MatrixSquareSymmetric*>(&m);
  if (matrixLoc != nullptr)
    _recopy(*matrixLoc);
  else
  {
    _allocate_();
    AMatrix::copyElements(m);
  }
}

MatrixSquareSymmetric& MatrixSquareSymmetric::operator= (const MatrixSquareSymmetric &r)
{
  if (this != &r)
  {
    AMatrixSquare::operator=(r);
    _recopy(r);
  }
  return *this;
}

MatrixSquareSymmetric::~MatrixSquareSymmetric()
{
}

/**
 * Converts a VectorVectorDouble into a Square Symmetric Matrix
 * Note: the input argument is stored by row (if coming from [] specification)
 * @param  X Input VectorVectorDouble argument
 * @param opt_eigen Option for use of Eigen Library
 * @return The returned square symmetric matrix
 *
 * @remark: the matrix is transposed implicitly while reading
 */
MatrixSquareSymmetric* MatrixSquareSymmetric::createFromVVD(const VectorVectorDouble& X, int opt_eigen)
{
  int nrow = (int) X.size();
  int ncol = (int) X[0].size();
  if (nrow != ncol)
  {
    messerr("The matrix does not seem to be square");
    return nullptr;
  }
  MatrixSquareSymmetric* mat = new MatrixSquareSymmetric(nrow, opt_eigen);
  mat->_fillFromVVD(X);
  return mat;
}

MatrixSquareSymmetric* MatrixSquareSymmetric::createFromVD(const VectorDouble &X,
                                                           int nrow,
                                                           int opt_eigen)
{
  int ncol = nrow;
  if (nrow * ncol != (int) X.size())
  {
    messerr("Inconsistency between arguments 'nrow'(%d) and 'ncol'(%d)", nrow, ncol);
    messerr("and the dimension of the input Vector (%d)", (int) X.size());
  }
  // Check symmetry
  MatrixRectangular* mattemp = MatrixRectangular::createFromVD(X, nrow, ncol);
  if (! mattemp->isSymmetric())
  {
    messerr("The input matrix does not seem to be Square and symmetric");
    delete mattemp;
    return nullptr;
  }
  delete mattemp;

  MatrixSquareSymmetric *mat = new MatrixSquareSymmetric(nrow, opt_eigen);

  int lec = 0;
  for (int irow = 0; irow < nrow; irow++)
    for (int icol = 0; icol < ncol; icol++)
      mat->setValue(irow, icol, X[lec++]);
  return mat;
}

double MatrixSquareSymmetric::_getValueByRank_(int irank) const
{
  return _squareSymMatrix[irank];
}

double& MatrixSquareSymmetric::_getValueRef_(int irow, int icol)
{
  int rank = _getIndexToRank_(irow, icol);
  return _squareSymMatrix[rank];
}

void MatrixSquareSymmetric::_setValueByRank_(int irank, double value)
{
  _squareSymMatrix[irank] = value;
}

void MatrixSquareSymmetric::_prodMatVecInPlacePtr_(const double *x, double *y, bool transpose) const
{
  _matrix_triangular_product(getNRows(),2,_squareSymMatrix.data(),x,y);
}

void MatrixSquareSymmetric::_prodVecMatInPlacePtr_(const double *x, double *y, bool transpose) const
{
  _matrix_triangular_product(getNRows(),2,_squareSymMatrix.data(),x,y);
}

/**
 * \warning : values is provided as a square complete matrix
 */
void MatrixSquareSymmetric::_setValues(const double* values, bool byCol)
{
  // Check that the input argument corresponds to a square symmetric matrix
  for (int icol = 0; icol < getNCols(); icol++)
    for (int irow = 0; irow < getNRows(); irow++)
    {
      double val1 = values[icol * getNRows() + irow];
      double val2 = values[irow * getNCols() + icol];
      if (ABS(val1 - val2) > EPSILON10)
      {
        messerr(
            "Argument 'values' must correspond to a Square Symmetric Matrix");
        messerr("- Element[%d,%d] = %lf", icol, irow, val1);
        messerr("- Element(%d,%d) = %lf", irow, icol, val2);
        messerr("Operation is aborted");
        return;
      }
    }

  AMatrixDense::_setValues(values, byCol);
}

int MatrixSquareSymmetric::_invert()
{
  if (isFlagEigen())
    return AMatrixDense::_invert();
  else
    return _matrix_invert_triangle(getNRows(),_squareSymMatrix.data());
}

void MatrixSquareSymmetric::_allocate_()
{
  _squareSymMatrix.resize(_getMatrixPhysicalSize_());
}

int MatrixSquareSymmetric::_getIndexToRank_(int irow, int icol) const
{
  int n = getNRows();
  if (irow >= icol)
    return (icol * n + irow - icol * (icol + 1) / 2);
  else
    return (irow * n + icol - irow * (irow + 1) / 2);
}

int MatrixSquareSymmetric::_getMatrixPhysicalSize_() const
{
  int n = getNRows();
  return (n * (n + 1) / 2);
}

int MatrixSquareSymmetric::_solve(const VectorDouble& b, VectorDouble& x) const
{
  if (isFlagEigen())
    return AMatrixDense::_solve(b, x);
  else
  {
    int size = (int) b.size();
    VectorDouble alocal = _squareSymMatrix;
    VectorDouble blocal = b;
    int pivot = _matrix_solve(alocal, blocal, x, size, 1);
    return (pivot != 0);
  }
}

bool MatrixSquareSymmetric::_isPhysicallyPresent(int irow, int icol) const
{
  if (icol >  irow) return false;
  return true;
}

/**
 * Perform the product: this = t(Y) %*% X %*% Y (T=false) or Y % X %*% t(Y) (T=true)
 * @param y: Matrix (possibly rectangular)
 * @param x: Square matrix (optional)
 * @param transpose: transposition flag (T in the description)
 * \remarks The number of rows of Y must be equal to the dimension of X
 * \remarks The output matrix is square with dimension equal to the number of columns of Y
 */
void MatrixSquareSymmetric::normMatrix(const AMatrix& y, const AMatrixSquare& x, bool transpose)
{
  bool xEmpty = x.empty();
  int n = 0;

  if (xEmpty)
  {
    if (transpose)
    {
      if (getNSize() != y.getNRows())
        my_throw("Incompatible matrix dimensions: y.nrows != this.size");
      n = y.getNCols();
    }
    else
    {
      if (getNSize() != y.getNCols())
        my_throw("Incompatible matrix dimensions: y.ncols != this.size");
      n = y.getNRows();
    }
  }
  else
  {
    if (transpose)
    {
      if (y.getNCols() != x.getNSize())
        my_throw("Incompatible matrix dimensions: y.ncols != x.nsize");
      n = x.getNSize();
    }
    else
    {
      if (y.getNRows() != x.getNSize())
        my_throw("Incompatible matrix dimensions: y.nrows != x.nsize");
      n = x.getNSize();
    }
  }

  int nout = getNSize();
  for (int irow = 0; irow < nout; irow++)
    for (int icol = 0; icol <= irow; icol++)
    {
      double value = 0.;

      if (xEmpty)
      {
        if (! transpose)
        {
          for (int k = 0; k < n; k++)
            value += y.getValue(k,irow) * y.getValue(k,icol);
        }
        else
        {
          for (int k = 0; k < n; k++)
            value += y.getValue(irow,k) * y.getValue(icol,k);
        }
      }
      else
      {
        if (!transpose)
        {
          for (int k = 0; k < n; k++)
            for (int l = 0; l < n; l++)
              value += y.getValue(k, irow) * x.getValue(k, l) * y.getValue(l, icol);
        }
        else
        {
          for (int k = 0; k < n; k++)
            for (int l = 0; l < n; l++)
              value += y.getValue(irow, k) * x.getValue(k, l) * y.getValue(icol, l);
        }
      }

      setValue(irow,icol,value);
    }
}

int MatrixSquareSymmetric::_terminateEigen(const VectorDouble &eigenValues,
                                           const VectorDouble &eigenVectors,
                                           bool optionPositive,
                                           bool changeOrder)
{
  int nrows = getNRows();

  _eigenValues = eigenValues;

  if (_eigenVectors != nullptr) delete _eigenVectors;

  if (changeOrder)
    std::reverse(_eigenValues.begin(), _eigenValues.end());

  _eigenVectors = MatrixSquareGeneral::createFromVD(eigenVectors, nrows, false,
                                                    0, changeOrder);

  if (optionPositive) _eigenVectors->makePositiveColumn();

  _flagEigenDecompose = true;

  return 0;
}

int MatrixSquareSymmetric::computeEigen(bool optionPositive)
{
  if (isFlagEigen())
  {
    return AMatrixDense::_computeEigen(optionPositive);
  }
  int nrows = getNRows();
  VectorDouble eigenValues(nrows, 0.);
  VectorDouble eigenVectors(nrows * nrows, 0);

  if (matrix_eigen(this->getValues().data(), nrows,
                   eigenValues.data(),
                   eigenVectors.data())) return 1;

  return _terminateEigen(eigenValues, eigenVectors, optionPositive, false);
}

int MatrixSquareSymmetric::computeGeneralizedEigen(const MatrixSquareSymmetric& b, bool optionPositive)
{
  if (isFlagEigen())
  {
    return AMatrixDense::_computeGeneralizedEigen(b, optionPositive);
  }
  int nrows = getNRows();
  VectorDouble eigenValues(nrows, 0.);
  VectorDouble eigenVectors(nrows * nrows, 0);

  if (_matrix_geigen(this->getValues().data(), b.getValues().data(), nrows,
                     eigenValues.data(), eigenVectors.data())) return 1;

  return _terminateEigen(eigenValues, eigenVectors, optionPositive, true);
}

/// =============================================================================
/// The subsequent methods rely on the specific local storage ('squareSymMatrix')
/// =============================================================================

void MatrixSquareSymmetric::_recopy(const MatrixSquareSymmetric& r)
{
  _squareSymMatrix = r._squareSymMatrix;
  _tl = r._tl;
  _xl = r._xl;
  _flagCholeskyDecompose = r._flagCholeskyDecompose;
  _flagCholeskyInverse   = r._flagCholeskyInverse;
  _flagEigenDecompose    = r._flagEigenDecompose;
  _factor                = r._factor;
}

double MatrixSquareSymmetric::_getValue(int irow, int icol) const
{
  int rank = _getIndexToRank_(irow,icol);
  return _squareSymMatrix[rank];
}

void MatrixSquareSymmetric::_setValue(int irow, int icol, double value)
{
  int irank = _getIndexToRank_(irow, icol);
  _squareSymMatrix[irank] = value;
}

void MatrixSquareSymmetric::_updValue(int irow, int icol, const EOperator& oper, double value)
{
  int irank = _getIndexToRank_(irow, icol);
  _squareSymMatrix[irank] = modifyOperator(oper, _squareSymMatrix[irank], value);
}

/*****************************************************************************/
/*!
 **  Calculates the generalized eigen value problem
 **           A X = B X l
 **
 ** \return  Return code:
 ** \return   0 no error
 ** \return   1 convergence problem
 **
 ** \param[in]  a     square symmetric matrix (dimension = neq * neq)
 ** \param[in]  b     square symmetric matrix (dimension = neq * neq)
 ** \param[in]  neq   matrix dimension
 **
 ** \param[out] value  matrix of the eigen values (dimension: neq)
 ** \param[out] vector matrix of the eigen vectors (dimension: neq*neq)
 **
 *****************************************************************************/
int MatrixSquareSymmetric::_matrix_geigen(const double *a,
                                          const double *b,
                                          int neq,
                                          double *value,
                                          double *vector) const
{
  // Compute eigen decomposition of B
  VectorDouble LB(neq);
  VectorDouble PhiB(neq * neq);
  if (matrix_eigen(b, neq, LB.data(), PhiB.data())) return 1;

  // Compute auxiliary terms
  VectorDouble PhiBm = PhiB;
  int ecr = 0;
  for (int i = 0; i < neq; i++)
    for (int j = 0; j < neq; j++, ecr++)
      PhiBm[ecr] /= sqrt(LB[i]);

  VectorDouble Am(neq * neq, 0.);
  for (int i = 0; i < neq; i++)
    for (int j = 0; j < neq; j++)
      for (int k = 0; k < neq; k++)
        for (int l = 0; l < neq; l++)
          Am[i * neq + j] += PhiBm[i * neq + k] * a[l * neq + k] * PhiBm[j * neq + l];

  // Compute eigen decomposition of Am
  VectorDouble LA(neq);
  VectorDouble PhiA(neq * neq);
  if (matrix_eigen(Am.data(), neq, LA.data(), PhiA.data())) return 1;

  VectorDouble Phi(neq * neq,0.);
  ecr = 0;
  for (int i = 0; i < neq; i++)
    for (int j = 0; j < neq; j++)
      for (int k = 0; k < neq; k++)
        Phi[j * neq + i] += PhiBm[k * neq + i] * PhiA[j * neq + k];

  // Sort the eigen values by increasing values
  VectorInt ranks = VH::sortRanks(LA, true, neq);

  // Ultimate assignments
  for (int i = 0; i < neq; i++)
    value[i] = LA[ranks[i]];

  for (int i = 0; i < neq; i++)
    for (int j = 0; j < neq; j++)
      vector[i * neq + j] = -Phi[ranks[i] * neq + j];

  return 0;
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
void MatrixSquareSymmetric::_matrix_triangular_product(int neq,
                                                       int mode,
                                                       const double *al,
                                                       const double *b,
                                                       double *x) const
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
        value += AT(j,i) * b[j];
      for (j = i + 1; j < neq; j++)
        value += AT(i,j) * b[j];
      x[i] = value;
    }
  }
  else
  {
    for (i = 0; i < neq; i++)
    {
      value = 0.;
      for (j = 0; j <= i; j++)
        value += AL(i,j) * b[j];
      for (j = i + 1; j < neq; j++)
        value += AL(j,i) * b[j];
      x[i] = value;
    }
  }
  return;
}

/****************************************************************************/
/*!
 **  Check if a matrix is definite positive
 **
 ** \return  True if the matrix is definite positive; False otherwise
 **
 *****************************************************************************/
bool MatrixSquareSymmetric::isDefinitePositive()
{
  /* Calculate the eigen values and vectors */

  if (computeEigen()) messageAbort("matrix_eigen");

  // Get the Eigen values

  VectorDouble valpro = getEigenValues();

  /* Check if the eigen values are all positive */

  for (int i = 0, n = (int) valpro.size(); i < n; i++)
  {
    if (valpro[i] < -1.0e-10)
    {
      messerr("The matrix is not definite positive: Eigen value #%d = %lf",
              i + 1, valpro[i]);
      return false;
    }
  }
  return true;
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
 ** \param[in]  eps  tolerance
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
int MatrixSquareSymmetric::_matrix_solve(VectorDouble& at,
                                         VectorDouble& b,
                                         VectorDouble& x,
                                         int neq,
                                         int nrhs,
                                         double eps) const

{
  double pivot, ratio;

  for (int k = 0; k < neq - 1; k++)
  {
    pivot = AT(k, k);
    if (ABS(pivot) < eps) return (k + 1);
    for (int i = k + 1; i < neq; i++)
    {
      ratio = AT(k,i)/ pivot;
      for (int j=i; j<neq; j++)  AT(i,j) -= AT(k,j) * ratio;
      for (int l=0; l<nrhs; l++) BS(i,l) -= BS(k,l) * ratio;
    }
  }

  pivot = AT(neq - 1, neq - 1);
  if (ABS(pivot) < eps) return (neq);

  for (int l = 0; l < nrhs; l++)
    XS(neq-1,l) = BS(neq-1,l) / pivot;

  for (int l = 0; l < nrhs; l++)
  {
    for (int k = neq - 2; k >= 0; k--)
    {
      ratio = BS(k, l);
      for (int j = k + 1; j < neq; j++)
        ratio -= AT(k,j) * XS(j,l);
      XS(k,l) = ratio / AT(k,k);
    }
  }
  return (0);
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
 **
 ** \remark  It is unnecessary to edit a message if inversion problem occurs
 **
 *****************************************************************************/
int MatrixSquareSymmetric::_matrix_invert_triangle(int neq, double *tl)
{
  VectorDouble a(neq * neq);
  _matrix_tri2sq(neq, tl, a.data());
  if (matrix_invert(a.data(), neq, -1)) return 1;
  _matrix_sq2tri(0, neq, a.data(), tl);
  return 0;
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
void MatrixSquareSymmetric::_matrix_tri2sq(int neq, const double *tl, double *a)
{
  for (int i = 0; i < neq; i++)
    for (int j = 0; j < neq; j++)
    {
      AS(i,j)= (j < i) ? TL(i,j) : TL(j,i);
    }
  }

/*****************************************************************************/
/*!
 **  Transform a square symmetric matrix into a triangular one
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
void MatrixSquareSymmetric::_matrix_sq2tri(int mode, int neq, const double *a, double *tl)
{
  for (int i = 0; i < neq; i++)
    for (int j = 0; j < neq; j++)
    {
      if (mode == 0)
      {
        if (j <= i) TL(i,j) = AS(i,j);
      }
      else
      {
        if (j >= i) TL(j,i) = AS(i,j);
      }
    }
  }

/*****************************************************************************/
/*!
 **  Create the Symmetric matrix as the product of 'tl' (lower triangle) by its transpose
 **
 ** \param[in]  neq    Number of rows or columns in the system
 ** \param[in]  tl     Lower triangular matrix defined by column (Dimension; neq*(neq+1)/2)
 ** \param[in]  opt_eigen Option for use of Eigen Library
 **
 *****************************************************************************/
MatrixSquareSymmetric* MatrixSquareSymmetric::createFromTLTU(int neq,
                                                             const VectorDouble &tl,
                                                             int opt_eigen)
{
  MatrixSquareSymmetric *mat = new MatrixSquareSymmetric(neq, opt_eigen);

  for (int i = 0; i < neq; i++)
    for (int j = 0; j < neq; j++)
    {
      double value = 0.;
      for (int k = 0; k < neq; k++)
      {
        if (k > i || k > j) continue;
        value += TL(i,k) * TL(j,k);
      }
      mat->setValue(i, j, value);
    }
  return mat;
}

/*****************************************************************************/
/*!
 **  Fill a square matrix with a triangular matrix
 **
 ** \param[in]  mode   0: TL (upper); 1: TL (lower)
 ** \param[in]  neq    number of equations in the system
 ** \param[in]  tl     Triangular matrix (lower part)
 ** \param[in]  opt_eigen Option for use of Eigen Library
 **
 *****************************************************************************/
MatrixSquareSymmetric* MatrixSquareSymmetric::createFromTriangle(int mode,
                                                                 int neq,
                                                                 const VectorDouble &tl,
                                                                 int opt_eigen)
{
  MatrixSquareSymmetric *mat = new MatrixSquareSymmetric(neq, opt_eigen);

  mat->fill(0.);

  for (int i = 0; i < neq; i++)
    for (int j = 0; j < neq; j++)
    {
      if (mode == 0)
      {
        if (j <= i) mat->setValue(i,j,TL(i,j));
      }
      else
      {
        if (j >= i) mat->setValue(i,j,TL(j,i));
      }
    }
  return mat;
}

int MatrixSquareSymmetric::getTriangleSize() const
{
  int neq = getNRows();
  int size = neq * (neq + 1) / 2;
  return size;
}

/*****************************************************************************/
/*!
 **  Performs the Cholesky triangular decomposition of a definite
 **  positive symmetric matrix
 **         A = t(TL) * TL
 **
 ** \return  Error return code
 **
 *****************************************************************************/
int MatrixSquareSymmetric::computeCholesky()
{
  _flagCholeskyDecompose = false;
  if (isFlagEigen())
  {
    _factor = _eigenMatrix.llt();
    int neq = getNRows();

    _tl.resize(getTriangleSize());
    Eigen::MatrixXd mymat = _factor.matrixL();
    for (int ip = 0; ip < neq; ip++)
      for (int jp = 0; jp <= ip; jp++)
        _TL(ip,jp) = mymat(ip,jp);
  }
  else
  {
    int neq = getNRows();
    _tl.resize(getTriangleSize());

    for (int ip = 0; ip < neq; ip++)
      for (int jp = 0; jp <= ip; jp++)
        _TL(ip,jp)= getValue(ip,jp);

    for (int ip = 0; ip < neq; ip++)
    {
      double prod = _TL(ip, ip);
      for (int kp = 0; kp < ip; kp++)
        prod -= _TL(ip,kp)* _TL(ip,kp);
      if (prod < 0.) return 1;
      _TL(ip,ip)= sqrt(prod);

      for (int jp = ip + 1; jp < neq; jp++)
      {
        prod = _TL(jp, ip);
        for (int kp = 0; kp < ip; kp++)
          prod -= _TL(ip,kp)* _TL(jp,kp);
        if (_TL(ip,ip)<= 0.) return 1;
        _TL(jp,ip)= prod / _TL(ip,ip);
      }
    }
  }
  _flagCholeskyDecompose = true;
  return 0;
}

double MatrixSquareSymmetric::computeCholeskyLogDeterminant() const
{
  if (! isFlagEigen())
  {
    messerr("LogDet from Cholesky is only coded for Eigen library");
    return TEST;
  }
  if (! _checkCholeskyAlreadyPerformed(1)) return TEST;

  auto diag = _factor.matrixLLT().diagonal();
  double det = 0.;
  for (int i = 0; i < _factor.rows(); i++)
    det += log(diag[i]);
  return det;
}

bool MatrixSquareSymmetric::_checkCholeskyAlreadyPerformed(int status) const
{
  if (status == 1 && ! _flagCholeskyDecompose)
  {
    messerr("This operation requires a previous call to choleskyDecompose()");
    return false;
  }
  if (status == 2 && ! _flagCholeskyInverse)
  {
    messerr("This operation requires a previous call to choleskyInvert()");
    return false;
  }
  return true;
}

VectorDouble MatrixSquareSymmetric::getCholeskyTL() const
{
  if (! _checkCholeskyAlreadyPerformed(1)) return VectorDouble();
  return _tl;
}

VectorDouble MatrixSquareSymmetric::getCholeskyXL() const
{
  if (! _checkCholeskyAlreadyPerformed(2)) return VectorDouble();
  return _xl;
}

/*****************************************************************************/
/*!
 **  Invert the Cholesky matrix
 **
 *****************************************************************************/
int MatrixSquareSymmetric::invertCholesky()
{
  if (! _checkCholeskyAlreadyPerformed(1)) return 1;

  int neq = getNRows();
  _xl.resize(getTriangleSize());
  _flagCholeskyInverse = false;

  for (int i = 0; i < neq; i++)
  {
    for (int j = 0; j < i; j++)
    {
      double sum = 0.;
      for (int l = j; l < i; l++)
        sum += _TL(i,l) * _XL(l,j);
      _XL(i,j)= - sum / _TL(i,i);
    }
    _XL(i,i) = 1. / _TL(i,i);
  }

  _flagCholeskyInverse = true;
  return 0;
}

int MatrixSquareSymmetric::solveCholeskyMat(const MatrixRectangular& b, MatrixRectangular& x)
{
  if (! isFlagEigen())
  {
    messerr("solve from Cholesky is only coded for Eigen library");
    return 1;
  }
  if (! _checkCholeskyAlreadyPerformed(1)) return 1;

  int nrows = b.getNRows();
  int ncols = b.getNCols();
  x.resize(nrows, ncols);

  VectorDouble xcol(nrows);
  for (int icol = 0; icol < ncols; icol++)
  {
    VectorDouble bcol = b.getColumn(icol);
    solveCholesky(bcol, xcol);
    x.setColumn(icol, xcol);
  }
  return 0;
}

int MatrixSquareSymmetric::solveCholesky(const VectorDouble& b, VectorDouble& x)
{
  if (! isFlagEigen())
  {
    messerr("solve from Cholesky is only coded for Eigen library");
    return 1;
  }
  if (! _checkCholeskyAlreadyPerformed(1)) return 1;

  int size = (int) b.size();
  Eigen::Map<const Eigen::VectorXd> bm(b.data(), size);
  Eigen::VectorXd xm = _factor.solve(bm);

  x.resize(size);
  Eigen::Map<Eigen::VectorXd>(&x[0], size) = xm;

  return 0;
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
 ** \param[in]  tl   Triangular matrix defined by column
 ** \param[in]  a    matrix (dimension neq * nrhs)
 **
 *****************************************************************************/
MatrixRectangular MatrixSquareSymmetric::productCholeskyInPlace(int mode,
                                                                int neq,
                                                                int nrhs,
                                                                const VectorDouble &tl,
                                                                const MatrixRectangular &a)
{
  MatrixRectangular x(neq, nrhs);

  double val = 0.;
  if (mode == 0)
  {
    for (int irhs = 0; irhs < nrhs; irhs++)
      for (int i = 0; i < neq; i++)
      {
        val = 0.;
        for (int j = i; j < neq; j++)
          val += TL(j,i) * a.getValue(j,irhs);
        x.setValue(i,irhs,val);
      }
    }
    else if (mode == 1)
    {
      for (int irhs=0; irhs<nrhs; irhs++)
        for (int i=0; i<neq; i++)
        {
          val = 0.;
          for (int j=0; j<=i; j++)
            val += TL(i,j) * a.getValue(j,irhs);
          x.setValue(i,irhs,val);
        }
    }
    else if (mode == 2)
    {
      for (int irhs=0; irhs<nrhs; irhs++)
        for (int i=0; i<neq; i++)
        {
          val = 0.;
          for (int j=0; j<=i; j++)
            val += a.getValue(irhs,j) * TL(i,j);
          x.setValue(irhs,i,val);
        }
    }
    else if (mode == 3)
    {
      for (int irhs=0; irhs<nrhs; irhs++)
        for (int i=0; i<neq; i++)
        {
          val = 0.;
          for (int j=i; j<neq; j++)
            val += a.getValue(irhs,j) * TL(j,i);
          x.setValue(irhs,i,val);
        }
    }
    else if (mode == 4)
    {
      for (int irhs=0; irhs<nrhs; irhs++)
        for (int i=0; i<neq; i++)
        {
          val = 0.;
          for (int j=0; j<=i; j++)
            val += a.getValue(irhs,j) * TL(i,j);
          x.setValue(irhs,i,val);
        }
    }
    else if (mode == 5)
    {
      for (int irhs=0; irhs<nrhs; irhs++)
        for (int i=0; i<neq; i++)
        {
          val = 0.;
          for (int j=i; j<neq; j++)
            val += a.getValue(irhs,j) * TL(j,i);
          x.setValue(irhs,i,val);
        }
    }
  return x;
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
 *****************************************************************************/
MatrixSquareSymmetric MatrixSquareSymmetric::normCholeskyInPlace(int mode,
                                                                 int neq,
                                                                 const VectorDouble &tl,
                                                                 const MatrixSquareSymmetric &a)
{
  MatrixSquareSymmetric b(neq);
  double vala;

  for (int i = 0; i < neq; i++)
    for (int j = 0; j < neq; j++)
    {
      double val = 0.;
      if (mode == 0)
      {
        for (int l = 0; l <= j; l++)
          for (int k = 0; k <= i; k++)
          {
            if (! a.empty())
              vala = a.getValue(k, l);
            else
              vala = (k == l);
            val += TL(i,k) * vala * TL(j,l);
          }
        }
        else
        {
          for (int l=j; l<neq; l++)
            for (int k=i; k<neq; k++)
            {
              if (! a.empty())
                vala = a.getValue(k,l);
              else
                vala = (k == l);
              val += TL(k,i) * vala * TL(l,j);
            }
        }
      b.setValue(i,j,val);
    }
  return b;
}

/*****************************************************************************/
/*!
 **  Solve a linear system: H %*% g = x
 **
 ** \return  Error return code
 **
 ** \param[in]  gmat    right-hand side vector (Dimension: neq)
 **
 ** \param[out] xmat    solution vector (Dimension: neq)
 **
 ** \remark In output, 'this' contains the inverse matrix
 **
 *****************************************************************************/
int MatrixSquareSymmetric::_matrix_qo(const VectorDouble& gmat, VectorDouble& xmat)
{
  if (computeGeneralizedInverse(*this)) return 1;
  prodMatVecInPlace(gmat, xmat);
  return 0;
}

/*****************************************************************************/
/*!
 **  Minimize 1/2 t(x) %*% H %*% x + t(g) %*% x under the constraints
 **  t(A) %*% x = b
 **
 ** \return  Error return code
 **
 ** \param[in]  flag_invert Tells if the inverse has already been calculated
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
int MatrixSquareSymmetric::_matrix_qoc(bool flag_invert,
                                       const VectorDouble& gmat,
                                       int na,
                                       const MatrixRectangular& amat,
                                       const VectorDouble& bmat,
                                       VectorDouble& xmat,
                                       VectorDouble& lambda)
{
  double value;

  /* Initializations */

  int error = 1;
  int neq = getNRows();

  /* Core allocation */

  double* ha   = (double*) mem_alloc(sizeof(double) * neq * na, 1);
  double* evec = (double*) mem_alloc(sizeof(double) * na, 1);
  MatrixSquareSymmetric temp(na);

  /* Preliminary solution of the linear system with no constraint */

  if (!flag_invert)
  {
    if (_matrix_qo(gmat, xmat)) goto label_end;
  }

  /* Product HA = H %*% A */

  for (int i = 0; i < neq; i++)
    for (int j = 0; j < na; j++)
    {
      value = 0.;
      for (int k = 0; k < neq; k++)
        value += getValue(i,k) * amat.getValue(k,j);
      HA(i,j) = value;
    }

    /* Product temp = t(A) %*% H %*% A */

  for (int i = 0; i < na; i++)
    for (int j = 0; j < na; j++)
    {
      value = 0.;
      for (int k = 0; k < neq; k++)
        value += amat.getValue(k,i) * HA(k,j);
      temp.setValue(i,j,value);
    }

    /* Generalized inverse of temp */

  if (temp.computeGeneralizedInverse(temp)) goto label_end;

  /* Evaluate evec = t(A) %*% x - b */

  for (int i = 0; i < na; i++)
  {
    value = 0.;
    for (int j = 0; j < neq; j++)
      value += amat.getValue(j,i) * xmat[j];
    evec[i] = value - bmat[i];
  }

  /* Evaluate lambda = temp %*% evec */

  for (int i = 0; i < na; i++)
  {
    value = 0.;
    for (int j = 0; j < na; j++)
      value += temp.getValue(i,j) * evec[j];
    lambda[i] = value;
  }

  /* Evaluate x = x - H %*% A %*% lambda */

  for (int i = 0; i < neq; i++)
  {
    value = 0.;
    for (int j = 0; j < na; j++)
      value += HA(i,j) * lambda[j];
    xmat[i] -= value;
  }

  /* Set the error return code */

  error = 0;

  label_end:
  ha   = (double*) mem_free((char* ) ha);
  evec = (double*) mem_free((char* ) evec);
  return (error);
}

/*****************************************************************************/
/*!
 **  Minimize 1/2 t(x) %*% H %*% x + t(g) %*% x under the constraints
 **  t(Ae) %*% x = be and
 **  t(Ai) %*% x = bi
 **
 ** \return  Error return code
 **
 ** \param[in]     gmat   right-hand side vector (Dimension: neq)
 ** \param[in]     aemat  Matrix rectangular for equalities (Dimension: neq * nae)
 ** \param[in]     bemat  right-hand side for equalities (Dimension: nae)
 ** \param[in]     aimat  Matrix rectangular for inequalities (Dimension: neq * nai)
 ** \param[in]     bimat  right-hand side for inequalities (Dimension: nai)
 **
 ** \param[in,out] xmat solution of the linear system with constraints (neq)
 **
 ** REMARKS:    The initial xmat has to be satisfied by all the constraints.
 **
 *****************************************************************************/
int MatrixSquareSymmetric::minimizeWithConstraintsInPlace(const VectorDouble& gmat,
                                                          const MatrixRectangular& aemat,
                                                          const VectorDouble& bemat,
                                                          const MatrixRectangular& aimat,
                                                          const VectorDouble& bimat,
                                                          VectorDouble& xmat)
{
  int ncur, first, lec;
  double omega, omin, value;

  /* Initializations */

  int neq = getNRows();
  int nae = aemat.getNCols();
  int nai = aimat.getNCols();
  int namax = nae + nai;

  /* Case when there is no equality nor inequality constraints */

  if (namax <= 0)
  {
    return _matrix_qo(gmat, xmat);
  }

  /* Core allocation */

  VectorInt emptyInt;
  VectorDouble emptyDouble;
  VectorInt active(nai);
  VectorDouble xcand(neq);
  VectorDouble lambda(namax);
  VectorDouble vmat(namax);
  VectorDouble beimat(namax);
  MatrixRectangular aeimat(neq, namax);

  /* We first perform the optimization with equality constraints only */

  if (_matrix_qoc(false, gmat, nae, aemat, bemat, xcand, lambda)) return 1;
  if (nai <= 0)
  {
    for (int i = 0; i < neq; i++)
      xmat[i] = xcand[i];
    return 0;
  }

  /* Evaluate the array active */

  if (_constraintsError(VectorInt(), aimat, bimat, xcand, emptyDouble, active) == 0)
  {
    for (int i = 0; i < neq; i++)
      xmat[i] = xcand[i];
    return 0;
  }

  /* Implicit loop */

  bool sortie = false;
  while (!sortie)
  {

    /* Construct the inequality matrices reduced to the active constraints */

    ncur = _constraintsConcatenateMat(nae, nai, neq, active, aemat, aimat, aeimat);
    ncur = _constraintsConcatenateVD(nae, nai, active, bemat, bimat, beimat);
    if (_matrix_qoc(true, gmat, ncur, aeimat, beimat, xcand, lambda)) return 1;

    if (_constraintsError(active, aimat, bimat, xcand, vmat, emptyInt) == 0)
    {
      for (int i = 0; i < neq; i++)
        xmat[i] = xcand[i];

      /* Look for the constraint that should not be used */

      first = -1;
      lec = nae;
      for (int i = 0; i < nai; i++)
      {
        if (active[i] == 0) continue;
        active[i] = lambda[lec] >= 0;
        if (active[i]) first = i;
        lec++;
      }

      if (_constraintsCount(nai, active) == 0)
      {
        /* If no constraint has been used, end of the implicit loop */
        sortie = true;
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
      for (int i = 0; i < nai; i++)
      {
        if (active[i]) continue;
        value = 0.;
        for (int j = 0; j < neq; j++)
          value += aimat.getValue(j,i)* (xcand[j] - xmat[j]);
        omega = vmat[i] / value;
        if (omega > omin) continue;
        first = i;
        omin = omega;
      }

      for (int i = 0; i < neq; i++)
        xmat[i] += omin * (xcand[i] - xmat[i]);
      active[first] = 1;
    }
  }
  return 0;
}

/*****************************************************************************/
/*!
 **  Calculate how constraints are fulfilled
 **
 **  \return Count of the constraints not fulfilled
 **
 ** \param[in]  active   Array of active/non active inequalities (optional)
 ** \param[in]  aimat    Inequality material (Dimension: neq * nai)
 ** \param[in]  bimat    right-hand side for inequalities (Dimension: nai)
 ** \param[out] xmat     solution of the linear system with no constraint (neq)
 **
 ** \param[out] vmat     matrix of errors (if not NULL)
 ** \param[out] flag     array specifying if constraint is active (if not NULL)
 **
 *****************************************************************************/
int MatrixSquareSymmetric::_constraintsError(const VectorInt& active,
                                             const MatrixRectangular& aimat,
                                             const VectorDouble& bimat,
                                             const VectorDouble& xmat,
                                             VectorDouble& vmat,
                                             VectorInt& flag)
{
  double eps = EPSILON10;

  int neq = getNRows();
  int nai = aimat.getNCols();
  int number = 0;
  int ecr = 0;
  for (int i = 0; i < nai; i++)
  {
    if (! active.empty() && active[i]) continue;

    /* Calculate: T(a) %*% x */

    double value = 0.;
    for (int j = 0; j < neq; j++)
      value += aimat.getValue(j,i) * xmat[j];

      /* Calculate: T(a) %*% x - b */

    double ecart = value - bimat[i];

    /* Store the results */

    if (! vmat.empty()) vmat[ecr] = ecart;
    bool flag_active = (ecart < -eps);
    if (! flag.empty()) flag[ecr] = flag_active;
    if (flag_active) number++;
    ecr++;
  }
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
 ** \param[in]  neq      First dimension of the array
 ** \param[in]  active   Array of active/non active inequalities
 ** \param[in]  tabemat  Equality material (Dimension: neq * nai)
 ** \param[in]  tabimat  Inequality material
 **
 ** \param[out] tabout   Output array
 **
 *****************************************************************************/
int MatrixSquareSymmetric::_constraintsConcatenateMat(int nae,
                                                      int nai,
                                                      int neq,
                                                      const VectorInt& active,
                                                      const MatrixRectangular &tabemat,
                                                      const MatrixRectangular &tabimat,
                                                      MatrixRectangular &tabout)
{
  /* Copy the equalities */

  int number = 0;
  for (int i = 0; i < nae; i++)
  {
    for (int j = 0; j < neq; j++)
    {
      tabout.setValue(j,number,tabemat.getValue(j,i));
    }
    number++;
  }

    /* Copy the active inequalities */

  for (int i = 0; i < nai; i++)
  {
    if (active[i] == 0) continue;
    for (int j = 0; j < neq; j++)
    {
      tabout.setValue(j,number,tabimat.getValue(j,i));
    }
    number++;
  }
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
 ** \param[in]  tabemat  Equality material (Dimension: neq * nai)
 ** \param[in]  tabimat  Inequality material
 **
 ** \param[out] tabout   Output array
 **
 *****************************************************************************/
int MatrixSquareSymmetric::_constraintsConcatenateVD(int nae,
                                                     int nai,
                                                     const VectorInt &active,
                                                     const VectorDouble &tabemat,
                                                     const VectorDouble &tabimat,
                                                     VectorDouble &tabout)
{
  /* Copy the equalities */

  int number = 0;
  for (int i = 0; i < nae; i++)
  {
    tabout[number] = tabemat[i];
    number++;
  }

    /* Copy the active inequalities */

  for (int i = 0; i < nai; i++)
  {
    if (active[i] == 0) continue;
    tabout[number] = tabimat[i];
    number++;
  }
  return (number);
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
int MatrixSquareSymmetric::_constraintsCount(int nai, VectorInt& active)
{
  int number = 0;
  for (int i = 0; i < nai; i++)
    if (active[i]) number++;
  return (number);
}

/****************************************************************************/
/*!
 **  Calculate the generalized inverse of the input square symmetric matrix
 **
 ** \return  Error returned code
 **

 ** \param[out] tabout    Inverted matrix (suqrae symmetric)
 ** \param[out] maxicond  Maximum value for the Condition Index (MAX(ABS(eigval)))
 ** \param[in]  eps       Tolerance
 **
 ** \remark The input and output matrices can match
 **
 *****************************************************************************/
int MatrixSquareSymmetric::computeGeneralizedInverse(MatrixSquareSymmetric &tabout,
                                                     double maxicond,
                                                     double eps)
{
  if (! isSameSize(tabout))
  {
    messerr("The argument 'tabout' must have same dimensions as input matrix");
    return 1;
  }

  // Calculate the Eigen vectors
  if (computeEigen()) return 1;
  VectorDouble eigval = getEigenValues();
  const MatrixSquareGeneral *eigvec = getEigenVectors();

  // Compute the conditioning

  double valcond = VH::maximum(eigval, true);
  if (valcond > maxicond)
    return 1;

  /* Calculate the generalized inverse */

  int neq = getNRows();
  for (int i = 0; i < neq; i++)
    for (int j = 0; j < neq; j++)
    {
      double value = 0.;
      for (int k = 0; k < neq; k++)
      {
        if (ABS(eigval[k]) > valcond * eps)
          value += eigvec->getValue(i,k) * eigvec->getValue(j,k) / eigval[k];
      }
      tabout.setValue(i, j, value);
    }
  return 0;
}
