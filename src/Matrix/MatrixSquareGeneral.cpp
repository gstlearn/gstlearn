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

#include "Matrix/MatrixSquareGeneral.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/AException.hpp"

#define SQ(i,j,neq)   ((j) * neq + (i))
#define A(i,j)         a[SQ(i,j,neq)]
#define B(i,j)         b[SQ(i,j,neq)]
#define C(i,j)         c[SQ(i,j,neqm1)]

MatrixSquareGeneral::MatrixSquareGeneral(int nrow, int opt_eigen)
  : AMatrixSquare(nrow, opt_eigen)
  , _squareMatrix()
{
  _allocate();
}

MatrixSquareGeneral::MatrixSquareGeneral(const MatrixSquareGeneral &r) 
  : AMatrixSquare(r),
    _squareMatrix()
{
  _recopyLocal(r);
}

MatrixSquareGeneral::MatrixSquareGeneral(const AMatrix &m)
  : AMatrixSquare(m),
    _squareMatrix()
{
  const MatrixSquareGeneral* matrixLoc = dynamic_cast<const MatrixSquareGeneral*>(&m);
  if (matrixLoc != nullptr)
    _recopyLocal(*matrixLoc);
  else
  {
    _allocate();
    AMatrix::copyElements(m);
  }
}

MatrixSquareGeneral& MatrixSquareGeneral::operator= (const MatrixSquareGeneral &r)
{
  if (this != &r)
  {
    AMatrixSquare::operator=(r);
    _recopyLocal(r);
  }
  return *this;
}

MatrixSquareGeneral::~MatrixSquareGeneral()
{
  _deallocate();
}

/**
 * Converts a VectorVectorDouble into a Matrix
 * Note: the input argument is stored by row (if coming from [] specification)
 * @param X Input VectorVectorDouble argument
 * @param opt_eigen Option for use of Eigen Library
 * @return The returned square matrix
 *
 * @remark: the matrix is transposed implicitly while reading
 */
MatrixSquareGeneral* MatrixSquareGeneral::createFromVVD(const VectorVectorDouble& X, int opt_eigen)
{
  int nrow = (int) X.size();
  int ncol = (int) X[0].size();
  if (nrow != ncol)
  {
    messerr("The matrix does not seem to be square");
    return nullptr;
  }

  MatrixSquareGeneral* mat = new MatrixSquareGeneral(nrow, opt_eigen);
  mat->_fillFromVVD(X);
  return mat;
}

MatrixSquareGeneral* MatrixSquareGeneral::createFromVD(const VectorDouble &X,
                                                       int nrow,
                                                       bool byCol,
                                                       int opt_eigen,
                                                       bool invertColumnOrder)
{
  int ncol = nrow;
  if (nrow * ncol != (int) X.size())
  {
    messerr("Inconsistency between arguments 'nrow'(%d) and 'ncol'(%d)", nrow, ncol);
    messerr("and the dimension of the input Vector (%d)", (int) X.size());
  }
  MatrixSquareGeneral* mat = new MatrixSquareGeneral(nrow, opt_eigen);

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

double MatrixSquareGeneral::_getValue(int irow, int icol) const
{
  if (isFlagEigen())
    return AMatrixDense::_getValue(irow, icol);
  else
    return _getValueLocal(irow, icol);
}

double MatrixSquareGeneral::_getValueByRank(int irank) const
{
  if (isFlagEigen())
    return AMatrixDense::_getValueByRank(irank);
  else
    return _getValueLocal(irank);
}

double& MatrixSquareGeneral::_getValueRef(int irow, int icol)
{
  if (isFlagEigen())
    return AMatrixDense::_getValueRef(irow, icol);
  else
    return _getValueRefLocal(irow, icol);
}

void MatrixSquareGeneral::_setValue(int irow, int icol, double value)
{
  if (isFlagEigen())
    AMatrixDense::_setValue(irow, icol, value);
  else
    _setValueLocal(irow, icol, value);
}

void MatrixSquareGeneral::_updValue(int irow, int icol, const EOperator& oper, double value)
{
  if (isFlagEigen())
    AMatrixDense::_updValue(irow, icol, oper, value);
  else
    _updValueLocal(irow, icol, oper, value);
}

void MatrixSquareGeneral::_setValueByRank(int irank, double value)
{
  if (isFlagEigen())
    AMatrixDense::_setValueByRank(irank, value);
  else
    _setValueLocal(irank, value);
}

/**
 * Returns 'y' = 'this' %*% 'x'
 * @param x  Input Vector
 * @param y Output Vector
 * @param transpose True if the matrix 'this' must be transposed
 */
void MatrixSquareGeneral::_prodMatVecInPlacePtr(const double *x, double *y, bool transpose) const
{
  if (isFlagEigen())
    AMatrixDense::_prodMatVecInPlacePtr(x, y, transpose);
  else
    _prodMatVecInPlacePtrLocal(x, y, transpose);
}

/**
 * Returns 'y' = 'this' %*% 'x'
 * @param x  Input Vector
 * @param y Output Vector
 * @param transpose True if the matrix 'this' must be transposed
 */
void MatrixSquareGeneral::_prodVecMatInPlacePtr(const double *x, double *y, bool transpose) const
{
  if (isFlagEigen())
    AMatrixDense::_prodVecMatInPlacePtr(x, y, transpose);
  else
    _prodVecMatInPlacePtrLocal(x, y, transpose);
}

void MatrixSquareGeneral::_transposeInPlace()
{
  if (isFlagEigen())
    AMatrixDense::_transposeInPlace();
  else
    _transposeInPlaceLocal();
}

int MatrixSquareGeneral::_invert()
{
  if (isFlagEigen())
    return AMatrixDense::_invert();
  else
    return _invertLocal();
}

void MatrixSquareGeneral::_deallocate()
{
  if (isFlagEigen())
    AMatrixDense::_deallocate();
  else
  {
    // Specific code for this class should be placed here
  }
}

void MatrixSquareGeneral::_allocate()
{
  if (isFlagEigen())
    AMatrixDense::_allocate();
  else
    _allocateLocal();
}

int MatrixSquareGeneral::_getMatrixPhysicalSize() const
{
  return(getNRows() * getNCols());
}

int MatrixSquareGeneral::_solve(const VectorDouble& /*b*/, VectorDouble& /*x*/) const
{
  my_throw("Invert method is limited to Square Symmetrical Matrices");
  return 0;
}

/// ==========================================================================
/// The subsequent methods rely on the specific local storage ('squareMatrix')
/// ==========================================================================

void MatrixSquareGeneral::_allocateLocal()
{
  _squareMatrix.resize(_getMatrixPhysicalSize());
}

void MatrixSquareGeneral::_recopyLocal(const MatrixSquareGeneral &r)
{
  _squareMatrix = r._squareMatrix;
}

double MatrixSquareGeneral::_getValueLocal(int irow, int icol) const
{
  if (!_isIndexValid(irow, icol)) return TEST;
  int rank = _getIndexToRank(irow, icol);
  return _squareMatrix[rank];
}

double MatrixSquareGeneral::_getValueLocal(int irank) const
{
  return _squareMatrix[irank];
}

double& MatrixSquareGeneral::_getValueRefLocal(int irow, int icol)
{
  int rank = _getIndexToRank(irow,icol);
  return _squareMatrix[rank];
}

void MatrixSquareGeneral::_setValueLocal(int irow, int icol, double value)
{
  if (!_isIndexValid(irow, icol)) return;
  int rank = _getIndexToRank(irow, icol);
  _squareMatrix[rank] = value;
}

void MatrixSquareGeneral::_updValueLocal(int irow, int icol, const EOperator& oper, double value)
{
  if (!_isIndexValid(irow, icol)) return;
  int rank = _getIndexToRank(irow, icol);
  _squareMatrix[rank] = modifyOperator(oper, _squareMatrix[rank], value);
}

void MatrixSquareGeneral::_setValueLocal(int irank, double value)
{
  if (!_isRankValid(irank)) return;
  _squareMatrix[irank] = value;
}

void MatrixSquareGeneral::_prodMatVecInPlacePtrLocal(const double *x, double *y, bool transpose) const
{
  if (transpose)
    matrix_product_safe(1, getNRows(), getNCols(), x, _squareMatrix.data(), y);
  else
    matrix_product_safe(getNRows(), getNCols(), 1, _squareMatrix.data(), x, y);
}

void MatrixSquareGeneral::_prodVecMatInPlacePtrLocal(const double *x, double *y, bool transpose) const
{
  if (transpose)
    matrix_product_safe(getNRows(), getNCols(), 1, _squareMatrix.data(), x, y);
  else
    matrix_product_safe(1, getNRows(), getNCols(), x, _squareMatrix.data(), y);
}

void MatrixSquareGeneral::_transposeInPlaceLocal()
{
  int nrow = getNRows();
  int ncol = getNCols();
  VectorDouble old = _squareMatrix;
  matrix_transpose(nrow, ncol, _squareMatrix, old);
  _squareMatrix = old;
  _setNCols(nrow);
  _setNRows(ncol);
}

int MatrixSquareGeneral::_invertLocal()
{
  if (getNRows() <= 3)
    return _matrix_invreal(_squareMatrix, getNRows());
  else
  {
    int error = _invertLU();
    return error;
  }
}

MatrixSquareGeneral* prodNormMatMat(const AMatrixDense &a,
                                    const AMatrixDense &m,
                                    bool transpose)
{
  int nrow = (transpose) ? a.getNCols() : a.getNRows();
  MatrixSquareGeneral *mat = new MatrixSquareGeneral(nrow, a.isFlagEigen());
  mat->prodNormMatMatInPlace(a, m, transpose);
  return mat;
}

MatrixSquareGeneral* prodNormMat(const AMatrixDense &a, const VectorDouble& vec, bool transpose)
{
  int nsym = (transpose) ? a.getNCols() : a.getNRows();
  MatrixSquareGeneral *mat = new MatrixSquareGeneral(nsym, a.isFlagEigen());
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
int MatrixSquareGeneral::_matrix_invreal(VectorDouble& mat, int neq)
{
  /* Calculate the determinant */

  double det = matrix_determinant(neq, mat);
  if (isZero(det)) return 1;

  if (neq > 1)
  {

    /* Core allocation */

    VectorDouble cofac(neq * neq);

    /* Calculate the cofactor */

    if (_matrix_cofactor(neq, mat, cofac)) return 1;

    /* Transpose the cofactor to obtain the adjoint matrix */

    matrix_transpose(neq, neq, cofac, mat);
  }

  /* Final normation */

  for (int i = 0; i < neq * neq; i++)
    mat[i] /= det;

  return 0;
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
int MatrixSquareGeneral::_matrix_cofactor(int neq, VectorDouble& a, VectorDouble& b)
{
  // Process the case when the matrix A is of dimension 1

  int neqm1 = neq - 1;
  if (neqm1 <= 0)
  {
    B(0,0)= 1.;
    return 0;
  }
  VectorDouble c(neqm1 * neqm1);

  /* Processing */

  for (int j = 0; j < neq; j++)
  {
    for (int i = 0; i < neq; i++)
    {

      /* Form the adjoint a_ij */

      int i1 = 0;
      for (int ii = 0; ii < neq; ii++)
      {
        if (ii == i) continue;
        int j1 = 0;
        for (int jj = 0; jj < neq; jj++)
        {
          if (jj == j) continue;
          C(i1,j1) = A(ii,jj);
          j1++;
        }
        i1++;
      }

      /* Calculate the determinate */
      double det = matrix_determinant(neqm1, c);

      /* Fill in the elements of the cofactor */
      B(i,j) = pow(-1.0, i+j+2.0) * det;
    }
  }
  return (0);
}

