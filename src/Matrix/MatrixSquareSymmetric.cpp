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

MatrixSquareSymmetric::MatrixSquareSymmetric(int nrow, int opt_eigen)
    : AMatrixSquare(nrow, opt_eigen),
      _squareSymMatrix()
{
  _allocate();
}

MatrixSquareSymmetric::MatrixSquareSymmetric(const MatrixSquareSymmetric &r) 
  : AMatrixSquare(r),
   _squareSymMatrix()
{
  _recopyLocal(r);
}

MatrixSquareSymmetric::MatrixSquareSymmetric(const AMatrix &m)
    : AMatrixSquare(m),
      _squareSymMatrix()
{
  if (!m.isSymmetric())
  {
    messerr("The input matrix should be Symmetric");
    _clear();
    return;
  }
  const MatrixSquareSymmetric* matrixLoc = dynamic_cast<const MatrixSquareSymmetric*>(&m);
  if (matrixLoc != nullptr)
    _recopyLocal(*matrixLoc);
  else
  {
    _allocate();
    AMatrix::copyElements(m);
  }
}

MatrixSquareSymmetric& MatrixSquareSymmetric::operator= (const MatrixSquareSymmetric &r)
{
  if (this != &r)
  {
    AMatrixSquare::operator=(r);
    _recopyLocal(r);
  }
  return *this;
}

MatrixSquareSymmetric::~MatrixSquareSymmetric()
{
  _deallocate();
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

double MatrixSquareSymmetric::_getValue(int irow, int icol) const
{
  if (isFlagEigen())
    return AMatrixDense::_getValue(irow, icol);
  else
    return _getValueLocal(irow, icol);
}

double MatrixSquareSymmetric::_getValueByRank(int irank) const
{
  if (isFlagEigen())
    return AMatrixDense::_getValueByRank(irank);
  else
    return _getValueLocal(irank);
}

double& MatrixSquareSymmetric::_getValueRef(int irow, int icol)
{
  if (isFlagEigen())
    return AMatrixDense::_getValueRef(irow, icol);
  else
    return _getValueRef(irow, icol);
}

void MatrixSquareSymmetric::_setValue(int irow, int icol, double value)
{
  if (isFlagEigen())
  {
    // Do not forget to make a symmetrical call (when stored in an Eigen format)
    AMatrixDense::_setValue(irow, icol, value);
    if (irow != icol) AMatrixDense::_setValue(icol, irow, value);
  }
  else
    _setValueLocal(irow, icol, value);
}

void MatrixSquareSymmetric::_updValue(int irow, int icol, const EOperator& oper, double value)
{
  if (isFlagEigen())
  {
    // Do not forget to make a symmetrical call (when stored in an Eigen format)
    AMatrixDense::_updValue(irow, icol, oper, value);
    if (irow != icol) AMatrixDense::_updValue(icol, irow, oper, value);
  }
  else
    _updValueLocal(irow, icol, oper, value);
}

void MatrixSquareSymmetric::_setValueByRank(int irank, double value)
{
  if (isFlagEigen())
    AMatrixDense::_setValueByRank(irank, value);
  else
    _setValueLocal(irank, value);
}

void MatrixSquareSymmetric::_prodMatVecInPlacePtr(const double *x, double *y, bool transpose) const
{
  if (isFlagEigen())
    AMatrixDense::_prodMatVecInPlacePtr(x, y, transpose);
  else
    _prodMatVecInPlacePtrLocal(x, y, transpose);
}

void MatrixSquareSymmetric::_prodVecMatInPlacePtr(const double *x, double *y, bool transpose) const
{
  if (isFlagEigen())
    AMatrixDense::_prodVecMatInPlacePtr(x, y, transpose);
  else
    _prodVecMatInPlacePtrLocal(x, y, transpose);
}

/**
 * \warning : values is provided as a square complete matrix
 */
void MatrixSquareSymmetric::_setValues(const double* values, bool byCol)
{
  if (isFlagEigen())
    AMatrixDense::_setValues(values, byCol);
  else
    _setValuesLocal(values, byCol);
}

int MatrixSquareSymmetric::_invert()
{
  if (isFlagEigen())
    return AMatrixDense::_invert();
  else
    return _invertLocal();
}

void MatrixSquareSymmetric::_deallocate()
{
  if (isFlagEigen())
    AMatrixDense::_deallocate();
  else
  {
    delete _eigenVectors;
  }
}

void MatrixSquareSymmetric::_allocate()
{
  if (isFlagEigen())
    AMatrixDense::_allocate();
  else
    _allocateLocal();
}

int MatrixSquareSymmetric::_getIndexToRank(int irow, int icol) const
{
  if (isFlagEigen())
    return AMatrixDense::_getIndexToRank(irow, icol);
  else
    return _getIndexToRankLocal(irow, icol);
}

int MatrixSquareSymmetric::_getMatrixPhysicalSize() const
{
  if (isFlagEigen())
    return AMatrixDense::_getMatrixPhysicalSize();
  else
    return _getMatrixPhysicalSizeLocal();
}

int MatrixSquareSymmetric::_solve(const VectorDouble& b, VectorDouble& x) const
{
  if (isFlagEigen())
    return AMatrixDense::_solve(b, x);
  else
    return _solveLocal(b, x);
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

MatrixSquareSymmetric* MatrixSquareSymmetric::createReduce(const VectorInt &validRows) const
{
  // Order and shrink the input vectors
  VectorInt localValidRows = VH::filter(validRows, 0, getNRows());
  int newNRows = (int) localValidRows.size();
  if (newNRows <= 0)
  {
    messerr("The new Matrix has no Row left");
    return nullptr;
  }

  MatrixSquareSymmetric* res = new MatrixSquareSymmetric(newNRows);
  res->copyReduce(this, localValidRows, localValidRows);

  return res;
}

int MatrixSquareSymmetric::computeEigen(bool optionPositive)
{
  if (isFlagEigen())
    return AMatrixDense::_computeEigen(optionPositive);
  else
    return _computeEigenLocal(optionPositive);
}

int MatrixSquareSymmetric::computeGeneralizedEigen(const MatrixSquareSymmetric& b, bool optionPositive)
{
  if (isFlagEigen())
    return AMatrixDense::_computeGeneralizedEigen(b, optionPositive);
  else
    return _computeGeneralizedEigenLocal(b, optionPositive);
}

/// =============================================================================
/// The subsequent methods rely on the specific local storage ('squareSymMatrix')
/// =============================================================================

void MatrixSquareSymmetric::_recopyLocal(const MatrixSquareSymmetric& r)
{
  _squareSymMatrix = r._squareSymMatrix;
  _flagEigenDecompose = r._flagEigenDecompose;
  _eigenValues = r._eigenValues;
  if (_eigenVectors != nullptr) _eigenVectors = r._eigenVectors->clone();
}

double MatrixSquareSymmetric::_getValueLocal(int irow, int icol) const
{
  if (! _isIndexValid(irow,icol)) return TEST;
  int rank = _getIndexToRank(irow,icol);
  return _squareSymMatrix[rank];
}

double MatrixSquareSymmetric::_getValueLocal(int irank) const
{
  return _squareSymMatrix[irank];
}

double& MatrixSquareSymmetric::_getValueRefLocal(int irow, int icol)
{
  int rank = _getIndexToRank(irow, icol);
  return _squareSymMatrix[rank];
}

void MatrixSquareSymmetric::_setValueLocal(int irow, int icol, double value)
{
  if (! _isIndexValid(irow, icol)) return;
  int irank = _getIndexToRank(irow, icol);
  _squareSymMatrix[irank] = value;
}

void MatrixSquareSymmetric::_updValueLocal(int irow, int icol, const EOperator& oper, double value)
{
  if (! _isIndexValid(irow, icol)) return;
  int irank = _getIndexToRank(irow, icol);
  _squareSymMatrix[irank] = modifyOperator(oper, _squareSymMatrix[irank], value);
}

void MatrixSquareSymmetric::_setValueLocal(int irank, double value)
{
  if (! _isRankValid(irank)) return;
  _squareSymMatrix[irank] = value;
}

void MatrixSquareSymmetric::_prodMatVecInPlacePtrLocal(const double *x, double *y, bool transpose) const
{
  _matrix_triangular_product(getNRows(),2,_squareSymMatrix.data(),x,y);
}

void MatrixSquareSymmetric::_prodVecMatInPlacePtrLocal(const double *x, double *y, bool transpose) const
{
  _matrix_triangular_product(getNRows(),2,_squareSymMatrix.data(),x,y);
}

/**
 * \warning : values is provided as a square complete matrix
 */
void MatrixSquareSymmetric::_setValuesLocal(const double *values, bool byCol)
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

  int ecr = 0;
  if (byCol)
  {
    for (int icol = 0; icol < getNCols(); icol++)
      for (int irow = 0; irow < getNRows(); irow++)
        setValue(irow, icol, values[ecr++]);
  }
  else
  {
    for (int irow = 0; irow < getNRows(); irow++)
      for (int icol = 0; icol < getNCols(); icol++)
        setValue(irow, icol, values[ecr++]);
  }
}

int MatrixSquareSymmetric::_invertLocal()
{
  return matrix_invert_triangle(getNRows(),_squareSymMatrix.data(), -1);
}

void MatrixSquareSymmetric::_allocateLocal()
{
  _squareSymMatrix.resize(_getMatrixPhysicalSize());
}

int MatrixSquareSymmetric::_getIndexToRankLocal(int irow, int icol) const
{
  int rank;
  int n = getNRows();
  if (irow >= icol)
    rank = icol * n + irow - icol * (icol + 1) / 2;
  else
    rank = irow * n + icol - irow * (irow + 1) / 2;
  return rank;
}

int MatrixSquareSymmetric::_getMatrixPhysicalSizeLocal() const
{
  int n = getNRows();
  int size = n * (n + 1) / 2;
  return size;
}

int MatrixSquareSymmetric::_solveLocal(const VectorDouble& b, VectorDouble& x) const
{
  int pivot;
  int size = (int) b.size();
  VectorDouble alocal = _squareSymMatrix;
  VectorDouble blocal = b;
  pivot = _matrix_solve(alocal, blocal, x, size, 1);
  return (pivot != 0);
}

int MatrixSquareSymmetric::_computeEigenLocal(bool optionPositive)
{
  int nrows = getNRows();
  _flagEigenDecompose = true;
  _eigenValues = VectorDouble(nrows, 0.);
  VectorDouble eigenVectors = VectorDouble(nrows * nrows, 0);

  int err = matrix_eigen(this->getValues().data(), nrows, _eigenValues.data(), eigenVectors.data());
  if (err == 0)
  {
    _eigenVectors = MatrixSquareGeneral::createFromVD(eigenVectors, nrows, false, 0, false);

    if (optionPositive) _eigenVectors->makePositiveColumn();
  }

  return err;
}

int MatrixSquareSymmetric::_computeGeneralizedEigenLocal(const MatrixSquareSymmetric& b, bool optionPositive)
{
  int nrows = getNRows();
  _flagEigenDecompose = true;
  _eigenValues = VectorDouble(nrows, 0.);
  VectorDouble eigenVectors = VectorDouble(nrows * nrows, 0);

  int err = _matrix_geigen(this->getValues().data(),b.getValues().data(), nrows, _eigenValues.data(), eigenVectors.data());
  if (err == 0)
  {
    std::reverse(_eigenValues.begin(), _eigenValues.end());
    _eigenVectors = MatrixSquareGeneral::createFromVD(eigenVectors, nrows, false, 0, true);

    if (optionPositive) _eigenVectors->makePositiveColumn();
  }
  return err;
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
        ratio -= AT(k,j)* XS(j,l);
      XS(k,l)= ratio / AT(k,k);
    }
  }
  return (0);
}
