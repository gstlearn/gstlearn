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

MatrixSquareSymmetric::MatrixSquareSymmetric(int nrow)
  : AMatrixSquare(nrow),
    _flagCholeskyDecompose(false),
    _flagCholeskyInverse(false),
    _tl(),
    _xl(),
    _factor()
{
}

MatrixSquareSymmetric::MatrixSquareSymmetric(const MatrixSquareSymmetric &m) 
  : AMatrixSquare(m),
   _flagCholeskyDecompose(m._flagCholeskyDecompose),
   _flagCholeskyInverse(m._flagCholeskyInverse),
   _tl(),
   _xl(),
   _factor()
{
  _recopy(m);
}

MatrixSquareSymmetric::MatrixSquareSymmetric(const AMatrix &m)
  : AMatrixSquare(m),
    _flagCholeskyDecompose(false),
    _flagCholeskyInverse(false),
    _tl(),
    _xl(),
    _factor()
{
  if (!m.isSymmetric())
  {
    messerr("The input matrix should be Symmetric");
    _clear();
    return;
  }
  copyElements(m);
}

MatrixSquareSymmetric& MatrixSquareSymmetric::operator= (const MatrixSquareSymmetric &m)
{
  if (this != &m)
  {
    AMatrixSquare::operator=(m);
    _recopy(m);
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
 * @return The returned square symmetric matrix
 *
 * @remark: the matrix is transposed implicitly while reading
 */
MatrixSquareSymmetric* MatrixSquareSymmetric::createFromVVD(const VectorVectorDouble& X)
{
  int nrow = (int) X.size();
  int ncol = (int) X[0].size();
  if (nrow != ncol)
  {
    messerr("The matrix does not seem to be square");
    return nullptr;
  }
  MatrixSquareSymmetric* mat = new MatrixSquareSymmetric(nrow);
  mat->_fillFromVVD(X);
  return mat;
}

MatrixSquareSymmetric* MatrixSquareSymmetric::createFromVD(const VectorDouble &X,
                                                           int nrow)
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

  MatrixSquareSymmetric *mat = new MatrixSquareSymmetric(nrow);

  int lec = 0;
  for (int irow = 0; irow < nrow; irow++)
    for (int icol = 0; icol < ncol; icol++)
      mat->setValue(irow, icol, X[lec++]);
  return mat;
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
  return AMatrixDense::_invert();
}

bool MatrixSquareSymmetric::_isPhysicallyPresent(int irow, int icol) const
{
  return (icol <= irow);
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

int MatrixSquareSymmetric::computeEigen(bool optionPositive)
{
  return AMatrixDense::_computeEigen(optionPositive);
}

int MatrixSquareSymmetric::computeGeneralizedEigen(const MatrixSquareSymmetric& b, bool optionPositive)
{
  return AMatrixDense::_computeGeneralizedEigen(b, optionPositive);
}

int MatrixSquareSymmetric::_terminateEigen(const VectorDouble &eigenValues,
                                           const VectorDouble &eigenVectors,
                                           bool optionPositive,
                                           bool changeOrder)
{
  int nrows = getNRows();

  _eigenValues = eigenValues;

  delete _eigenVectors;

  if (changeOrder)
    std::reverse(_eigenValues.begin(), _eigenValues.end());

  _eigenVectors = MatrixSquareGeneral::createFromVD(eigenVectors, nrows, false, changeOrder);

  if (optionPositive) _eigenVectors->makePositiveColumn();

  _flagEigenDecompose = true;

  return 0;
}

void MatrixSquareSymmetric::_recopy(const MatrixSquareSymmetric& r)
{
  _tl = r._tl;
  _xl = r._xl;
  _flagCholeskyDecompose = r._flagCholeskyDecompose;
  _flagCholeskyInverse   = r._flagCholeskyInverse;
  _flagEigenDecompose    = r._flagEigenDecompose;
  _factor                = r._factor;
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

  if (computeEigen() != 0) messageAbort("matrix_eigen");

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
 **  Create the Symmetric matrix as the product of 'tl' (lower triangle) by its transpose
 **
 ** \param[in]  neq    Number of rows or columns in the system
 ** \param[in]  tl     Lower triangular matrix defined by column (Dimension; neq*(neq+1)/2)
 **
 *****************************************************************************/
MatrixSquareSymmetric* MatrixSquareSymmetric::createFromTLTU(int neq,
                                                             const VectorDouble &tl)
{
  MatrixSquareSymmetric *mat = new MatrixSquareSymmetric(neq);

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
 ** \param[in]  tl     Triangular matrix (any part)
 **
 *****************************************************************************/
MatrixSquareSymmetric* MatrixSquareSymmetric::createFromTriangle(int mode,
                                                                 int neq,
                                                                 const VectorDouble &tl)
{
  MatrixSquareSymmetric* mat = new MatrixSquareSymmetric(neq);

  mat->fill(0.);

  for (int i = 0; i < neq; i++)
    for (int j = 0; j < neq; j++)
    {
      if (mode == 0)
      {
        if (j <= i) mat->setValue(i, j, TL(i, j));
      }
      else
      {
        if (j >= i) mat->setValue(i, j, TL(j, i));
      }
    }
  return mat;
}

int MatrixSquareSymmetric::_getTriangleSize() const
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

  _factor = _eigenMatrix.llt();
  int neq = getNRows();

  _tl.resize(_getTriangleSize());
  Eigen::MatrixXd mymat = _factor.matrixL();
  for (int ip = 0; ip < neq; ip++)
    for (int jp = 0; jp <= ip; jp++)
      _TL(ip,jp) = mymat(ip,jp);

  _flagCholeskyDecompose = true;
  return 0;
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

double MatrixSquareSymmetric::getCholeskyTL(int i, int j) const
{
  if (!_checkCholeskyAlreadyPerformed(1)) return TEST;
  int neq = getNRows();
  return (i >= j) ? _TL(i,j) : 0.;
}

VectorDouble MatrixSquareSymmetric::getCholeskyXL() const
{
  if (! _checkCholeskyAlreadyPerformed(2)) return VectorDouble();
  return _xl;
}

double MatrixSquareSymmetric::getCholeskyXL(int i, int j) const
{
  if (!_checkCholeskyAlreadyPerformed(2)) return TEST;
  int neq = getNRows();
  return (i >= j) ? _XL(i, j) : 0.;
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
  _xl.resize(_getTriangleSize());
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
  if (computeGeneralizedInverse(*this) != 0) return 1;
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

  int neq = getNRows();

  /* Core allocation */

  VectorDouble ha(neq * na);
  VectorDouble evec(na);
  MatrixSquareSymmetric temp(na);

  /* Preliminary solution of the linear system with no constraint */

  if (!flag_invert)
  {
    if (_matrix_qo(gmat, xmat) != 0) return 1;
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

  if (temp.computeGeneralizedInverse(temp) != 0) return 1;

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

  return 0;
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

  if (_matrix_qoc(false, gmat, nae, aemat, bemat, xcand, lambda) != 0) return 1;
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
    if (_matrix_qoc(true, gmat, ncur, aeimat, beimat, xcand, lambda) != 0) return 1;

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
        active[i] = (int) (lambda[lec] >= 0);
        if (active[i] != 0) first = i;
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
        if (active[i] != 0) continue;
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
    if (! active.empty() && active[i] != 0) continue;

    /* Calculate: T(a) %*% x */

    double value = 0.;
    for (int j = 0; j < neq; j++)
      value += aimat.getValue(j,i) * xmat[j];

      /* Calculate: T(a) %*% x - b */

    double ecart = value - bimat[i];

    /* Store the results */

    if (! vmat.empty()) vmat[ecr] = ecart;
    bool flag_active = (ecart < -eps);
    if (! flag.empty()) flag[ecr] = (int) flag_active;
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
    if (active[i] != 0) number++;
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
  if (computeEigen() != 0) return 1;
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

/**
 * @brief Create an output Square Symmetric Matrix by selecting some rows (and
columns)
 *        of the Input matrix 'A'
 *
 * @param A        Input Square Symmetric Matrix
 * @param rowKeep  Set of Rows (same for columns) to be kept
 * @param flagInvert when True, transform 'rowKeep' into 'rowDrop' 
 * @return Pointer to the newly created Square Symmetric Matrix
 */
MatrixSquareSymmetric*
MatrixSquareSymmetric::sample(const MatrixSquareSymmetric* A,
                              const VectorInt& rowKeep,
                              bool flagInvert)
{
  int ntotal = A->getNRows();
  VectorInt rows = rowKeep;
  if (rows.empty()) rows = VH::sequence(ntotal);
  if (flagInvert) rows = VH::complement(VH::sequence(ntotal), rows);

  int nrows = (int)rows.size();
  if (nrows <= 0) return nullptr;

  for (int irow = 0; irow < nrows; irow++)
  {
    if (!checkArg("Selected Row index", rows[irow], ntotal)) return nullptr;
  }

  MatrixSquareSymmetric* mat = new MatrixSquareSymmetric(nrows);
  for (int irow = 0; irow < nrows; irow++)
    for (int icol = 0; icol <= irow; icol++)
      mat->setValue(irow, icol, A->getValue(rows[irow], rows[icol]));
  return mat;
}

MatrixSquareSymmetric* MatrixSquareSymmetric::createRandomDefinitePositive(int neq, int seed)
{
  MatrixSquareSymmetric local(neq);
  local.fillRandom(seed);
  MatrixSquareSymmetric* mat = new MatrixSquareSymmetric(neq);
  mat->prodMatMatInPlace(&local, &local, true);
  return mat;
}
