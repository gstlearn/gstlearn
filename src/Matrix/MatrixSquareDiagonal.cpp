/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "Matrix/MatrixSquareDiagonal.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Basic/AException.hpp"

MatrixSquareDiagonal::MatrixSquareDiagonal(int nrows)
  : MatrixSquareSymmetric(nrows)
  , _diagMatrix()
{
  _allocate();
}

MatrixSquareDiagonal::MatrixSquareDiagonal(const MatrixSquareDiagonal &r) 
  : MatrixSquareSymmetric(r)
{
  _recopy(r);
}

MatrixSquareDiagonal& MatrixSquareDiagonal::operator= (const MatrixSquareDiagonal &r)
{
  if (this != &r)
  {
    _deallocate();
    MatrixSquareSymmetric::operator=(r);
    _recopy(r);
  }
  return *this;
}

MatrixSquareDiagonal::~MatrixSquareDiagonal()
{
  _deallocate();
}

int MatrixSquareDiagonal::_getIndexToRank(int irow, int icol) const
{
  return icol;
}

double MatrixSquareDiagonal::_getValue(int irow, int icol) const
{
  if (irow == icol)
    return _diagMatrix[irow];
  else
    return 0.;
}

double MatrixSquareDiagonal::_getValue(int irank) const
{
  _isRankValid(irank);
  return _diagMatrix[irank];
}

void MatrixSquareDiagonal::_setValue(int irow, int icol, double value)
{
  if (irow == icol)
  {
    _diagMatrix[irow] = value;
  }
  else
  {
    if (ABS(value) > EPSILON10)
    {
      messerr("Attempt to assign a non-zero value to a non-diagonal term");
      messerr("- Element(%d,%d) = %lf", irow, icol, value);
      messerr("Operation is aborted");
      return;
    }
  }
}

void MatrixSquareDiagonal::_setValue(int irank, double value)
{
  _isRankValid(irank);
  _diagMatrix[irank] = value;
}

void MatrixSquareDiagonal::_prodVector(const double *inv, double *outv) const
{
  int nrow = getNRows();
  for (int irow=0; irow<nrow; irow++)
    outv[irow] = inv[irow] * _diagMatrix[irow];
}

void MatrixSquareDiagonal::transposeInPlace()
{
  // Nothing should be done
  return;
}

/**
 * Fill 'this' with the contents of 'values'.
 * Note that 'values' is dimensioned to 'nrows' by 'ncols'
 * @param values Input array
 * @param byCol True is the values are sorted by Column
 */
void MatrixSquareDiagonal::_setValues(const double* values, bool /*byCol*/)
{
  int ecr = 0;
  for (int icol = 0; icol < getNCols(); icol++)
    for (int irow = 0; irow < getNRows(); irow++, ecr++)
    {
      if (irow != icol)
      {
        if (ABS(values[ecr]) > EPSILON10)
        {
          messerr("Input 'values' does not correspond to a Diagonal matrix");
          messerr("- Element(%d,%d) = %lf", irow, icol, values[ecr]);
          messerr("Operation is aborted");
          return;
        }
      }
      else
      {
        setValue(irow, icol, values[ecr]);
      }
    }
}

double MatrixSquareDiagonal::determinant() const
{
  double deter = 1.;
  for (int irow = 0; irow < getNRows(); irow++)
    deter *= _diagMatrix[irow];
  return deter;
}

int MatrixSquareDiagonal::_invert()
{
  for (int irow=0; irow<getNRows(); irow++) 
  {
    if (_diagMatrix[irow] == 0) return 1;
    _diagMatrix[irow] = 1. / _diagMatrix[irow];
  }
  return 0;
}

double& MatrixSquareDiagonal::_getValueRef(int irow, int icol)
{
  if (! _isIndexValid(irow,icol))
    throw("Function aborted due to previous errors");
  return _diagMatrix[irow];
}

void MatrixSquareDiagonal::_deallocate()
{
}

void MatrixSquareDiagonal::_recopy(const MatrixSquareDiagonal &r)
{
  _diagMatrix = r._diagMatrix;
}

void MatrixSquareDiagonal::_allocate()
{
  _diagMatrix.resize(_getMatrixSize(),0.);
}

bool MatrixSquareDiagonal::_isIndexValid(int irow, int icol) const
{
  AMatrix::_isIndexValid(irow,icol);
  if (irow != icol)
  {
    messerr("Argument 'irow' and 'icol' should be equal for Diagonal Matrix");
    return false;
  }
  return true;
}

int MatrixSquareDiagonal::_getMatrixSize() const
{
  return (getNRows());
}

bool MatrixSquareDiagonal::isValid(int irow, int icol, bool printWhyNot) const
{
  if (! AMatrix::isValid(irow,icol,printWhyNot)) return false;
  if (irow != icol)
  {
    if (printWhyNot)
      messerr("Arguments 'irow'(%d) and 'icol'(%d) should be equal",irow,icol);
    return false;
  }
  return true;
}

void MatrixSquareDiagonal::addScalar(double v)
{
  for (int irow = 0; irow < getNRows(); irow++)
    _setValue(irow, irow, _getValue(irow, irow) + v);
}

int MatrixSquareDiagonal::_solve(const VectorDouble& b, VectorDouble& x) const
{
  for (int irow=0; irow<getNRows(); irow++)
  {
    double pivot = getValue(irow,irow);
    if (ABS(pivot) < EPSILON10) return 1;
    x[irow] = b[irow] / pivot;
  }
  return 0;
}

/*! Set the contents of a Column */
void MatrixSquareDiagonal::setColumn(int /*icol*/, const VectorDouble& /*tab*/)
{
  my_throw("This function does not make sense for Diagonal Matrix");
}

/*! Set the contents of a Row */
void MatrixSquareDiagonal::setRow(int /*irow*/, const VectorDouble& /*tab*/)
{
  my_throw("This function does not make sense for Diagonal Matrix");
}

String MatrixSquareDiagonal::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  sstr << "- Number of rows    = " <<  getNRows() << std::endl;
  sstr << "- Number of columns = " <<  getNCols() << std::endl;
  sstr << toMatrixDiagonal(String(), VectorString(), VectorString(),
                           getNCols(), getValues());
  return sstr.str();
}

bool MatrixSquareDiagonal::_isPhysicallyPresent(int irow, int icol) const
{
  if (irow != icol) return false;
  return true;
}

/**
 * Converts a VectorVectorDouble into a Matrix
 * Note: the input argument is stored by row (if coming from [] specification)
 * @param X Input VectorVectorDouble argument
 * @return The returned matrix
 *
 * @remark: the matrix is transposed implicitly while reading
 */
MatrixSquareDiagonal* MatrixSquareDiagonal::createFromVVD(const VectorVectorDouble& X)
{
  int nrow = (int) X.size();
  int ncol = (int) X[0].size();
  MatrixRectangular* mattemp = new MatrixRectangular(nrow, ncol);
  if (mattemp->isDiagonal())
  {
    messerr("The matrix does not seem to be Square and Diagonal");
    delete mattemp;
    return nullptr;
  }
  delete mattemp;

  MatrixSquareDiagonal* mat = new MatrixSquareDiagonal(nrow);
  mat->_fillFromVVD(X);
  return mat;
}

