/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "Matrix/MatrixSquareDiagonalCst.hpp"
#include "Matrix/AMatrixSquare.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"

#include <math.h>

MatrixSquareDiagonalCst::MatrixSquareDiagonalCst(int nrow, bool sparse)
  : AMatrixSquare(nrow, sparse)
  , _cstDiagMatrix(0.)
{
}

MatrixSquareDiagonalCst::MatrixSquareDiagonalCst(const MatrixSquareDiagonalCst &r) 
  : AMatrixSquare(r)
  , _cstDiagMatrix(r._cstDiagMatrix)
{
}

MatrixSquareDiagonalCst& MatrixSquareDiagonalCst::operator= (const MatrixSquareDiagonalCst &r)
{
  if (this != &r)
  {
    AMatrixSquare::operator=(r);
    _cstDiagMatrix = r._cstDiagMatrix;
  }
  return *this;
}

MatrixSquareDiagonalCst::~MatrixSquareDiagonalCst()
{
}

double MatrixSquareDiagonalCst::_getValue(int irow, int icol) const
{
  if (irow == icol)
    return _cstDiagMatrix;
  else
    return 0.;
}

double MatrixSquareDiagonalCst::_getValue(int /*irank*/) const
{
  return _cstDiagMatrix;
}

void MatrixSquareDiagonalCst::_setValue(int irow, int icol, double value)
{
  if (irow == icol)
    _cstDiagMatrix = value;
  else
  {
    if (ABS(value) > EPSILON10)
    {
      messerr("Attempt to assign a non-zero value to a non-diagonal term");
      messerr("Element(%d,%d) = %lf",irow,icol,value);
      messerr("Operation is aborted");
      return;
    }
  }
}

void MatrixSquareDiagonalCst::_setValue(int /*irank*/, double value)
{
  _cstDiagMatrix = value;
}

double MatrixSquareDiagonalCst::determinant() const
{
  int nrow = getNRows();
  double deter = pow(_cstDiagMatrix, (double) nrow);
  return deter;
}

void MatrixSquareDiagonalCst::_prodVector(const double *inv, double *outv) const
{
  int nrow = getNRows();
  for (int irow = 0; irow < nrow; irow++)
    outv[irow] = inv[irow] * _cstDiagMatrix;
}

void MatrixSquareDiagonalCst::transposeInPlace()
{
  // Nothing should be done
  return;
}

int MatrixSquareDiagonalCst::_invert()
{
  if (_cstDiagMatrix == 0) return 1;
  _cstDiagMatrix = 1. / _cstDiagMatrix;
  return 0;
}

double& MatrixSquareDiagonalCst::_getValueRef(int irow, int icol)
{
  if (! _isValidIndex(irow,icol))
    my_throw("Impossible to return the Address of the non-existant element");
  return _cstDiagMatrix;
}

bool MatrixSquareDiagonalCst::_isValidIndex(int irow, int icol) const
{
  AMatrix::_isIndexValid(irow,icol);
  if (irow != icol)
  {
    messerr("Argument 'irow' and 'icol' should be equal for Diagonal Matrix");
    return false;
  }
  return true;
}

void MatrixSquareDiagonalCst::_setValues(const double* values, bool /*byCol*/)
{
  double refval = TEST;
  int ecr = 0;
  for (int icol = 0; icol < getNCols(); icol++)
    for (int irow = 0; irow < getNRows(); irow++, ecr++)
    {
      if (irow != icol)
      {
        if (! FFFF(values[ecr]) && ABS(values[ecr]) > EPSILON10)
        {
          messerr("Input 'values' does not correspond to a Diagonal matrix");
          messerr("- Element(%d,%d) = %lf",irow,icol,values[ecr]);
          messerr("Operation is aborted");
          return;
        }
      }
      else
      {
        if (FFFF(refval))
          refval = values[ecr];
        else
        {
          if (ABS(refval - values[ecr]) > EPSILON10)
          {
            messerr("Input 'values' does not correspond to a Diagonal matrix");
            messerr("- Element(%d,%d) = %lf",irow,icol,values[ecr]);
            messerr("Operation is aborted");
            return;
          }
        }
        setValue(irow, icol, values[ecr]);
      }
    }
}

bool MatrixSquareDiagonalCst::isValid(int irow, int icol, bool printWhyNot) const
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

void MatrixSquareDiagonalCst::addScalar(double /*v*/)
{
  my_throw("This function does not make sense for Diagonal Matrix");
}

void MatrixSquareDiagonalCst::addScalarDiag(double v)
{
  _cstDiagMatrix += v;
}

int MatrixSquareDiagonalCst::_solve(const VectorDouble& b, VectorDouble& x) const
{
  if (ABS(_cstDiagMatrix) < EPSILON10) return 1;
  for (int rank=0; rank<(int) b.size(); rank++)
    x[rank] = b[rank] / _cstDiagMatrix;
  return 0;
}

/*! Set the contents of a Column */
void MatrixSquareDiagonalCst::setColumn(int /*icol*/, const VectorDouble& /*tab*/)
{
  my_throw("This function does not make sense for Diagonal Matrix");
}

/*! Set the contents of a Row */
void MatrixSquareDiagonalCst::setRow(int /*irow*/, const VectorDouble& /*tab*/)
{
  my_throw("This function does not make sense for Diagonal Constant Matrix");
}

void MatrixSquareDiagonalCst::setDiagonal(const VectorDouble& /*tab*/)
{
  my_throw("This function does not make sense for Diagonal Constant Matrix");
}

String MatrixSquareDiagonalCst::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  if (isSparse())
  {
    sstr << AMatrix::toString(strfmt);
  }
  else
  {
    sstr << "- Number of rows    = " << getNRows() << std::endl;
    sstr << "- Number of columns = " << getNCols() << std::endl;
    sstr << toMatrixDiagCst(String(), VectorString(), VectorString(), getNCols(),
                            getValues());
  }
  return sstr.str();
}

bool MatrixSquareDiagonalCst::_isPhysicallyPresent(int irow, int icol) const
{
  if (irow != 0 || icol != 0) return false;
  return true;
}
