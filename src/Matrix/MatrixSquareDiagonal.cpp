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
#include "Matrix/MatrixSquareDiagonal.hpp"
#include "Matrix/AMatrixSquare.hpp"
#include "geoslib_f.h"

#include "Basic/AException.hpp"

MatrixSquareDiagonal::MatrixSquareDiagonal(int nrows, bool sparse)
  : AMatrixSquare(nrows, sparse)
  , _diagMatrix()
{
  _allocate();
}

MatrixSquareDiagonal::MatrixSquareDiagonal(const MatrixSquareDiagonal &r) 
  : AMatrixSquare(r)
{
  _recopy(r);
}

MatrixSquareDiagonal& MatrixSquareDiagonal::operator= (const MatrixSquareDiagonal &r)
{
  if (this != &r)
  {
    _deallocate();
    AMatrixSquare::operator=(r);
    _recopy(r);
  }
  return *this;
}

MatrixSquareDiagonal::~MatrixSquareDiagonal()
{
  _deallocate();
}

IClonable* MatrixSquareDiagonal::clone() const
{
  return new MatrixSquareDiagonal(*this);
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
      my_throw("Attempt to assign a non-zero value to a non-diagonal term");
  }
}

void MatrixSquareDiagonal::_setValue(int irank, double value)
{
  _isRankValid(irank);
  _diagMatrix[irank] = value;
}

void MatrixSquareDiagonal::_prodVector(const double *in, double *out) const
{
  int nrow = getNRows();
  for (int irow=0; irow<nrow; irow++)
    out[irow] = in[irow] * _diagMatrix[irow];
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
 * @param byCol  True if values in 'tab' are sorted by Column
 */
void MatrixSquareDiagonal::_setValues(const double* values, bool byCol)
{
  int ecr = 0;
  for (int icol = 0; icol < getNCols(); icol++)
    for (int irow = 0; irow < getNRows(); irow++, ecr++)
    {
      if (irow != icol)
      {
        if (ABS(values[ecr]) > EPSILON10)
          my_throw("Input 'values' does not correspond to a Diagonal matrix");
      }
      else
      {
        setValue(irow, icol, values[ecr]);
      }
    }
}

double MatrixSquareDiagonal::_determinant() const
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
  _diagMatrix.resize(_getMatrixSize());
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
  my_throw("This function does not make sense for Diagonal Matrix");
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
void MatrixSquareDiagonal::setColumn(int icol, const VectorDouble& tab)
{
  my_throw("This function does not make sense for Diagonal Matrix");
}

/*! Set the contents of a Row */
void MatrixSquareDiagonal::setRow(int irow, const VectorDouble& tab)
{
  my_throw("This function does not make sense for Diagonal Matrix");
}

String MatrixSquareDiagonal::toString(int level) const
{
  std::stringstream sstr;

   if (isSparse())
   {
     sstr << AMatrix::toString(level);
   }
   else
   {
     sstr << "- Number of rows    = " <<  getNRows() << std::endl;
     sstr << "- Number of columns = " <<  getNCols() << std::endl;
    sstr << toMatrixDiagonal(String(), VectorString(), VectorString(),
                             getNCols(), getValues());
   }
  return sstr.str();
}
