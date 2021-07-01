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
#include "MatrixC/MatrixCSDiag.hpp"
#include "MatrixC/AMatrixCSquare.hpp"
#include "geoslib_f.h"

#include "Basic/AException.hpp"

MatrixCSDiag::MatrixCSDiag(int nrows, bool sparse)
  : AMatrixCSquare(nrows, sparse)
  , _diagMatrix()
{
  _allocate();
}

MatrixCSDiag::MatrixCSDiag(const MatrixCSDiag &r) 
  : AMatrixCSquare(r)
{
  _recopy(r);
}

MatrixCSDiag& MatrixCSDiag::operator= (const MatrixCSDiag &r)
{
  if (this != &r)
  {
    _deallocate();
    AMatrixCSquare::operator=(r);
    _recopy(r);
  }
  return *this;
}

MatrixCSDiag::~MatrixCSDiag()
{
  _deallocate();
}

IClonable* MatrixCSDiag::clone() const
{
  return new MatrixCSDiag(*this);
}

double MatrixCSDiag::_getValue(int irow, int icol) const
{
  if (irow == icol)
    return _diagMatrix[irow];
  else
    return 0.;
}

double MatrixCSDiag::_getValue(int irank) const
{
  _isRankValid(irank);
  return _diagMatrix[irank];
}

void MatrixCSDiag::_setValue(int irow, int icol, double value)
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

void MatrixCSDiag::_setValue(int irank, double value)
{
  _isRankValid(irank);
  _diagMatrix[irank] = value;
}

void MatrixCSDiag::_prodVector(const double *in, double *out) const
{
  int nrow = getNRows();
  for (int irow=0; irow<nrow; irow++)
    out[irow] = in[irow] * _diagMatrix[irow];
}

void MatrixCSDiag::transposeInPlace()
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
void MatrixCSDiag::_setValues(const double* values, bool byCol)
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

double MatrixCSDiag::_determinant() const
{
  double deter = 1.;
  for (int irow = 0; irow < getNRows(); irow++)
    deter *= _diagMatrix[irow];
  return deter;
}

int MatrixCSDiag::_invert()
{
  for (int irow=0; irow<getNRows(); irow++) 
  {
    if (_diagMatrix[irow] == 0) return 1;
    _diagMatrix[irow] = 1. / _diagMatrix[irow];
  }
  return 0;
}

double& MatrixCSDiag::_getValueRef(int irow, int icol)
{
  if (! _isIndexValid(irow,icol))
    throw("Function aborted due to previous errors");
  return _diagMatrix[irow];
}

void MatrixCSDiag::_deallocate()
{
}

void MatrixCSDiag::_recopy(const MatrixCSDiag &r)
{
  _diagMatrix = r._diagMatrix;
}

void MatrixCSDiag::_allocate()
{
  _diagMatrix.resize(_getMatrixSize());
}

bool MatrixCSDiag::_isIndexValid(int irow, int icol) const
{
  AMatrixC::_isIndexValid(irow,icol);
  if (irow != icol)
  {
    messerr("Argument 'irow' and 'icol' should be equal for Diagonal Matrix");
    return false;
  }
  return true;
}

int MatrixCSDiag::_getMatrixSize() const
{
  return (getNRows());
}

bool MatrixCSDiag::isValid(int irow, int icol, bool printWhyNot) const
{
  if (! AMatrixC::isValid(irow,icol,printWhyNot)) return false;
  if (irow != icol)
  {
    if (printWhyNot)
      messerr("Arguments 'irow'(%d) and 'icol'(%d) should be equal",irow,icol);
    return false;
  }
  return true;
}

void MatrixCSDiag::addScalar(double v)
{
  my_throw("This function does not make sense for Diagonal Matrix");
}

int MatrixCSDiag::_solve(const VectorDouble& b, VectorDouble& x) const
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
void MatrixCSDiag::setColumn(int icol, const VectorDouble& tab)
{
  my_throw("This function does not make sense for Diagonal Matrix");
}

/*! Set the contents of a Row */
void MatrixCSDiag::setRow(int irow, const VectorDouble& tab)
{
  my_throw("This function does not make sense for Diagonal Matrix");
}

String MatrixCSDiag::toString(int level) const
{
  std::stringstream sstr;

   if (isSparse())
   {
     sstr << AMatrixC::toString(level);
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
