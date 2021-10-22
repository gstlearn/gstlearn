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
#include "Matrix/MatrixSDiagCst.hpp"
#include "Matrix/AMatrixSquare.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "geoslib_f.h"

MatrixSDiagCst::MatrixSDiagCst(int nrow, bool sparse)
  : AMatrixSquare(nrow, sparse)
  , _cstDiagMatrix(0.)
{
}

MatrixSDiagCst::MatrixSDiagCst(const MatrixSDiagCst &r) 
  : AMatrixSquare(r)
  , _cstDiagMatrix(r._cstDiagMatrix)
{
}

MatrixSDiagCst& MatrixSDiagCst::operator= (const MatrixSDiagCst &r)
{
  if (this != &r)
  {
    AMatrixSquare::operator=(r);
    _cstDiagMatrix = r._cstDiagMatrix;
  }
  return *this;
}

MatrixSDiagCst::~MatrixSDiagCst()
{
}

IClonable* MatrixSDiagCst::clone() const
{
  return new MatrixSDiagCst(*this);
}

double MatrixSDiagCst::_getValue(int irow, int icol) const
{
  if (irow == icol)
    return _cstDiagMatrix;
  else
    return 0.;
}

double MatrixSDiagCst::_getValue(int irank) const
{
  return _cstDiagMatrix;
}

void MatrixSDiagCst::_setValue(int irow, int icol, double value)
{
  if (irow == icol)
    _cstDiagMatrix = value;
  else
  {
    if (abs(value) > EPSILON10)
      my_throw("Attempt to assign a non-zero value to a non-diagonal term");
  }
}

void MatrixSDiagCst::_setValue(int irank, double value)
{
  _cstDiagMatrix = value;
}

double MatrixSDiagCst::_determinant() const
{
  double deter = 1.;
  for (int irow=0; irow<getNRows(); irow++)
    deter *= _cstDiagMatrix;
  return deter;
}

void MatrixSDiagCst::_prodVector(const double *in, double *out) const
{
  int nrow = getNRows();
  for (int irow = 0; irow < nrow; irow++)
    out[irow] = in[irow] * _cstDiagMatrix;
}

void MatrixSDiagCst::transposeInPlace()
{
  // Nothing should be done
  return;
}

int MatrixSDiagCst::_invert()
{
  if (_cstDiagMatrix == 0) return 1;
  _cstDiagMatrix = 1. / _cstDiagMatrix;
  return 0;
}

double& MatrixSDiagCst::_getValueRef(int irow, int icol)
{
  if (! _isValidIndex(irow,icol))
    my_throw("Impossible to return the Address of the non-existant element");
  return _cstDiagMatrix;
}

bool MatrixSDiagCst::_isValidIndex(int irow, int icol) const
{
  AMatrix::_isIndexValid(irow,icol);
  if (irow != icol)
  {
    messerr("Argument 'irow' and 'icol' should be equal for Diagonal Matrix");
    return false;
  }
  return true;
}

void MatrixSDiagCst::_setValues(const double* values, bool byCol)
{
  double refval = TEST;
  int ecr = 0;
  for (int icol = 0; icol < getNCols(); icol++)
    for (int irow = 0; irow < getNRows(); irow++, ecr++)
    {
      if (irow != icol)
      {
        if (! FFFF(values[ecr]) && ABS(values[ecr]) > EPSILON10)
          my_throw("Input 'values' does not correspond to a Diagonal matrix");
      }
      else
      {
        if (FFFF(refval))
          refval = values[ecr];
        else
        {
          if (ABS(refval - values[ecr]) > EPSILON10)
            my_throw("Input 'values' does not correspond to a Diagonal matrix");
        }
        setValue(irow, icol, values[ecr]);
      }
    }
}

bool MatrixSDiagCst::isValid(int irow, int icol, bool printWhyNot) const
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

void MatrixSDiagCst::addScalar(double v)
{
  my_throw("This function does not make sense for Diagonal Matrix");
}

void MatrixSDiagCst::addScalarDiag(double v)
{
  _cstDiagMatrix += v;
}

int MatrixSDiagCst::_solve(const VectorDouble& b, VectorDouble& x) const
{
  if (ABS(_cstDiagMatrix) < EPSILON10) return 1;
  for (int rank=0; rank<(int) b.size(); rank++)
    x[rank] = b[rank] / _cstDiagMatrix;
  return 0;
}

/*! Set the contents of a Column */
void MatrixSDiagCst::setColumn(int icol, const VectorDouble& tab)
{
  my_throw("This function does not make sense for Diagonal Matrix");
}

/*! Set the contents of a Row */
void MatrixSDiagCst::setRow(int irow, const VectorDouble& tab)
{
  my_throw("This function does not make sense for Diagonal Constant Matrix");
}

void MatrixSDiagCst::setDiagonal(const VectorDouble& tab)
{
  my_throw("This function does not make sense for Diagonal Constant Matrix");
}

String MatrixSDiagCst::toString(int level) const
{
  std::stringstream sstr;

  if (isSparse())
  {
    sstr << AMatrix::toString(level);
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
