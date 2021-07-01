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
#include "MatrixC/MatrixCSDiagCst.hpp"
#include "MatrixC/AMatrixCSquare.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "geoslib_f.h"

MatrixCSDiagCst::MatrixCSDiagCst(int nrow, bool sparse)
  : AMatrixCSquare(nrow, sparse)
  , _cstDiagMatrix(0.)
{
}

MatrixCSDiagCst::MatrixCSDiagCst(const MatrixCSDiagCst &r) 
  : AMatrixCSquare(r)
  , _cstDiagMatrix(r._cstDiagMatrix)
{
}

MatrixCSDiagCst& MatrixCSDiagCst::operator= (const MatrixCSDiagCst &r)
{
  if (this != &r)
  {
    AMatrixCSquare::operator=(r);
    _cstDiagMatrix = r._cstDiagMatrix;
  }
  return *this;
}

MatrixCSDiagCst::~MatrixCSDiagCst()
{
}

IClonable* MatrixCSDiagCst::clone() const
{
  return new MatrixCSDiagCst(*this);
}

double MatrixCSDiagCst::_getValue(int irow, int icol) const
{
  if (irow == icol)
    return _cstDiagMatrix;
  else
    return 0.;
}

double MatrixCSDiagCst::_getValue(int irank) const
{
  return _cstDiagMatrix;
}

void MatrixCSDiagCst::_setValue(int irow, int icol, double value)
{
  if (irow == icol)
    _cstDiagMatrix = value;
  else
  {
    if (abs(value) > EPSILON10)
      my_throw("Attempt to assign a non-zero value to a non-diagonal term");
  }
}

void MatrixCSDiagCst::_setValue(int irank, double value)
{
  _cstDiagMatrix = value;
}

double MatrixCSDiagCst::_determinant() const
{
  double deter = 1.;
  for (int irow=0; irow<getNRows(); irow++)
    deter *= _cstDiagMatrix;
  return deter;
}

void MatrixCSDiagCst::_prodVector(const double *in, double *out) const
{
  int nrow = getNRows();
  for (int irow = 0; irow < nrow; irow++)
    out[irow] = in[irow] * _cstDiagMatrix;
}

void MatrixCSDiagCst::transposeInPlace()
{
  // Nothing should be done
  return;
}

int MatrixCSDiagCst::_invert()
{
  if (_cstDiagMatrix == 0) return 1;
  _cstDiagMatrix = 1. / _cstDiagMatrix;
  return 0;
}

double& MatrixCSDiagCst::_getValueRef(int irow, int icol)
{
  _checkValidIndex(irow,icol);
  return _cstDiagMatrix;
}

void MatrixCSDiagCst::_checkValidIndex(int irow, int icol) const
{
  AMatrixC::_isIndexValid(irow,icol);
  if (irow != icol)
    my_throw("Argument 'irow' and 'icol' should be equal for Diagonal Matrix");
}

void MatrixCSDiagCst::_setValues(const double* values, bool byCol)
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

bool MatrixCSDiagCst::isValid(int irow, int icol, bool printWhyNot) const
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

void MatrixCSDiagCst::addScalar(double v)
{
  my_throw("This function does not make sense for Diagonal Matrix");
}

void MatrixCSDiagCst::addScalarDiag(double v)
{
  _cstDiagMatrix += v;
}

int MatrixCSDiagCst::_solve(const VectorDouble& b, VectorDouble& x) const
{
  if (ABS(_cstDiagMatrix) < EPSILON10) return 1;
  for (int rank=0; rank<(int) b.size(); rank++)
    x[rank] = b[rank] / _cstDiagMatrix;
  return 0;
}

/*! Set the contents of a Column */
void MatrixCSDiagCst::setColumn(int icol, const VectorDouble& tab)
{
  my_throw("This function does not make sense for Diagonal Matrix");
}

/*! Set the contents of a Row */
void MatrixCSDiagCst::setRow(int irow, const VectorDouble& tab)
{
  my_throw("This function does not make sense for Diagonal Constant Matrix");
}

void MatrixCSDiagCst::setDiagonal(const VectorDouble& tab)
{
  my_throw("This function does not make sense for Diagonal Constant Matrix");
}

String MatrixCSDiagCst::toString(int level) const
{
  std::stringstream sstr;

  if (isSparse())
  {
    sstr << AMatrixC::toString(level);
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
