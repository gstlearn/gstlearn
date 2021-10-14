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
#include "MatrixC/MatrixCSSym.hpp"
#include "MatrixC/AMatrixCSquare.hpp"
#include "Basic/AException.hpp"
#include "geoslib_f.h"
#include "geoslib_old_f.h"

MatrixCSSym::MatrixCSSym(int nrow, bool sparse)
: AMatrixCSquare(nrow, sparse)
, _squareSymMatrix()
{
  _allocate();
}

MatrixCSSym::MatrixCSSym(const MatrixCSSym &r) 
  : AMatrixCSquare(r)
{
  _recopy(r);
}

MatrixCSSym::MatrixCSSym(const AMatrixC &m)
    : AMatrixCSquare(m.getNRows(), m.isSparse()),
      _squareSymMatrix()
{
  if (m.isEmpty())
    my_throw("The input matrix should be non-empty");
  if (!m.isSquare())
    my_throw("The input matrix should be Square");
  if (!m.isSymmetric())
    my_throw("The input matrix should be Symmetric");

  _setNRows(m.getNRows());
  _setNCols(m.getNCols());
  _allocate();
  for (int icol = 0; icol < m.getNCols(); icol++)
    for (int irow = 0; irow < m.getNRows(); irow++)
    {
      setValue(irow,icol,m.getValue(irow,icol));
    }
}

MatrixCSSym& MatrixCSSym::operator= (const MatrixCSSym &r)
{
  if (this != &r)
  {
    _deallocate();
    AMatrixCSquare::operator=(r);
    _recopy(r);
  }
  return *this;
}

MatrixCSSym::~MatrixCSSym()
{
  _deallocate();
}

IClonable* MatrixCSSym::clone() const
{
  return new MatrixCSSym(*this);
}

double MatrixCSSym::_getValue(int irow, int icol) const
{
  _isIndexValid(irow,icol);
  int rank = _getIndexToRank(irow,icol);
  return _squareSymMatrix[rank];
}

double MatrixCSSym::_getValue(int irank) const
{
  return _squareSymMatrix[irank];
}

void MatrixCSSym::_setValue(int irow, int icol, double value)
{
  _isIndexValid(irow, icol);
  int irank = _getIndexToRank(irow, icol);
  _squareSymMatrix[irank] = value;
}

void MatrixCSSym::_setValue(int irank, double value)
{
  _isRankValid(irank);
  _squareSymMatrix[irank] = value;
}

/**
 * Loading values from an Upper Triangular matrix
 * @param nsize
 * @param tab
 */
void MatrixCSSym::initMatTri(int     nsize,
                             double* tab)
{
  _isNumberValid(nsize,nsize);
  _setNSize(nsize);
  _allocate();
  for (int i=0; i<_getMatrixSize(); i++) _squareSymMatrix[i] = tab[i];
}

void MatrixCSSym::_prodVector(const double *in, double *out) const
{
  matrix_triangular_product(getNRows(),2,_squareSymMatrix.data(),in,out);
}

/**
 * \warning : values is provided as a square complete matrix
 */
void MatrixCSSym::_setValues(const double* values, bool byCol)
{
  // Check that the input argument corresponds to a square symmetric matrix
  for (int icol = 0; icol < getNCols(); icol++)
    for (int irow = 0; irow < getNRows(); irow++)
    {
      double val1 = values[icol * getNRows() + irow];
      double val2 = values[irow * getNCols() + icol];
      if (ABS(val1 - val2) > EPSILON10)
        my_throw("'values' must correspond to a Square Symmetric Matrix");
    }

  int ecr = 0;
  for (int icol = 0; icol < getNCols(); icol++)
    for (int irow = 0; irow < getNRows(); irow++, ecr++)
    {
      setValue(irow, icol, values[ecr]);
    }
}

int MatrixCSSym::_invert()
{
  return matrix_invert_triangle(getNRows(),_squareSymMatrix.data(), -1);
}

double& MatrixCSSym::_getValueRef(int irow, int icol)
{
  _isIndexValid(irow,icol);
  int rank = _getIndexToRank(irow,icol);
  return _squareSymMatrix[rank];
}

void MatrixCSSym::_deallocate()
{
}

void MatrixCSSym::_recopy(const MatrixCSSym &r)
{
  _squareSymMatrix = r._squareSymMatrix;
}

void MatrixCSSym::_allocate()
{
  _squareSymMatrix.resize(_getMatrixSize());
}

int MatrixCSSym::_getIndexToRank(int irow,
                                 int icol) const
{
  int rank;

  int n = getNRows();
  if (irow >= icol)
    rank = icol * n + irow - icol * (icol + 1) / 2;
  else
    rank = irow * n + icol - irow * (irow + 1) / 2;
  return rank;
}

int MatrixCSSym::_getMatrixSize() const
{
  int n = getNRows();
  int size = n * (n + 1) / 2;
  return size;
}

int MatrixCSSym::_solve(const VectorDouble& b, VectorDouble& x) const
{
  int pivot;
  return matrix_solve(1,_squareSymMatrix.data(),b.data(),x.data(),
                      static_cast<int> (b.size()),1,&pivot);
}

double MatrixCSSym::_determinant() const
{
  return matrix_determinant(getNRows(),_squareSymMatrix.data());
}

String MatrixCSSym::toString(int level) const
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
     sstr << toMatrixSymmetric(String(), VectorString(), VectorString(),
                               true, getNCols(), getValues());
   }
  return sstr.str();
}
