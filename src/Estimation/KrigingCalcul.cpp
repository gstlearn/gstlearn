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
#include "Estimation/KrigingCalcul.hpp"

KrigingCalcul::KrigingCalcul(const MatrixSquareSymmetric* C,
                             const MatrixRectangular* X)
  : _C(C)
  , _X(X)
  , _C0(nullptr)
  , _X0(nullptr)
  , _Z()
  , _beta()
  , _lambda()
  , _Cm1(nullptr)
  , _neq(0)
  , _nbfl(0)
  , _nrhs(0)
{
}

KrigingCalcul::~KrigingCalcul()
{
  delete _Cm1;
}

int KrigingCalcul::setC(const MatrixSquareSymmetric* C)
{
  if (C == nullptr)
  {
    _C = nullptr;
  }
  else
  {
    _C   = C;
    _neq = C->getNRows();
  }
  return 0;
}

bool KrigingCalcul::_matchDimensions(const AMatrix* mat, int nrowsRef, int ncolsRef)
{
  int nrows = mat->getNRows();
  int ncols = mat->getNCols();
  if (nrowsRef > 0 && nrows != nrowsRef)
  {
    messerr("Number of Rows of 'mat' (%d) should match number of rows of 'C' (%d)",
      nrows, nrowsRef);
    return false;
  }
  if (ncolsRef > 0 && ncols != ncolsRef)
  {
    messerr("Number of Columns of 'mat' (%d) should match number of columns of 'C' (%d)",
            ncols, ncolsRef);
    return false;
  }
  return true;
}

int KrigingCalcul::setX(const MatrixRectangular* X)
{
  if (X == nullptr)
  {
    _X = nullptr;
  }
  else
  {
    if (!_matchDimensions(X, _neq, _nbfl)) return 1;
    _X    = X;
    _nbfl = X->getNCols();
  }
  return 0;
}

int KrigingCalcul::setC0(const MatrixRectangular* C0)
{
  if (C0 == nullptr)
  {
    _C0 = nullptr;
  }
  else
  {
    if (!_matchDimensions(C0, _neq, _nrhs)) return 1;
    _C0   = C0;
    _nrhs = C0->getNCols();
  }
  return 0;
}

int KrigingCalcul::setX0(const MatrixRectangular* X0)
{
  if (X0 == nullptr)
  {
    _X0 = nullptr;
  }
  else
  {
    if (!_matchDimensions(X0, _neq, _nrhs)) return 1;
    _X0   = X0;
    _nrhs = X0->getNCols();
  }
  return 0;
}

int KrigingCalcul::setZ(const VectorDouble& Z)
{
  int size = (int)Z.size();
  if (_neq > 0 && _neq != size)
  {
    messerr("Dimension of 'Z' (%d) should match (%d)", size, _neq);
    return 1;
  }
  _Z = Z;
  return 0;
}

int KrigingCalcul::setBeta(const VectorDouble& beta)
{
  int size = (int)beta.size();
  if (_nrhs > 0 && _nrhs != size)
  {
    messerr("Dimension of 'Beta' (%d) should match (%d)", size, _nrhs);
    return 1;
  }
  _beta = beta;
  return 0;
}

int KrigingCalcul::_invertC()
{
  if (_C == nullptr) return 1;
  _Cm1 = _C->clone();
  _Cm1->invert();
  return 0;
}
