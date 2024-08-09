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
#include "Matrix/AMatrixDense.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"

KrigingCalcul::KrigingCalcul(const MatrixSquareSymmetric* C00,
                             const MatrixSquareSymmetric* C,
                             const MatrixRectangular* X)
  : _C00(C00)
  , _C(C)
  , _X(X)
  , _C0(nullptr)
  , _X0(nullptr)
  , _Z()
  , _Xi()
  , _Beta()
  , _Zstar()
  , _Cm1(nullptr)
  , _lambdaSK(nullptr)
  , _lambdaUK(nullptr)
  , _XtCm1(nullptr)
  , _LambdaSKtX(nullptr)
  , _X0mLambdaSKtX(nullptr)
  , _var(nullptr)
  , _varZ(nullptr)
  , _varBeta(nullptr)
  , _neq(0)
  , _nbfl(0)
  , _nrhs(0)
{
}

KrigingCalcul::~KrigingCalcul()
{
  delete _Cm1;
  delete _lambdaSK;
  delete _lambdaUK;
  delete _XtCm1;
  delete _LambdaSKtX;
  delete _X0mLambdaSKtX;
  delete _var;
  delete _varZ;
  delete _varBeta;
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
  _Beta = beta;
  return 0;
}

int KrigingCalcul::setBayes(const VectorDouble& mBayes,
                            const MatrixSquareSymmetric* SBayes)
{
  if (_nbfl > 0)
  {
    int size = (int)mBayes.size();
    if (_nbfl != size)
    {
      messerr("Dimension of 'mBayes' (%d) should match (%d)", size, _nbfl);
      return 1;
    }
    int sizeM = SBayes->getNRows();
    if (_nbfl != sizeM)
    {
      messerr("Dimension of 'SBayes' (%d) should match (%d)", sizeM, _nbfl);
      return 1;
    }
  }
  _mBayes = mBayes;
  _SBayes = SBayes;
  return 0;
}

int KrigingCalcul::_computeCm1()
{
  if (_C == nullptr) return 1;
  _Cm1 = _C->clone();
  _Cm1->invert();
  return 0;
}

/**
 * @brief Compute the drift value at data locations
 */
int KrigingCalcul::_computeXi()
{
  if (_Beta.empty()) return 1;
  _Xi = _X->prodMatVec(_Beta);
  return 0;
}

int KrigingCalcul::_needCm1()
{
  if (_Cm1 != nullptr) return 0;
  return _computeCm1();
}

int KrigingCalcul::_computeLambdaSK()
{
  if (_C0 == nullptr) return 1;
  if (_needCm1()) return 1;
  if (_lambdaSK == nullptr) _lambdaSK = new MatrixRectangular(_neq, _nrhs);
  _lambdaSK->prodMatMatInPlace(_Cm1, _C0);
  return 0;
}

int KrigingCalcul::_needLambdaSK()
{
  if (_lambdaSK != nullptr) return 0;
  return _computeLambdaSK();
}

int KrigingCalcul::_needLambdaUK()
{
  if (_lambdaUK != nullptr) return 0;
  return _computeLambdaUK();
}

int KrigingCalcul::_computeZstar(bool flagSK)
{
  if (_Z.empty()) return 1;
  if (flagSK)
  {
    if (_needLambdaSK()) return 1;
    _Zstar = _lambdaSK->prodMatVec(_Z);
  }
  else
  {
    if (_needLambdaUK()) return 1;
    _Zstar = _lambdaUK->prodMatVec(_Z);
  }
  return 0;
}

int KrigingCalcul::_needLambdaSKtX()
{
  if (_LambdaSKtX != nullptr) return 0;
  return _computeLambdaSKtX();
}

int KrigingCalcul::_needVarBeta()
{
  if (_varBeta != nullptr) return 0;
  return _computeBeta();
}

int KrigingCalcul::_computeVarianceZ(bool flagSK)
{
  if (_C0 == nullptr) return 1;
  if (_needLambdaSK()) return 1;

  if (_varZ == nullptr) _varZ = new MatrixSquareSymmetric(_nrhs);
  _varZ->prodMatMatInPlace(_lambdaSK, _C0, true, false);

  if (flagSK) return 0;
  if (_X0 == nullptr) return 0;

  if (_needLambdaSKtX()) return 1;
  if (_needVarBeta()) return 1;

  MatrixSquareSymmetric* p1 = new MatrixSquareSymmetric(_nrhs);
  p1->prodNormMatMatInPlace(*_LambdaSKtX, *_varBeta);
  MatrixSquareSymmetric* p2 = new MatrixSquareSymmetric(_nrhs);
  p2->prodNormMatMatInPlace(*_X0, *_varBeta);

  _varZ->linearCombination(1., _varZ, -1., p1, +1., p2);

  delete p1;
  delete p2;
  return 0;
}

int KrigingCalcul::_needVarZ(bool flagSK)
{
  if (_varZ != nullptr) return 0;
  return _computeVarianceZ(flagSK);
}

int KrigingCalcul::_needX0mLambdaSKtX()
{
  if (_X0mLambdaSKtX != nullptr) return 0;
  return _computeX0mLambdaSKtX();
}

int KrigingCalcul::_computeVariance(bool flagSK)
{
  if (_C00 == nullptr) return 1;
  if (_needVarZ(false)) return 1;

  if (_var == nullptr) _var = new MatrixSquareSymmetric(_nrhs);
  _var->linearCombination(1., _C00, -1., _varZ);

  if (flagSK) return 0;

  if (_needX0mLambdaSKtX()) return 1;
  if (_needVarBeta()) return 1;

  MatrixSquareSymmetric* p1 = new MatrixSquareSymmetric(_nrhs);
  p1->prodNormMatMatInPlace(*_X0mLambdaSKtX, *_varBeta);

  _var->linearCombination(1., _var, 1., p1);

  delete p1;
  return 0;
}

int KrigingCalcul::_computeXtCm1()
{
  if (_X == nullptr) return 1;
  if (_needCm1()) return 1;

  if (_XtCm1 == nullptr) _XtCm1 = new MatrixRectangular(_nbfl, _neq);
  _XtCm1->prodMatMatInPlace(_X, _Cm1, true, false);
  return 0;
}

int KrigingCalcul::_needXtCm1()
{
  if (_XtCm1 != nullptr) return 0;
  return _computeXtCm1();
}

int KrigingCalcul::_computeVarianceBeta()
{
  if (_X == nullptr) return 1;
  if (_needXtCm1()) return 1;

  if (_varBeta == nullptr) _varBeta = new MatrixSquareSymmetric(_nbfl);
  _varBeta->prodMatMatInPlace(_XtCm1, _X);
  if (_varBeta->invert()) return 1;
  return 0;
}

int KrigingCalcul::_computeBeta()
{
  if (_Z.empty()) return 1;
  if (_needVarBeta()) return 1;
  if (_needXtCm1()) return 1;

  VectorDouble XtCm1Z = _XtCm1->prodMatVec(_Z);
  _Beta               = _varBeta->prodMatVec(XtCm1Z);
  return 0;
}

int KrigingCalcul::_computeLambdaSKtX()
{
  if (_X == nullptr) return 1;
  if (_needLambdaSK()) return 1;

  if (_LambdaSKtX == nullptr) _LambdaSKtX = new MatrixRectangular(_nbfl, _nrhs);
  _LambdaSKtX->prodMatMatInPlace(_lambdaSK, _X, true, false);
  return 0;
}

int KrigingCalcul::_computeX0mLambdaSKtX()
{
  if (_X0 == nullptr) return 1;
  if (_needLambdaSKtX()) return 1;

  if (_X0mLambdaSKtX == nullptr)
    _X0mLambdaSKtX = new MatrixRectangular(_nbfl, _nrhs);
  _X0mLambdaSKtX->linearCombination(1., _X0, -1., _LambdaSKtX);
  return 0;
}

int KrigingCalcul::_computeLambdaUK()
{
  if (_needXtCm1()) return 1;
  if (_needVarBeta()) return 1;
  if (_needX0mLambdaSKtX()) return 1;

  MatrixRectangular* p1 = new MatrixRectangular(_neq, _nbfl);
  p1->prodMatMatInPlace(_XtCm1, _varBeta, true, false);
  MatrixRectangular* p2 = new MatrixRectangular(_neq, _nbfl);
  p2->prodMatMatInPlace(p1, _X0mLambdaSKtX, false, true);

  if (_lambdaUK == nullptr) _lambdaUK = new MatrixRectangular(_neq, _nrhs);
  _lambdaUK->linearCombination(1., _lambdaSK, 1., p2);

  delete p1;
  delete p2;
  return 0;
}
