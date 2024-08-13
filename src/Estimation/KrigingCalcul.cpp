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
#include "Basic/VectorHelper.hpp"

KrigingCalcul::KrigingCalcul(const VectorDouble& Z,
                             const MatrixSquareSymmetric* Sigma,
                             const MatrixRectangular* X,
                             const MatrixSquareSymmetric* C00,
                             const VectorDouble& Beta)
  : _C00(nullptr)
  , _Sigma(nullptr)
  , _X(nullptr)
  , _Sigma0(nullptr)
  , _X0(nullptr)
  , _PriorCov()
  , _Z()
  , _PriorMean()
  , _Zstar()
  , _Beta()
  , _LambdaSK(nullptr)
  , _LambdaUK(nullptr)
  , _Stdv(nullptr)
  , _VarSK(nullptr)
  , _VarZSK(nullptr)
  , _VarZUK(nullptr)
  , _XtInvSigma(nullptr)
  , _X0mLambdaSKtX(nullptr)
  , _InvSigma(nullptr)
  , _Sigmac(nullptr)
  , _InvPriorCov(nullptr)
  , _neq(0)
  , _nbfl(0)
  , _nrhs(0)
  , _flagSK(true)
  , _flagBayes(false)
{
  setZ(Z);
  setSigma(Sigma);
  setX(X);
  setC00(C00);
  setBeta(Beta);
}

KrigingCalcul::~KrigingCalcul()
{
  delete _LambdaSK;
  delete _LambdaUK;
  delete _Stdv;
  delete _VarSK;
  delete _VarZSK;
  delete _VarZUK;
  delete _XtInvSigma;
  delete _X0mLambdaSKtX;
  delete _InvSigma;
  delete _Sigmac;
  delete _InvPriorCov;
}

bool KrigingCalcul::_checkDimensionVector(const String& name,
                                          const VectorDouble& vec,
                                          int *sizeRef)
{
  int size = (int)vec.size();
  if (*sizeRef > 0 && size != *sizeRef)
  {
    messerr(
      "Dimension of %s (%d) incorrect: it should be (%d)",
      name.c_str(),size, *sizeRef);
    return false;
  }
  if (size > 0) *sizeRef = size;
  return true;
}

bool KrigingCalcul::_checkDimensionMatrix(const String& name,
                                          const AMatrix* mat,
                                          int* nrowsRef,
                                          int* ncolsRef)
{
  int nrows = mat->getNRows();
  int ncols = mat->getNCols();
  if (*nrowsRef > 0 && nrows != *nrowsRef)
  {
    messerr(
      "Number of Rows of %s (%d) incorrect: it should be (%d)",
      name.c_str(), nrows, *nrowsRef);
    return false;
  }
  if (*ncolsRef > 0 && ncols != *ncolsRef)
  {
    messerr("Number of Columns of %s (%d) incorrect: it should be (%d)",
            name.c_str(), ncols, *ncolsRef);
    return false;
  }
  if (nrows > 0) *nrowsRef = nrows;
  if (ncols > 0) *ncolsRef = ncols;
  return true;
}

int KrigingCalcul::setC00(const MatrixSquareSymmetric* C00)
{
  if (C00 == nullptr)
  {
    _C00 = C00;
  }
  else
  {
    if (!_checkDimensionMatrix("C00", C00, &_nrhs, &_nrhs)) return 1;
    _C00 = C00;
  }
  return 0;
}

int KrigingCalcul::setSigma(const MatrixSquareSymmetric* Sigma)
{
  if (Sigma == nullptr)
  {
    _Sigma = nullptr;
  }
  else
  {
    _Sigma   = Sigma;
    _neq = _Sigma->getNRows();
  }
  return 0;
}

int KrigingCalcul::setX(const MatrixRectangular* X)
{
  if (X == nullptr || X->getNRows() <= 0 || X->getNCols() <= 0)
  {
    _X = nullptr;
    _flagSK = true;
  }
  else
  {
    if (!_checkDimensionMatrix("X", X, &_neq, &_nbfl)) return 1;
    _X    = X;
    _flagSK = (_nbfl <= 0);
  }
  return 0;
}

int KrigingCalcul::setSigma0(const MatrixRectangular* Sigma0)
{
  if (Sigma0 == nullptr)
  {
    _Sigma0 = nullptr;
  }
  else
  {
    if (!_checkDimensionMatrix("Sigma0", Sigma0, &_neq, &_nrhs)) return 1;
    _Sigma0   = Sigma0;
  }
  return 0;
}

int KrigingCalcul::setX0(const MatrixRectangular* X0)
{
  if (X0 == nullptr || X0->getNRows() <= 0 || X0->getNCols() <= 0)
  {
    _X0 = nullptr;
  }
  else
  {
    if (!_checkDimensionMatrix("X0", X0, &_nrhs, &_nbfl)) return 1;
    _X0   = X0;
  }
  return 0;
}

int KrigingCalcul::setZ(const VectorDouble& Z)
{
  if (! _checkDimensionVector("Z", Z, &_neq)) return 1;
  _Z = Z;
  return 0;
}

int KrigingCalcul::setBeta(const VectorDouble& beta)
{
  if (beta.empty()) return 0;
  if (!_checkDimensionVector("Beta", beta, &_nbfl)) return 1;
  _Beta = beta;
  return 0;
}

int KrigingCalcul::setBayes(const VectorDouble& PriorMean,
                            const MatrixSquareSymmetric* PriorCov)
{
  if (!_checkDimensionVector("PriorMean", PriorMean, &_nbfl)) return 1;
  if (!_checkDimensionMatrix("PriorCov", PriorCov, &_nbfl, &_nbfl)) return 1;
  _PriorMean = PriorMean;
  _PriorCov  = PriorCov;
  _flagBayes = true;
  _Beta.clear(); // Beta is cleared in order to provoke its calculation
  return 0;
}

VectorDouble KrigingCalcul::getEstimation()
{
  if (_needZstar()) return VectorDouble();
  return _Zstar;
}

VectorDouble KrigingCalcul::getPostMean()
{
  if (_needBeta()) return VectorDouble();
  return _Beta;
}

VectorDouble KrigingCalcul::getStdv()
{
  if (_needStdv()) return VectorDouble();
  return _Stdv->getDiagonal();
}

const MatrixSquareSymmetric* KrigingCalcul::getStdvMat()
{
  if (_needStdv()) return nullptr;
  return _Stdv;
}

VectorDouble KrigingCalcul::getVarianceZstar()
{
  if (_flagSK)
  {
    if (_needVarZSK()) return VectorDouble();
    return _VarZSK->getDiagonal();
  }
  if (_needVarZUK()) return VectorDouble();
  return _VarZUK->getDiagonal();
}

const MatrixSquareSymmetric* KrigingCalcul::getVarianceZstarMat()
{
  if (_flagSK)
  {
    if (_needVarZSK()) return nullptr;
    return _VarZSK;
  }
  if (_needVarZUK()) return nullptr;
  return _VarZUK;
}

const MatrixRectangular* KrigingCalcul::getLambdaSK()
{
  if (_needLambdaSK()) return nullptr;
  return _LambdaSK;
}

const MatrixRectangular* KrigingCalcul::getLambdaUK()
{
  if (_needLambdaUK()) return nullptr;
  return _LambdaUK;
}

const MatrixSquareSymmetric* KrigingCalcul::getPostCov()
{
  if (_needSigmac()) return nullptr;
  return _Sigmac;
}

int KrigingCalcul::_computeInvSigma()
{
  if (!_isPresentMatrix("Sigma", _Sigma)) return 1;
  _InvSigma = _Sigma->clone();
  _InvSigma->invert();
  return 0;
}

int KrigingCalcul::_computeInvPriorCov()
{
  if (!_isPresentMatrix("PriorCov", _PriorCov)) return 1;
  _InvPriorCov = _PriorCov->clone();
  _InvPriorCov->invert();
  return 0;
}

int KrigingCalcul::_needZstar()
{
  if (!_Zstar.empty()) return 0;
  return _computeZstar();
}

int KrigingCalcul::_needVarSK()
{
  if (_VarSK != nullptr) return 0;
  return _computeVarSK();
}

int KrigingCalcul::_needBeta()
{
  if (!_Beta.empty()) return 0;
  return _computeBeta();
}

int KrigingCalcul::_needInvSigma()
{
  if (_InvSigma != nullptr) return 0;
  return _computeInvSigma();
}

int KrigingCalcul::_computeLambdaSK()
{
  if (!_isPresentMatrix("Sigma0", _Sigma0)) return 1;
  if (_needInvSigma()) return 1;
  if (_LambdaSK == nullptr) _LambdaSK = new MatrixRectangular(_neq, _nrhs);
  _LambdaSK->prodMatMatInPlace(_InvSigma, _Sigma0);
  return 0;
}

int KrigingCalcul::_needLambdaSK()
{
  if (_LambdaSK != nullptr) return 0;
  return _computeLambdaSK();
}

int KrigingCalcul::_needLambdaUK()
{
  if (_LambdaUK != nullptr) return 0;
  return _computeLambdaUK();
}

int KrigingCalcul::_computeZstar()
{
  if (!_isPresentVector("Z", _Z)) return 1;
  if (_flagSK || _flagBayes)
  {
    if (_needLambdaSK()) return 1;
    _Zstar = _LambdaSK->prodMatVec(_Z, true);

    // Adding Mean per Variable
    if (_needBeta()) return 1;
    if (_flagSK)
    {
      VH::linearCombinationInPlace(1., _Zstar, 1., _Beta, _Zstar);
    }

    if (_flagBayes)
    {
      if (_needX0mLambdaSKtX()) return 1;
      VectorDouble mean = _X0mLambdaSKtX->prodMatVec(_Beta);
      VH::linearCombinationInPlace(1., _Zstar, 1., mean, _Zstar);
    }
  }
  else
  {
    if (_needLambdaUK()) return 1;
    _Zstar = _LambdaUK->prodMatVec(_Z, true);
  }
  return 0;
}

int KrigingCalcul::_needSigmac()
{
  if (_Sigmac != nullptr) return 0;
  return _computeSigmac();
}

int KrigingCalcul::_needStdv()
{
  if (_Stdv != nullptr) return 0;
  return _computeStdv();
}

int KrigingCalcul::_needInvPriorCov()
{
  if (_InvPriorCov != nullptr) return 0;
  return _computeInvPriorCov();
}

int KrigingCalcul::_needVarZUK()
{
  if (_VarZUK != nullptr) return 0;
  return _computeVarZUK();
}

int KrigingCalcul::_needVarZSK()
{
  if (_VarZSK != nullptr) return 0;
  return _computeVarZSK();
}

int KrigingCalcul::_needX0mLambdaSKtX()
{
  if (_X0mLambdaSKtX != nullptr) return 0;
  return _computeX0mLambdaSKtX();
}

int KrigingCalcul::_computeVarZSK()
{
  if (!_isPresentMatrix("Sigma0", _Sigma0)) return 1;

  if (_needLambdaSK()) return 1;
  if (_VarZSK == nullptr) _VarZSK = new MatrixSquareSymmetric(_nrhs);
  _VarZSK->prodMatMatInPlace(_LambdaSK, _Sigma0, true, false);
  return 0;
}

int KrigingCalcul::_computeVarZUK()
{
  if (!_isPresentMatrix("Sigma0", _Sigma0)) return 1;
  if (_needLambdaUK()) return 1;
  if (_VarZUK == nullptr) _VarZUK = new MatrixSquareSymmetric(_nrhs);
  _VarZUK->prodNormMatMatInPlace(*_LambdaUK, *_Sigma, true);
  return 0;
}

int KrigingCalcul::_computeVarSK()
{
  if (!_isPresentMatrix("C00", _C00)) return 1;
  if (_needVarZSK()) return 1;
  if (_VarSK == nullptr) _VarSK = new MatrixSquareSymmetric(_nrhs);
  _VarSK->linearCombination(1., _C00, -1., _VarZSK);
  return 0;
}

int KrigingCalcul::_computeStdv()
{
  if (!_isPresentMatrix("C00", _C00)) return 1;

  if (_Stdv == nullptr) _Stdv = new MatrixSquareSymmetric(_nrhs);
  if (_needVarSK()) return 1;
  _Stdv->linearCombination(1., _VarSK);

  if (!_flagSK)
  {
    if (_needX0mLambdaSKtX()) return 1;
    if (_needSigmac()) return 1;
    MatrixSquareSymmetric p1(_nrhs);
    p1.prodNormMatMatInPlace(*_X0mLambdaSKtX, *_Sigmac);
    _Stdv->linearCombination(1., _Stdv, 1., &p1);
  }

  // Transform variance into standard deviation

  for (int irow = 0; irow < _nrhs; irow++)
    for (int icol = 0; icol < _nrhs; icol++)
      if (irow >= icol)
      {
        double value = _Stdv->getValue(irow, icol);
        value        = (value > 0.) ? sqrt(value) : 0.;
        _Stdv->setValue(irow, icol, value);
      }
  return 0;
}

int KrigingCalcul::_computeXtInvSigma()
{
  if (!_isPresentMatrix("X", _X)) return 1;
  if (_needInvSigma()) return 1;

  if (_XtInvSigma == nullptr) _XtInvSigma = new MatrixRectangular(_nbfl, _neq);
  _XtInvSigma->prodMatMatInPlace(_X, _InvSigma, true, false);
  return 0;
}

int KrigingCalcul::_needXtInvSigma()
{
  if (_XtInvSigma != nullptr) return 0;
  return _computeXtInvSigma();
}

int KrigingCalcul::_computeSigmac()
{
  if (!_isPresentMatrix("X", _X)) return 1;
  if (_needXtInvSigma()) return 1;

  if (_Sigmac == nullptr) _Sigmac = new MatrixSquareSymmetric(_nbfl);
  _Sigmac->prodMatMatInPlace(_XtInvSigma, _X);

  if (_flagBayes)
  {
    if (_needInvPriorCov()) return 1;
    _Sigmac->linearCombination(1., _Sigmac, 1., _InvPriorCov);
  }
  if (_Sigmac->invert()) return 1;
  return 0;
}

int KrigingCalcul::_computeBeta()
{
  if (!_isPresentVector("Z", _Z)) return 1;
  if (_needSigmac()) return 1;
  if (_needXtInvSigma()) return 1;

  VectorDouble XtInvCZ = _XtInvSigma->prodMatVec(_Z);

  if (_flagBayes)
  {
    if (!_isPresentVector("PriorMean", _PriorMean)) return 1;
    if (_needInvPriorCov()) return 1;
    VectorDouble InvSMBayes = _InvPriorCov->prodMatVec(_PriorMean);
    VH::linearCombinationInPlace(1., XtInvCZ, 1., InvSMBayes, XtInvCZ);
  }

  _Beta = _Sigmac->prodMatVec(XtInvCZ);
  return 0;
}

int KrigingCalcul::_computeX0mLambdaSKtX()
{
  if (!_isPresentMatrix("X", _X)) return 1;
  if (!_isPresentMatrix("X0", _X0)) return 1;
  if (_needLambdaSK()) return 1;

  if (_X0mLambdaSKtX == nullptr)
    _X0mLambdaSKtX = new MatrixRectangular(_nrhs, _nbfl);

  MatrixRectangular LambdaSKtX(_nrhs, _nbfl);
  LambdaSKtX.prodMatMatInPlace(_LambdaSK, _X, true, false);
  _X0mLambdaSKtX->linearCombination(1., _X0, -1., &LambdaSKtX);
  return 0;
}

int KrigingCalcul::_computeLambdaUK()
{
  if (_needXtInvSigma()) return 1;
  if (_needSigmac()) return 1;
  if (_needX0mLambdaSKtX()) return 1;

  MatrixRectangular p1(_neq, _nbfl);
  p1.prodMatMatInPlace(_XtInvSigma, _Sigmac, true, false);
  MatrixRectangular p2(_neq, _nrhs);
  p2.prodMatMatInPlace(&p1, _X0mLambdaSKtX, false, true);

  if (_LambdaUK == nullptr) _LambdaUK = new MatrixRectangular(_neq, _nrhs);
  _LambdaUK->linearCombination(1., _LambdaSK, 1., &p2);

  return 0;
}

void KrigingCalcul::_printMatrix(const String& name, const AMatrix* mat)
{
  if (mat == nullptr || mat->empty()) return;
  message(" - %s (%d, %d)\n", name.c_str(), mat->getNRows(), mat->getNCols());
}

void KrigingCalcul::_printVector(const String& name, const VectorDouble& vec)
{
  if (vec.empty()) return;
  message(" - %s (%d)\n", name.c_str(), (int)vec.size());
}

void KrigingCalcul::printStatus() const
{
  message("\nGeneral Parameters\n");
  message("Number of Covariance Rows = %d\n", _neq);
  message("Number of Drift equations = %d\n", _nbfl);
  message("Number of Right_Hand sides = %d\n", _nrhs);
  if (_flagSK)
    message("Working with Known Mean(s)\n");
  else
    message("Working with Unknown Mean(s)\n");

  message("\nExternal Pointers\n");
  _printMatrix("C00", _C00);
  _printMatrix("Sigma", _Sigma);
  _printMatrix("Sigma0", _Sigma0);
  _printMatrix("X", _X);
  _printMatrix("X0", _X0);
  _printMatrix("PriorCov", _PriorCov);
  _printVector("Z", _Z);
  _printVector("PriorMean", _PriorMean);

  message("\nInternal Memory\n");
  _printVector("Zstar", _Zstar);
  _printVector("Beta", _Beta);
  _printMatrix("LambdaSK", _LambdaSK);
  _printMatrix("LambdaUK", _LambdaUK);
  _printMatrix("Stdv", _Stdv);
  _printMatrix("VarSK", _VarSK);
  _printMatrix("VarZSK", _VarZSK);
  _printMatrix("VarZUK", _VarZUK);
  _printMatrix("XtInvSigma", _XtInvSigma);
  _printMatrix("X0mLambdaSKtX", _X0mLambdaSKtX);
  _printMatrix("InvSigma", _InvSigma);
  _printMatrix("Sigmac", _Sigmac);
  _printMatrix("InvPriorCov", _InvPriorCov);
}

bool KrigingCalcul::_isPresentMatrix(const String& name, const AMatrix* mat)
{
  if (mat != nullptr) return true;
  messerr(">>> Matrix %s is missing (required)", name.c_str());
  return false;
}

bool KrigingCalcul::_isPresentVector(const String& name,
                                     const VectorDouble& vec)
{
  if (!vec.empty()) return true;
  messerr(">>> Vector %s is missing (required)", name.c_str());
  return false;
}