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
                             const MatrixSquareSymmetric* Sigma00,
                             const VectorDouble& Beta)
  : _Sigma00(nullptr)
  , _Sigma(nullptr)
  , _Sigma0(nullptr)
  , _X(nullptr)
  , _X0(nullptr)
  , _PriorCov()

  , _Z()
  , _PriorMean()
  , _Zstar()
  , _Beta()
  , _LambdaSK(nullptr)
  , _LambdaUK(nullptr)
  , _Mu(nullptr)
  , _Stdv(nullptr)
  , _VarSK(nullptr)
  , _VarZSK(nullptr)
  , _VarZUK(nullptr)
  , _XtInvSigma(nullptr)
  , _Y0(nullptr)
  , _InvSigma(nullptr)
  , _Sigmac(nullptr)
  , _InvPriorCov(nullptr)

  , _Sigma0p(nullptr)
  , _X0p(nullptr)
  , _Y0p(nullptr)
  , _InvCCK(nullptr)
  , _X0ptInvCCK(nullptr)
  , _Sigma0pInvCCK(nullptr)

  , _neq(0)
  , _nbfl(0)
  , _nrhs(0)
  , _ncck(0)
  , _flagSK(true)
  , _flagBayes(false)
  , _rankColCok()
{
  setZ(Z);
  setSigma(Sigma);
  setX(X);
  setSigma00(Sigma00);
  setBeta(Beta);
}

KrigingCalcul::~KrigingCalcul()
{
  delete _LambdaSK;
  delete _LambdaUK;
  delete _Mu;
  delete _Stdv;
  delete _VarSK;
  delete _VarZSK;
  delete _VarZUK;
  delete _XtInvSigma;
  delete _Y0;
  delete _InvSigma;
  delete _Sigmac;
  delete _InvPriorCov;

  delete _Sigma0p;
  delete _X0p;
  delete _Y0p;
  delete _InvCCK;
  delete _X0ptInvCCK;
  delete _Sigma0pInvCCK;
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

int KrigingCalcul::setSigma00(const MatrixSquareSymmetric* Sigma00)
{
  if (Sigma00 == nullptr)
  {
    _Sigma00 = Sigma00;
  }
  else
  {
    if (!_checkDimensionMatrix("Sigma00", Sigma00, &_nrhs, &_nrhs)) return 1;
    _Sigma00 = Sigma00;
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
    _Sigma0 = Sigma0;
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

int KrigingCalcul::setColCok(const VectorInt& rankColCok)
{
  _ncck = 0;
  if (rankColCok.empty()) return 0;
  _rankColCok = VH::unique(rankColCok);
  int ncol = (int) _rankColCok.size();
  if (ncol >= _nrhs)
  {
    messerr("All variables may be collocated");
    return 1;
  }
  for (int icol = 0; icol < ncol; icol++)
  {
    int jcol = _rankColCok[icol];
    if (jcol < 0 || jcol >= _nrhs)
    {
      messerr("The rank of Collocated variable (%d) should be smaller than %d",
              jcol, _nrhs);
      return 1;
    }
  }
  _ncck = (int)_rankColCok.size();
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

const MatrixRectangular* KrigingCalcul::getMu()
{
  if (_needMu()) return nullptr;
  return _Mu;
}

const MatrixSquareSymmetric* KrigingCalcul::getPostCov()
{
  if (_needSigmac()) return nullptr;
  return _Sigmac;
}

int KrigingCalcul::_needInvSigma()
{
  if (_InvSigma != nullptr) return 0;
  if (_needSigma()) return 1;
  _InvSigma = _Sigma->clone();

  // Collocated case
  if (_ncck > 0)
  {
    if (_needSigma0p()) return 1;
    if (_needSigma0pInvCCK()) return 1;
    MatrixSquareSymmetric Sigma0pInvCCKSigma0pt(_neq);
    Sigma0pInvCCKSigma0pt.prodMatMatInPlace(_Sigma0pInvCCK, _Sigma0p);
    _InvSigma->linearCombination(1., _InvSigma, -1., &Sigma0pInvCCKSigma0pt);
  }
  
  if (_InvSigma->invert()) return 1;
  return 0;
}

int KrigingCalcul::_needInvPriorCov()
{
  if (_InvPriorCov != nullptr) return 0;
  if (!_isPresentMatrix("PriorCov", _PriorCov)) return 1;
  _InvPriorCov = _PriorCov->clone();
  if (_InvPriorCov->invert()) return 1;
  return 0;
}

int KrigingCalcul::_needZstar()
{
  if (!_Zstar.empty()) return 0;
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
      if (_needY0()) return 1;
      VectorDouble mean = _Y0->prodMatVec(_Beta);
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

int KrigingCalcul::_needVarZSK()
{
  if (_VarZSK != nullptr) return 0;
  if (_needSigma0()) return 1;
  if (_needLambdaSK()) return 1;
  if (_VarZSK == nullptr) _VarZSK = new MatrixSquareSymmetric(_nrhs);
  _VarZSK->prodMatMatInPlace(_LambdaSK, _Sigma0, true, false);
  return 0;
}

int KrigingCalcul::_needVarZUK()
{
  if (_VarZUK != nullptr) return 0;
  if (_needSigma0()) return 1;
  if (_needLambdaUK()) return 1;
  if (_VarZUK == nullptr) _VarZUK = new MatrixSquareSymmetric(_nrhs);
  _VarZUK->prodNormMatMatInPlace(*_LambdaUK, *_Sigma, true);
  return 0;
}

int KrigingCalcul::_needVarSK()
{
  if (_VarSK != nullptr) return 0;
  if (!_isPresentMatrix("Sigma00", _Sigma00)) return 1;
  if (_needVarZSK()) return 1;
  if (_VarSK == nullptr) _VarSK = new MatrixSquareSymmetric(_nrhs);
  _VarSK->linearCombination(1., _Sigma00, -1., _VarZSK);
  return 0;
}

int KrigingCalcul::_needStdv()
{
  if (_Stdv != nullptr) return 0;
  if (!_isPresentMatrix("Sigma00", _Sigma00)) return 1;

  if (_Stdv == nullptr) _Stdv = new MatrixSquareSymmetric(_nrhs);
  if (_needVarSK()) return 1;
  _Stdv->linearCombination(1., _VarSK);

  if (!_flagSK)
  {
    if (_needY0()) return 1;
    if (_needSigmac()) return 1;
    MatrixSquareSymmetric p1(_nrhs);
    p1.prodNormMatMatInPlace(*_Y0, *_Sigmac);
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

int KrigingCalcul::_needX()
{
  if (!_isPresentMatrix("X", _X)) return 1;
  return 0;
}

int KrigingCalcul::_needX0()
{
  if (!_isPresentMatrix("X0", _X0)) return 1;
  return 0;
}

int KrigingCalcul::_needSigma()
{
  if (!_isPresentMatrix("Sigma", _Sigma)) return 1;
  return 0;
}

int KrigingCalcul::_needSigma0()
{
  if (!_isPresentMatrix("Sigma0", _Sigma0)) return 1;
  return 0;
}

int KrigingCalcul::_needXtInvSigma()
{
  if (_XtInvSigma != nullptr) return 0;
  if (_needX()) return 1;
  if (_needInvSigma()) return 1;

  if (_XtInvSigma == nullptr) _XtInvSigma = new MatrixRectangular(_nbfl, _neq);
  _XtInvSigma->prodMatMatInPlace(_X, _InvSigma, true, false);
  return 0;
}

int KrigingCalcul::_needSigmac()
{
  if (_Sigmac != nullptr) return 0;
  if (_needX()) return 1;
  if (_needXtInvSigma()) return 1;

  if (_Sigmac == nullptr) _Sigmac = new MatrixSquareSymmetric(_nbfl);
  _Sigmac->prodMatMatInPlace(_XtInvSigma, _X);

  // Bayesian case
  if (_flagBayes)
  {
    if (_needInvPriorCov()) return 1;
    _Sigmac->linearCombination(1., _Sigmac, 1., _InvPriorCov);
  }

  // Compute the inverse matrix
  if (_Sigmac->invert()) return 1;
  return 0;
}

int KrigingCalcul::_needInvCCK()
{
  if (_ncck <= 0) return 1;
  if (!_isPresentMatrix("Sigma00", _Sigma00)) return 1;
  if (!_isPresentIVector("rankColCok", _rankColCok)) return 1;
  if (_InvCCK == nullptr)
  {
    _InvCCK = MatrixSquareSymmetric::sample(_Sigma00, _rankColCok);
    if (_InvCCK->invert()) return 1;
  }
  return 0;
}

int KrigingCalcul::_needSigma0p()
{
  if (_ncck <= 0) return 1;
  if (_needSigma0()) return 1;
  if (!_isPresentIVector("rankColCok", _rankColCok)) return 1;

  if (_Sigma0p == nullptr)
    _Sigma0p = MatrixRectangular::sample(_Sigma0, VectorInt(), _rankColCok);
  return 0;
}

int KrigingCalcul::_needSigma0pInvCCK()
{
  if (_needSigma0p()) return 1;
  if (_needInvCCK()) return 1;

  if (_Sigma0pInvCCK == nullptr) _Sigma0pInvCCK = new MatrixRectangular(_neq, _ncck);
  _Sigma0pInvCCK->prodMatMatInPlace(_Sigma0p, _InvCCK);
  return 0;
}

int KrigingCalcul::_needX0p()
{
  if (_ncck <= 0) return 1;
  if (_needX0()) return 1;
  if (!_isPresentIVector("rankColCok", _rankColCok)) return 1;

  if (_X0p == nullptr)
    _X0p = MatrixRectangular::sample(_X0, _rankColCok, VectorInt());
  return 0;
}

int KrigingCalcul::_needX0ptInvCCK()
{
  if (_needX0p()) return 1;
  if (_needInvCCK()) return 1;

  if (_X0ptInvCCK == nullptr) _X0ptInvCCK = new MatrixRectangular(_nbfl, _ncck);
  _X0ptInvCCK->prodNormMatMatInPlace(*_X0p, *_InvCCK, true);
  return 0;
}

int KrigingCalcul::_needBeta()
{
  if (!_Beta.empty()) return 0;
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

int KrigingCalcul::_needY0()
{
  if (_needX0()) return 1;
  if (_needLambdaSK()) return 1;

  if (_Y0 == nullptr)
    _Y0 = new MatrixRectangular(_nrhs, _nbfl);

  MatrixRectangular LambdaSKtX(_nrhs, _nbfl);
  LambdaSKtX.prodMatMatInPlace(_LambdaSK, _X, true, false);
  _Y0->linearCombination(1., _X0, -1., &LambdaSKtX);
  return 0;
}

int KrigingCalcul::_needMu()
{
  if (_flagSK) return 0;
  if (_needSigmac()) return 1;
  if (_needY0()) return 1;
  if (_Mu == nullptr) _Mu = new MatrixRectangular(_nbfl, _nrhs);
  _Mu->prodMatMatInPlace(_Sigmac, _Y0);
  return 0;
}

int KrigingCalcul::_needLambdaSK()
{
  if (_LambdaSK != nullptr) return 0;
  if (_needSigma0()) return 1;
  if (_needInvSigma()) return 1;
  if (_LambdaSK == nullptr) _LambdaSK = new MatrixRectangular(_neq, _nrhs);
  _LambdaSK->prodMatMatInPlace(_InvSigma, _Sigma0);
  return 0;
}

int KrigingCalcul::_needLambdaUK()
{
  if (_LambdaUK != nullptr) return 0;
  if (_needXtInvSigma()) return 1;
  if (_needMu()) return 1;

  MatrixRectangular p1(_neq, _nrhs);
  p1.prodMatMatInPlace(_XtInvSigma, _Mu, true, false);

  if (_LambdaUK == nullptr) _LambdaUK = new MatrixRectangular(_neq, _nrhs);
  _LambdaUK->linearCombination(1., _LambdaSK, 1., &p1);

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
  _printMatrix("Sigma00", _Sigma00);
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
  _printMatrix("_Y0", _Y0);
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

bool KrigingCalcul::_isPresentIVector(const String& name,
                                      const VectorInt& vec)
{
  if (!vec.empty()) return true;
  messerr(">>> Vector %s is missing (required)", name.c_str());
  return false;
}

int KrigingCalcul::_needY0p()
{
  if (_needX()) return 1;
  if (_needX0()) return 1;
  if (_needXtInvSigma()) return 1;

  if (_Y0p == nullptr) _Y0p = new MatrixRectangular(_ncck, _nbfl);
  _Y0p->prodMatMatInPlace(_Sigma0p, _XtInvSigma, true, true);
  _Y0p->linearCombination(1., _X0p, -1., _Y0p);
  return 0;
}

int KrigingCalcul::_needLambda0()
{
  if (_ncck <= 0) return 1;
  if (!_isPresentMatrix("Sigma00", _Sigma00)) return 1;
  if (!_isPresentIVector("rankColCok", _rankColCok)) return 1;
  if (_needInvSigma()) return 1;
  if (_needSigmac()) return 1;

  MatrixSquareSymmetric* C22 = MatrixSquareSymmetric::sample(_Sigma00, _rankColCok);
  MatrixRectangular* C12 = MatrixRectangular::sample(_Sigma00, _rankColCok, VectorInt());

  MatrixRectangular Sigma0ptInvSigma(_ncck, _neq);
  Sigma0ptInvSigma.prodMatMatInPlace(_Sigma0p, _InvSigma, true);

  MatrixRectangular Y0pSigmac(_ncck, _nbfl);
  Y0pSigmac.prodMatMatInPlace(_Y0p, _Sigmac);

  return 0;
  }