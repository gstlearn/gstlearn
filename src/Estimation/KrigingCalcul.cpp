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

KrigingCalcul::KrigingCalcul(const VectorDouble* Z,
                             const MatrixSquareSymmetric* Sigma,
                             const MatrixRectangular* X,
                             const MatrixSquareSymmetric* Sigma00,
                             const VectorDouble* BetaRef)
  : _Sigma00(nullptr)
  , _Sigma(nullptr)
  , _Sigma0(nullptr)
  , _X(nullptr)
  , _X0(nullptr)
  , _PriorCov(nullptr)
  , _Z(nullptr)
  , _PriorMean(nullptr)
  , _BetaRef(nullptr)
  , _Zp(nullptr)
  , _rankColCok(nullptr)

  , _Zstar()
  , _Beta()
  , _LambdaSK(nullptr)
  , _LambdaUK(nullptr)
  , _MuUK(nullptr)
  , _Stdv(nullptr)
  , _VarZSK(nullptr)
  , _VarZUK(nullptr)
  , _XtInvSigma(nullptr)
  , _Y0(nullptr)
  , _InvSigmaSigma0(nullptr)
  , _InvSigma(nullptr)
  , _Sigmac(nullptr)
  , _InvPriorCov(nullptr)

  , _Sigma00pp(nullptr)
  , _Sigma00p(nullptr)
  , _Sigma0p(nullptr)
  , _X0p(nullptr)
  , _Y0p(nullptr)
  , _Z0p()
  , _Lambda0(nullptr)

  , _neq(0)
  , _nbfl(0)
  , _nrhs(0)
  , _ncck(0)
  , _flagSK(true)
  , _flagBayes(false)
{
  (void) setData(Z, Sigma, X, Sigma00, BetaRef);
}

KrigingCalcul::~KrigingCalcul()
{
  _resetMemory();
}

void KrigingCalcul::_resetMemory()
{
  _Zstar.clear();
  _Beta.clear();

  delete _LambdaSK;
  _LambdaSK = nullptr;
  delete _LambdaUK;
  _LambdaUK = nullptr;
  delete _MuUK;
  _MuUK = nullptr;
  delete _Stdv;
  _Stdv = nullptr;
  delete _VarZSK;
  _VarZSK = nullptr;
  delete _VarZUK;
  _VarZUK = nullptr;
  delete _XtInvSigma;
  _XtInvSigma = nullptr;
  delete _Y0;
  _Y0 = nullptr;
  delete _InvSigmaSigma0;
  _InvSigmaSigma0 = nullptr;
  delete _InvSigma;
  _InvSigma = nullptr;
  delete _Sigmac;
  _Sigmac = nullptr;
  delete _InvPriorCov;
  _InvPriorCov = nullptr;

  delete _Sigma00pp;
  _Sigma00pp = nullptr;
  delete _Sigma0p;
  _Sigma0p = nullptr;
  delete _X0p;
  _X0p = nullptr;
  delete _Y0p;
  _Y0p = nullptr;
  delete _Lambda0;
  _Lambda0 = nullptr;
}

/**
 * @brief Modify the Data-dependent elements
 *
 * @param Z Data flattened vector (possibly multivariate) 
 * @param Sigma Data-Data Covariance matrix
 * @param X     Data Drift Matrix
 * @param Sigma00  Target-target Variance Matrix
 * @param BetaRef  Vector of known Drift coefficients
 * @return int
 *
 * @note If one element is not provided, its address (if already defined) is
 * @note kept unchanged (even if its contents may have been updated)
 */
int KrigingCalcul::setData(const VectorDouble* Z,
                           const MatrixSquareSymmetric* Sigma,
                           const MatrixRectangular* X,
                           const MatrixSquareSymmetric* Sigma00,
                           const VectorDouble* BetaRef)
{
  _resetMemory();

  // Argument Z
  if (Z != nullptr)
  {
    if (!_checkDimensionVector("Z", Z, &_neq)) return 1;
    _Z = Z;
  }

  // Argument Sigma
  if (Sigma != nullptr)
  {
    if (!_checkDimensionMatrix("Sigma", Sigma, &_neq, &_neq)) return 1;
    _Sigma = Sigma;
  }

  // Argument X
  if (X == nullptr || X->getNRows() <= 0 || X->getNCols() <= 0)
  {
    _X      = nullptr;
    _flagSK = true;
  }
  else
  {
    if (!_checkDimensionMatrix("X", X, &_neq, &_nbfl)) return 1;
    _X      = X;
    _flagSK = (_nbfl <= 0);
  }

  // Argument Sigma00
  if (Sigma00 != nullptr)
  {
    if (!_checkDimensionMatrix("Sigma00", Sigma00, &_nrhs, &_nrhs)) return 1;
    _Sigma00 = Sigma00;
  }

  // Argument BetaRef
  if (BetaRef != nullptr)
  {
    if (!_checkDimensionVector("Beta", BetaRef, &_nbfl)) return 1;
    _BetaRef = BetaRef;
  }
  return 0;
}

int KrigingCalcul::setTarget(const MatrixRectangular* Sigma0,
                             const MatrixRectangular* X0)
{
  _resetMemory();

  // Argument Sigma0
  if (Sigma0 == nullptr)
  {
    _Sigma0 = nullptr;
  }
  else
  {
    if (!_checkDimensionMatrix("Sigma0", Sigma0, &_neq, &_nrhs)) return 1;
    _Sigma0 = Sigma0;
  }

  // Argument X0
  if (X0 == nullptr || X0->getNRows() <= 0 || X0->getNCols() <= 0)
  {
    _X0 = nullptr;
  }
  else
  {
    if (!_checkDimensionMatrix("X0", X0, &_nrhs, &_nbfl)) return 1;
    _X0 = X0;
  }
  return 0;
}

bool KrigingCalcul::_checkDimensionVector(const String& name,
                                          const VectorDouble* vec,
                                          int *sizeRef)
{
  int size = (int)vec->size();
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

int KrigingCalcul::setColCok(const VectorDouble* Zp,
                             const VectorInt* rankColCok)
{
  _resetMemory();

  if (Zp == nullptr || rankColCok == nullptr)
  {
    _ncck = 0;
    return 0;
  }

  // Argument Zp
  if (!_checkDimensionVector("Zp", Zp, &_nrhs)) return 1;

  // Argument rankColCok
  int ncck = (int) rankColCok->size();
  if (ncck >= _nrhs)
  {
    messerr("All variables may be collocated");
    return 1;
  }

  _rankColCok = rankColCok;
  _Zp         = Zp;
  _ncck       = ncck;

  return 0;
}

int KrigingCalcul::setBayes(const VectorDouble* PriorMean,
                            const MatrixSquareSymmetric* PriorCov)
{
  _resetMemory();

  if (PriorMean == nullptr || PriorCov == nullptr)
  {
    _flagBayes = false;
    return 0;
  }
  
  if (!_checkDimensionVector("PriorMean", PriorMean, &_nbfl)) return 1;
  if (!_checkDimensionMatrix("PriorCov", PriorCov, &_nbfl, &_nbfl)) return 1;
  _PriorMean = PriorMean;
  _PriorCov  = PriorCov;
  _flagBayes = true;
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

const MatrixRectangular* KrigingCalcul::getLambda0()
{
  if (_needLambda0()) return nullptr;
  return _Lambda0;
}

const MatrixRectangular* KrigingCalcul::getMuUK()
{
  if (_needMuUK()) return nullptr;
  return _MuUK;
}

const MatrixSquareSymmetric* KrigingCalcul::getPostCov()
{
  if (_needSigmac()) return nullptr;
  return _Sigmac;
}

const MatrixRectangular* KrigingCalcul::getX0()
{
  if (_needX0()) return nullptr;
  return _X0;
}

const MatrixRectangular* KrigingCalcul::getX0p()
{
  if (_needX0p()) return nullptr;
  return _X0p;
}

const MatrixRectangular* KrigingCalcul::getY0()
{
  if (_needY0()) return nullptr;
  return _Y0;
}

const MatrixRectangular* KrigingCalcul::getY0p()
{
  if (_needY0p()) return nullptr;
  return _Y0p;
}

const MatrixRectangular* KrigingCalcul::getSigma0()
{
  if (_needSigma0()) return nullptr;
  return _Sigma0;
}

const MatrixRectangular* KrigingCalcul::getSigma0p()
{
  if (_needSigma0p()) return nullptr;
  return _Sigma0p;
}

int KrigingCalcul::_needInvSigma()
{
  if (_needSigma()) return 1;
  _InvSigma = _Sigma->clone();

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
  if (! _Zstar.empty()) return 0;
  if (!_isPresentVector("Z", _Z)) return 1;
  if (_flagSK || _flagBayes)
  {
    if (_needLambdaSK()) return 1;
    _Zstar = _LambdaSK->prodMatVec(*_Z, true);

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
    _Zstar = _LambdaUK->prodMatVec(*_Z, true);
  }

  // Collocated case
  if (_ncck > 0)
  {
    if (_needZ0p()) return 1;
    VectorDouble Zstar0 = _Lambda0->prodMatVec(_Z0p, true);
    VH::linearCombinationInPlace(1., _Zstar, 1., Zstar0, _Zstar);
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
  _VarZUK->prodNormMatMatInPlace(_LambdaUK, _Sigma, true);

  if (_ncck > 0)
  {
    if (_needLambda0()) return 1;
    if (_needSigma0p()) return 1;
    if (_needSigma00pp()) return 1;
    MatrixRectangular* p1 = new MatrixRectangular(_nrhs, _ncck);
    p1->prodMatMatInPlace(_Lambda0, _Sigma00pp, true);
    MatrixSquareSymmetric* LCLt = new MatrixSquareSymmetric(_nrhs);
    LCLt->prodMatMatInPlace(p1, _Lambda0);
    delete p1;

    MatrixRectangular* p2 = new MatrixRectangular(_nrhs, _ncck);
    p2->prodMatMatInPlace(_Lambda0, _Sigma0p, true, true);
    MatrixSquareSymmetric* LC0tLUK = new MatrixSquareSymmetric(_nrhs);
    LC0tLUK->prodMatMatInPlace(p2, _LambdaUK);
    delete p2;

    _VarZUK->linearCombination(1., _VarZUK, 2., LC0tLUK, 1., LCLt);
    delete LCLt;
    delete LC0tLUK;
  }
  return 0;
}

int KrigingCalcul::_needStdv()
{
  if (_Stdv != nullptr) return 0;
  if (!_isPresentMatrix("Sigma00", _Sigma00)) return 1;

  if (_Stdv == nullptr) _Stdv = new MatrixSquareSymmetric(_nrhs);

  if (_flagSK)
  {
    if (_needVarZSK()) return 1;
    _Stdv->linearCombination(1., _Sigma00, -1., _VarZSK);
  }
  else
  {
    if (_needLambdaUK()) return 1;
    if (_needSigma0()) return 1;
    if (_needMuUK()) return 1;
    if (_needX0()) return 1;
    _Stdv = _Sigma00->clone();
    MatrixRectangular* p1 = new MatrixRectangular(_nrhs, _nrhs);
    p1->prodMatMatInPlace(_LambdaUK, _Sigma0, true);
    MatrixRectangular* p2 = new MatrixRectangular(_nrhs, _nrhs);
    p2->prodMatMatInPlace(_MuUK, _X0, true, true);
    _Stdv->linearCombination(1, _Stdv, -1., p1, +1., p2);
    delete p1;
    delete p2;
  }

  if (_ncck > 0)
  {
    if (_needSigma00p()) return 1;
    MatrixSquareSymmetric* p1 = new MatrixSquareSymmetric(_nrhs);
    p1->prodMatMatInPlace(_Sigma00p, _Lambda0, true);
    _Stdv->linearCombination(1., _Stdv, -1., p1);
    delete p1;
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

int KrigingCalcul::_needSigma00()
{
  if (!_isPresentMatrix("Sigma00", _Sigma00)) return 1;
  return 0;
}

int KrigingCalcul::_needSigma0()
{
  if (!_isPresentMatrix("Sigma0", _Sigma0)) return 1;
  return 0;
}

int KrigingCalcul::_needZ0p()
{
  if (! _Z0p.empty()) return 0;
  if (!_isPresentVector("Zp", _Zp)) return 1;

  // Sample the active values for collocated information
  _Z0p = VH::sample(*_Zp, *_rankColCok);
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

int KrigingCalcul::_needSigma00p()
{
  if (_ncck <= 0) return 1;
  if (_needSigma00()) return 1;
  _Sigma00p = MatrixRectangular::sample(_Sigma00, *_rankColCok, VectorInt());
  return 0;
}

int KrigingCalcul::_needSigma00pp()
{
  if (_ncck <= 0) return 1;
  if (_needSigma00()) return 1;
  _Sigma00pp = MatrixSquareSymmetric::sample(_Sigma00, *_rankColCok);
  return 0;
}

int KrigingCalcul::_needSigma0p()
{
  if (_ncck <= 0) return 1;
  if (_needSigma0()) return 1;

  if (_Sigma0p == nullptr)
    _Sigma0p = MatrixRectangular::sample(_Sigma0, VectorInt(), *_rankColCok);
  return 0;
}

int KrigingCalcul::_needX0p()
{
  if (_ncck <= 0) return 1;
  if (_needX0()) return 1;

  if (_X0p == nullptr)
    _X0p = MatrixRectangular::sample(_X0, *_rankColCok, VectorInt());
  return 0;
}

int KrigingCalcul::_needBeta()
{
  if (!_Beta.empty()) return 0;
  if (!_isPresentVector("Z", _Z)) return 1;
  if (_needSigmac()) return 1;
  if (_needXtInvSigma()) return 1;

  VectorDouble XtInvCZ = _XtInvSigma->prodMatVec(*_Z);

  if (_flagBayes)
  {
    if (!_isPresentVector("PriorMean", _PriorMean)) return 1;
    if (_needInvPriorCov()) return 1;
    VectorDouble InvSMBayes = _InvPriorCov->prodMatVec(*_PriorMean);
    VH::linearCombinationInPlace(1., XtInvCZ, 1., InvSMBayes, XtInvCZ);
  }

  _Beta = _Sigmac->prodMatVec(XtInvCZ);
  return 0;
}

int KrigingCalcul::_needY0()
{
  if (_needX0()) return 1;
  if (_needInvSigmaSigma0()) return 1;

  if (_Y0 == nullptr)
    _Y0 = new MatrixRectangular(_nrhs, _nbfl);

  MatrixRectangular* LambdaSKtX = new MatrixRectangular(_nrhs, _nbfl);
  LambdaSKtX->prodMatMatInPlace(_InvSigmaSigma0, _X, true, false);
  _Y0->linearCombination(1., _X0, -1., LambdaSKtX);
  delete LambdaSKtX;
  return 0;
}

int KrigingCalcul::_needY0p()
{
  if (_needX0p()) return 1;
  if (_needSigma0p()) return 1;
  if (_needXtInvSigma()) return 1;

  if (_Y0p == nullptr) _Y0p = new MatrixRectangular(_ncck, _nbfl);
  _Y0p->prodMatMatInPlace(_Sigma0p, _XtInvSigma, true, true);
  _Y0p->linearCombination(1., _X0p, -1., _Y0p);
  return 0;
}

int KrigingCalcul::_needMuUK()
{
  if (_flagSK) return 0;
  if (_needSigmac()) return 1;
  if (_needY0()) return 1;
  if (_MuUK == nullptr) _MuUK = new MatrixRectangular(_nbfl, _nrhs);

  if (_ncck > 0)
  {
    if (_needY0p()) return 1;
    if (_needLambda0()) return 1;

    MatrixRectangular* LtY = new MatrixRectangular(_nrhs, _nbfl);
    LtY->prodMatMatInPlace(_Lambda0, _Y0p, true);
    LtY->linearCombination(1., _Y0, -1., LtY);

    _MuUK->prodMatMatInPlace(_Sigmac, LtY, false, true);
    delete LtY;
  }
  else
  {
    _MuUK->prodMatMatInPlace(_Sigmac, _Y0, false, true);
  }

  return 0;
}

int KrigingCalcul::_needInvSigmaSigma0()
{
  if (_needSigma0()) return 1;
  if (_needInvSigma()) return 1;
  if (_InvSigmaSigma0 == nullptr)
    _InvSigmaSigma0 = new MatrixRectangular(_neq, _nrhs);
  _InvSigmaSigma0->prodMatMatInPlace(_InvSigma, _Sigma0);
  return 0;
}

int KrigingCalcul::_needLambdaSK()
{
  if (_ncck > 0)
  {
    if (_needInvSigma()) return 1;
    if (_needSigma0p()) return 1;
    if (_needLambda0()) return 1;

    MatrixRectangular* S = new MatrixRectangular(_neq, _nrhs);
    S->prodMatMatInPlace(_Sigma0p, _Lambda0);
    S->linearCombination(1., _Sigma0, -1., S);
    if (_LambdaSK == nullptr) _LambdaSK = new MatrixRectangular(_neq, _nrhs);
    _LambdaSK->prodMatMatInPlace(_InvSigma, S);
    delete S;
  }
  else
  {
    if (_needInvSigmaSigma0()) return 1;
    _LambdaSK = _InvSigmaSigma0->clone();
  }
  return 0;
}

int KrigingCalcul::_needLambdaUK()
{
  if (_needXtInvSigma()) return 1;
  if (_needLambdaSK()) return 1;
  if (_needMuUK()) return 1;

  MatrixRectangular* p1 = new MatrixRectangular(_neq, _nrhs);
  p1->prodMatMatInPlace(_XtInvSigma, _MuUK, true, false);

  if (_LambdaUK == nullptr) _LambdaUK = new MatrixRectangular(_neq, _nrhs);
  _LambdaUK->linearCombination(1., _LambdaSK, 1., p1);
  delete p1;

  return 0;
}

void KrigingCalcul::_printMatrix(const String& name, const AMatrix* mat)
{
  if (mat == nullptr || mat->empty()) return;
  message(" - %s (%d, %d)\n", name.c_str(), mat->getNRows(), mat->getNCols());
}

void KrigingCalcul::_printVector(const String& name, const VectorDouble* vec)
{
  if (vec == nullptr) return;
  if (vec->size() <= 0) return;
  message(" - %s (%d)\n", name.c_str(), (int)vec->size());
}

void KrigingCalcul::printStatus() const
{
  mestitle(1, "List of arrays used in 'KrigingCalcul'");
  message("\nGeneral Parameters\n");
  message("Number of Covariance Rows = %d\n", _neq);
  message("Number of Drift equations = %d\n", _nbfl);
  message("Number of Right_Hand sides = %d\n", _nrhs);
  if (_ncck > 0)
  {
    message("Number of Collocated Variables = %d\n", _ncck);
    VH::display("Rank of Collocated Variables", *_rankColCok, false);
  }
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
  if (_Z != nullptr) _printVector("Z", _Z);
  _printVector("PriorMean", _PriorMean);
  _printVector("BetaRef", _BetaRef);
  _printVector("Zp", _Zp);

  message("\nInternal Memory (retrievable)\n");
  _printVector("Zstar", &_Zstar);
  _printVector("Beta", &_Beta);
  _printMatrix("LambdaSK", _LambdaSK);
  _printMatrix("LambdaUK", _LambdaUK);
  _printMatrix("MuUK", _MuUK);
  _printMatrix("Stdv", _Stdv);
  _printMatrix("VarZSK", _VarZSK);
  _printMatrix("VarZUK", _VarZUK);

  message("\nInternal Memory (hidden)\n");
  _printMatrix("XtInvSigma", _XtInvSigma);
  _printMatrix("Y0", _Y0);
  _printMatrix("InvSigmaSigma0", _InvSigmaSigma0);
  _printMatrix("InvSigma", _InvSigma);
  _printMatrix("Sigmac", _Sigmac);
  _printMatrix("InvPriorCov", _InvPriorCov);

  if (_ncck > 0)
  {
    message("\nInternal Memory (for Collocated Cokriging only)\n");
    _printMatrix("Sigma00pp", _Sigma00pp);
    _printMatrix("Sigma00p", _Sigma00p);
    _printMatrix("Sigma0p", _Sigma0p);
    _printMatrix("X0p", _X0p);
    _printMatrix("Y0p", _Y0p);
    _printVector("Z0p", &_Z0p);
    _printMatrix("Lambda0", _Lambda0);
  }
}

bool KrigingCalcul::_isPresentMatrix(const String& name, const AMatrix* mat)
{
  if (mat != nullptr) return true;
  messerr(">>> Matrix %s is missing (required)", name.c_str());
  return false;
}

bool KrigingCalcul::_isPresentVector(const String& name,
                                     const VectorDouble* vec)
{
  if (vec != nullptr) return true;
  messerr(">>> Vector %s is missing (required)", name.c_str());
  return false;
}

bool KrigingCalcul::_isPresentIVector(const String& name,
                                      const VectorInt* vec)
{
  if (vec != nullptr) return true;
  messerr(">>> Vector %s is missing (required)", name.c_str());
  return false;
}

int KrigingCalcul::_needLambda0()
{
  if (_ncck <= 0) return 1;
  if (!_isPresentMatrix("Sigma00", _Sigma00)) return 1;
  if (_needInvSigma()) return 1;
  if (_needSigma00p()) return 1;
  if (_needSigma0p()) return 1;
  if (_needSigmac()) return 1;
  if (_needY0p()) return 1;
  if (_needY0()) return 1;
  if (_needSigma00pp()) return 1;

  MatrixRectangular* Sigma0ptInvSigma = new MatrixRectangular(_ncck, _neq);
  Sigma0ptInvSigma->prodMatMatInPlace(_Sigma0p, _InvSigma, true);

  MatrixRectangular* Y0pSigmac = new MatrixRectangular(_ncck, _nbfl);
  Y0pSigmac->prodMatMatInPlace(_Y0p, _Sigmac);

  MatrixSquareSymmetric* bot1 = new MatrixSquareSymmetric(_ncck);
  bot1->prodMatMatInPlace(Sigma0ptInvSigma, _Sigma0p);
  MatrixSquareSymmetric* bot2 = new MatrixSquareSymmetric(_ncck);
  bot2->prodMatMatInPlace(Y0pSigmac, _Y0p, false, true);
  MatrixSquareSymmetric* bot = new MatrixSquareSymmetric(_ncck);
  bot = _Sigma00pp->clone();
  bot->linearCombination(1., bot, -1., bot1, +1., bot2);
  delete bot1;
  delete bot2;

  if (bot->invert())
  {
    delete bot;
    delete Sigma0ptInvSigma;
    delete Y0pSigmac;
    return 1;
  }

  MatrixRectangular* top1 = new MatrixRectangular(_ncck, _nrhs);
  top1->prodMatMatInPlace(Sigma0ptInvSigma, _Sigma0);
  MatrixRectangular* top2 = new MatrixRectangular(_ncck, _nrhs);
  top2->prodMatMatInPlace(Y0pSigmac, _Y0, false, true);
  MatrixRectangular* top = _Sigma00p->clone();
  top->linearCombination(1., top, -1., top1, +1., top2);
  delete top1;
  delete top2;

  _Lambda0 = new MatrixRectangular(_ncck, _nrhs);
  _Lambda0->prodMatMatInPlace(bot, top);

  delete bot;
  delete top;
  delete Sigma0ptInvSigma;
  delete Y0pSigmac;

  return 0;
}