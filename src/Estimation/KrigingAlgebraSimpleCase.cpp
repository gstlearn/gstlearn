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
#include "Estimation/KrigingAlgebraSimpleCase.hpp"
#include "Basic/VectorNumT.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/String.hpp"
#include "Basic/AStringable.hpp"
#include "geoslib_define.h"
#include <cmath>
#include <memory>
#include <omp.h>

KrigingAlgebraSimpleCase::KrigingAlgebraSimpleCase(bool flagDual,
                                                   const VectorVectorInt* sampleRanks,
                                                   const VectorDouble* Z,
                                                   const MatrixSquareSymmetric* Sigma,
                                                   const MatrixDense* X,
                                                   const MatrixSquareSymmetric* Sigma00,
                                                   const VectorDouble& Means,
                                                   int flagchol,
                                                   bool neighUnique)
  : _Z(nullptr)
  , _sampleRanks(nullptr)
  , _X(nullptr)
  , _Sigma(nullptr)
  , _Sigma00(nullptr)
  , _Sigma0(nullptr)
  , _X0(nullptr)
  , _Means()
  , _InvSigma(nullptr)
  , _cholSigma(nullptr)
  , _XtInvSigma(nullptr)
  , _invSigmaX(nullptr)
  , _XtInvSigmaZ(nullptr)
  , _invSigmac(nullptr)
  , _Beta(nullptr)
  , _LambdaSK(nullptr)
  , _bDual(nullptr)
  , _invSigmaXBeta(nullptr)
  , _Zstar()
  , _LambdaSKtX()
  , _LambdaUK()
  , _MuUK()
  , _Stdv()
  , _VarZSK()
  , _VarZUK()
  , _Y0()
  , _LambdaUKtSigma0()
  , _MuUKtX0t()
  , _invSigmaXMuUK()
  , _X0Beta()
  , _nvar(0)
  , _neq(0)
  , _nbfl(0)
  , _nrhs(0)
  , _flagSK(true)
  , _flagDual(flagDual)
  , _neighUnique(neighUnique)
  , _flagCholesky(flagchol == -1 ? (flagDual || _neighUnique) : flagchol)
  , _dualHasChanged(true)
  , _invSigmaHasChanged(true)
  , _XtInvSigmaHasChanged(true)
{
  (void)setData(Z, sampleRanks, Means);
  (void)setLHS(Sigma, X);
  (void)setVariance(Sigma00);

  if (_flagCholesky)
  {
    _cholSigma = std::make_shared<CholeskyDense>(_Sigma.get());
    _invSigmaX = std::make_shared<MatrixDense>();
  }
  else
  {
    _InvSigma   = std::make_shared<MatrixSquareSymmetric>();
    _XtInvSigma = std::make_shared<MatrixDense>();
  }

  _XtInvSigmaZ = std::make_shared<VectorDouble>();
  _invSigmac      = std::make_shared<MatrixSquareSymmetric>();
  _Beta        = std::make_shared<VectorDouble>();
  _LambdaSK    = std::make_shared<MatrixDense>(); // Weights for SK (Dim: _neq * _nrhs)

  // Following elements are defined for Dual programming
  _bDual         = std::make_shared<VectorDouble>(); // Fake Covariance part in Dual (Dim: _neq)
  _invSigmaXBeta = std::make_shared<VectorDouble>();
}

KrigingAlgebraSimpleCase::KrigingAlgebraSimpleCase(KrigingAlgebraSimpleCase& r)
{

  // Quantities which doesn't evolve with the target in unique Neighborhood

  if (r._neighUnique)
    _copyPtrForUniqueNeigh(r);
  else
    _copyContentForMovingNeigh(r);

  // Quantities which always evolve with the target
  _copyOtherContent(r); // Ptrs
  _copyMatsAndVecs(r);  // Matrices and vectors

  _copyFlags(r);           // various flags
  _copyModelQuantities(r); // Quantities related to the model
}

void KrigingAlgebraSimpleCase::_copyModelQuantities(const KrigingAlgebraSimpleCase& r)
{
  // Quantities related to the model
  _Means = r._Means;
  _nvar  = r._nvar;
  _nbfl  = r._nbfl;

  // Unclear Status but copy can be done
  _neq  = r._neq;
  _nrhs = r._nrhs;
}

void KrigingAlgebraSimpleCase::_copyPtrForUniqueNeigh(KrigingAlgebraSimpleCase& r)
{
  _Z             = r._Z;
  _sampleRanks   = r._sampleRanks;
  _X             = r._X;
  _Sigma         = r._Sigma;
  _InvSigma      = r._InvSigma;
  _cholSigma     = r._cholSigma;
  _XtInvSigma    = r._XtInvSigma;
  _invSigmaX     = r._invSigmaX;
  _XtInvSigmaZ   = r._XtInvSigmaZ;
  _invSigmac        = r._invSigmac;
  _Beta          = r._Beta;
  _bDual         = r._bDual;
  _invSigmaXBeta = r._invSigmaXBeta;
}

void KrigingAlgebraSimpleCase::_copyContentForMovingNeigh(const KrigingAlgebraSimpleCase& r)
{
  if (r._Z == nullptr)
    _Z = std::make_shared<VectorDouble>();
  else
    _Z = std::shared_ptr<VectorDouble>(new VectorDouble(*r._Z));

  if (r._sampleRanks == nullptr)
    _sampleRanks = std::make_shared<VectorVectorInt>();
  else
    _sampleRanks = std::shared_ptr<VectorVectorInt>(new VectorVectorInt(*r._sampleRanks));

  if (r._X == nullptr)
    _X = std::make_shared<MatrixDense>();
  else
    _X = std::shared_ptr<MatrixDense>(r._X->clone());

  if (r._Sigma == nullptr)
    _Sigma = std::make_shared<MatrixSquareSymmetric>();
  else
    _Sigma = std::shared_ptr<MatrixSquareSymmetric>(r._Sigma->clone());

  if (r._flagCholesky)
  {
    if (r._cholSigma == nullptr)
      _cholSigma = std::make_shared<CholeskyDense>();
    else
      _cholSigma = std::shared_ptr<CholeskyDense>(new CholeskyDense(r._Sigma.get()));

    if (r._invSigmaX == nullptr)
      _invSigmaX = std::make_shared<MatrixDense>();
    else
      _invSigmaX = std::shared_ptr<MatrixDense>(_invSigmaX->clone());
  }
  else
  {
    if (r._InvSigma == nullptr)
      _InvSigma = std::make_shared<MatrixSquareSymmetric>();
    else
      _InvSigma = std::shared_ptr<MatrixSquareSymmetric>(r._InvSigma->clone());

    if (r._XtInvSigma == nullptr)
      _XtInvSigma = std::make_shared<MatrixDense>();
    else
      _XtInvSigma = std::shared_ptr<MatrixDense>(r._XtInvSigma->clone());
  }

  if (r._XtInvSigmaZ == nullptr)
    _XtInvSigmaZ = std::make_shared<VectorDouble>();
  else
    _XtInvSigmaZ = std::shared_ptr<VectorDouble>(new VectorDouble(*r._XtInvSigmaZ));

  if (r._invSigmac == nullptr)
    _invSigmac = std::make_shared<MatrixSquareSymmetric>();
  else
    _invSigmac = std::shared_ptr<MatrixSquareSymmetric>(r._invSigmac->clone());

  if (r._Beta == nullptr)
    _Beta = std::make_shared<VectorDouble>();
  else
    _Beta = std::shared_ptr<VectorDouble>(new VectorDouble(*r._Beta));

  if (r._bDual == nullptr)
    _bDual = std::make_shared<VectorDouble>();
  else
    _bDual = std::shared_ptr<VectorDouble>(new VectorDouble(*r._bDual));

  if (r._invSigmaXBeta == nullptr)
    _invSigmaXBeta = std::make_shared<VectorDouble>();
  else
    _invSigmaXBeta = std::shared_ptr<VectorDouble>(new VectorDouble(*r._invSigmaXBeta));
}

void KrigingAlgebraSimpleCase::_copyOtherContent(const KrigingAlgebraSimpleCase& r)
{
  if (r._Sigma00 == nullptr)
    _Sigma00 = std::make_shared<MatrixSquareSymmetric>();
  else
    _Sigma00 = std::shared_ptr<MatrixSquareSymmetric>(r._Sigma00->clone());

  if (r._Sigma0 == nullptr)
    _Sigma0 = std::make_shared<MatrixDense>();
  else
    _Sigma0 = std::shared_ptr<MatrixDense>(r._Sigma0->clone());

  if (r._X0 == nullptr)
    _X0 = std::make_shared<MatrixDense>();
  else
    _X0 = std::shared_ptr<MatrixDense>(r._X0->clone());

  if (r._LambdaSK == nullptr)
    _LambdaSK = std::make_shared<MatrixDense>();
  else
    _LambdaSK = std::shared_ptr<MatrixDense>(r._LambdaSK->clone());
}

void KrigingAlgebraSimpleCase::_copyMatsAndVecs(const KrigingAlgebraSimpleCase& r)
{
  // Vector Double or matrices are simply copied
  _Zstar      = r._Zstar;
  _LambdaSKtX = r._LambdaSKtX;
  _LambdaUK   = r._LambdaUK;
  _MuUK       = r._MuUK;
  _Stdv       = r._Stdv;
  _VarZSK     = r._VarZSK;
  _VarZUK     = r._VarZUK;
  _Y0         = r._Y0;
  _LambdaUKtSigma0         = r._LambdaUKtSigma0;
  _MuUKtX0t         = r._MuUKtX0t;
  _invSigmaXMuUK         = r._invSigmaXMuUK;
  _X0Beta         = r._X0Beta;
}

void KrigingAlgebraSimpleCase::_copyFlags(const KrigingAlgebraSimpleCase& r)
{
  _flagSK               = r._flagSK;
  _flagDual             = r._flagDual;
  _neighUnique          = r._neighUnique;
  _flagCholesky         = r._flagCholesky;
  _dualHasChanged       = r._dualHasChanged;
  _invSigmaHasChanged   = r._invSigmaHasChanged;
  _XtInvSigmaHasChanged = r._XtInvSigmaHasChanged;
}

KrigingAlgebraSimpleCase::~KrigingAlgebraSimpleCase()
{

}

void KrigingAlgebraSimpleCase::_resetAll()
{
  _resetLinkedToZ();
  _resetLinkedToLHS();
  _resetLinkedToRHS();
  _resetLinkedToSigma00();
}

void KrigingAlgebraSimpleCase::setDual(bool status)
{
  _resetAll();
  _flagDual = status;
}

void KrigingAlgebraSimpleCase::_resetLinkedToZ()
{
  _deleteDual();
  _deleteXtInvSigmaZ();
}
void KrigingAlgebraSimpleCase::_resetLinkedToLHS()
{
  _resetLinkedToSigma();
  _resetLinkedToX();
}

void KrigingAlgebraSimpleCase::_resetLinkedToRHS()
{
  _resetLinkedToSigma0();
  _resetLinkedToX0();
}

void KrigingAlgebraSimpleCase::_resetLinkedToX()
{
  _deleteXtInvSigma();
}

void KrigingAlgebraSimpleCase::_resetLinkedToX0()
{
  _deleteZstar();
  _deleteMuUK();
}
void KrigingAlgebraSimpleCase::_resetLinkedToSigma()
{
  _deleteInvSigma();
}
void KrigingAlgebraSimpleCase::_resetLinkedToSigma0()
{
  _deleteStdv();
  _deleteLambdaSK();
  _deleteVarZSK();
  _deleteVarZUK();
}
void KrigingAlgebraSimpleCase::_resetLinkedToSigma00()
{
  _deleteStdv();
}
void KrigingAlgebraSimpleCase::_deleteBeta()
{
  _deleteDual();
  if (_Beta == nullptr) return;
  _Beta->clear();
}
void KrigingAlgebraSimpleCase::_deleteInvSigma()
{
  _deleteLambdaSK();
  _deleteXtInvSigma();
  _invSigmaHasChanged = true;
}
void KrigingAlgebraSimpleCase::_deleteLambdaSK()
{

  _deleteLambdaUK();
  _deleteZstar();
  _deleteVarZSK();
  _deleteMuUK();
  if (_LambdaSK == nullptr) return;
  _LambdaSK->clear();
}
void KrigingAlgebraSimpleCase::_deleteLambdaUK()
{
  _deleteVarZUK();
  _deleteStdv();
  _deleteZstar();
  _LambdaUK.clear();
}
void KrigingAlgebraSimpleCase::_deleteMuUK()
{
  _deleteLambdaUK();
  _deleteStdv();
  _MuUK.clear();
}
void KrigingAlgebraSimpleCase::_deleteSigmac()
{
  _deleteBeta();
  _deleteMuUK();
  _deleteDual();
  if (_invSigmac == nullptr) return;
  _invSigmac->clear();
}
void KrigingAlgebraSimpleCase::_deleteZstar()
{
  _Zstar.clear();
}


void KrigingAlgebraSimpleCase::_deleteXtInvSigmaZ()
{

  _deleteLambdaUK();
  _deleteBeta();
  _deleteDual();

  if (_XtInvSigmaZ == nullptr) return;
  _XtInvSigmaZ->clear();
}
void KrigingAlgebraSimpleCase::_deleteXtInvSigma()
{
  _deleteLambdaUK();
  _deleteBeta();

  _XtInvSigmaHasChanged = true;
}
void KrigingAlgebraSimpleCase::_deleteStdv()
{
  _Stdv.clear();
}
void KrigingAlgebraSimpleCase::_deleteVarZSK()
{
  _deleteStdv();

  _VarZSK.clear();
}
void KrigingAlgebraSimpleCase::_deleteVarZUK()
{
  _VarZUK.clear();
}

void KrigingAlgebraSimpleCase::_deleteDual()
{
  _deleteZstar();
  _dualHasChanged = true;
}

/**
 * @brief Method to be used when the data has changed (e.g. Moving Neighborhood)
 */
void KrigingAlgebraSimpleCase::resetNewData()
{
  _neq = 0;
}

/**
 * @brief Modify the Data Values (and Means)
 *
 * @param Z Data flattened vector (possibly multivariate)
 * @param indices Vector Of Vector of sample ranks
 * @param Means  Vector of known Drift coefficients (optional)
 * @return int
 *
 * @note If one element is not provided, its address (if already defined) is
 * @note kept unchanged (even if its contents may have been updated)
 */
int KrigingAlgebraSimpleCase::setData(const VectorDouble* Z,
                                      const VectorVectorInt* indices,
                                      const VectorDouble& Means)
{
  _resetLinkedToZ();

  // Argument Z
  if (Z != nullptr)
  {
    if (!_checkDimensionVD("Z", Z, &_neq)) return 1;
    _Z = std::make_shared<VectorDouble>(*Z);
  }

  // Argument indices
  if (indices != nullptr)
  {
    if (!_checkDimensionVVI("sampleRanks", indices, &_nvar, &_neq)) return 1;
    _sampleRanks = std::make_shared<VectorVectorInt>(*indices);
  }

  // Argument Means

  if (!_checkDimensionVD("Means", &Means, &_nrhs)) return 1;
  _Means = Means;

  return 0;
}

/**
 * @brief Modify the elements linked to the LHS
 *
 * @param Sigma Data-Data Covariance matrix
 * @param X     Data Drift Matrix
 * @return int
 *
 * @note If one element is not provided, its address (if already defined) is
 * @note kept unchanged (even if its contents may have been updated)
 */
int KrigingAlgebraSimpleCase::setLHS(const MatrixSquareSymmetric* Sigma,
                                     const MatrixDense* X)
{
  _resetLinkedToLHS();

  // Argument Sigma
  if (Sigma != nullptr)
  {
    if (!_checkDimensionMatrix("Sigma", Sigma, &_neq, &_neq)) return 1;
    _Sigma = std::make_shared<MatrixSquareSymmetric>(*Sigma);
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
    _X      = std::make_shared<MatrixDense>(*X);
    _flagSK = (_nbfl <= 0);
  }

  return 0;
}

int KrigingAlgebraSimpleCase::setVariance(const MatrixSquareSymmetric* Sigma00)
{
  if (Sigma00 != nullptr)
  {
    if (!_checkDimensionMatrix("Sigma00", Sigma00, &_nrhs, &_nrhs)) return 1;
    _Sigma00 = std::make_shared<MatrixSquareSymmetric>(*Sigma00);
  }
  return 0;
}

int KrigingAlgebraSimpleCase::setRHS(MatrixDense* Sigma0,
                                     MatrixDense* X0)
{
  _resetLinkedToRHS();

  // Argument Sigma0
  if (Sigma0 == nullptr)
  {
    _Sigma0 = nullptr;
  }
  else
  {
    if (!_checkDimensionMatrix("Sigma0", Sigma0, &_neq, &_nrhs)) return 1;
    _Sigma0 = std::make_shared<MatrixDense>(*Sigma0);
  }

  _X0 = std::make_shared<MatrixDense>(*X0);
  return 0;
}

int KrigingAlgebraSimpleCase::updateRHS()
{
  _resetLinkedToRHS();
  return 0;
}
bool KrigingAlgebraSimpleCase::_checkDimensionVD(const String& name,
                                                 const VectorDouble* vec,
                                                 int* sizeRef)
{
  int size = (int)vec->size();
  if (*sizeRef > 0 && size > 0 && size != *sizeRef)
  {
    messerr("Dimension of %s (%d) incorrect: it should be (%d)",
            name.c_str(), size, *sizeRef);
    return false;
  }
  if (size > 0) *sizeRef = size;
  return true;
}

bool KrigingAlgebraSimpleCase::_checkDimensionVI(const String& name,
                                                 const VectorInt* vec,
                                                 int* sizeRef)
{
  int size = (int)vec->size();
  if (*sizeRef > 0 && size != *sizeRef)
  {
    messerr("Dimension of %s (%d) incorrect: it should be (%d)", name.c_str(), size,
            *sizeRef);
    return false;
  }
  if (size > 0) *sizeRef = size;
  return true;
}

bool KrigingAlgebraSimpleCase::_checkDimensionVVI(const String& name,
                                                  const VectorVectorInt* vec,
                                                  int* size1Ref,
                                                  int* size2Ref)
{
  int count = (int)vec->size();
  if (*size1Ref > 0 && count != *size1Ref)
  {
    messerr("First dimension of %s (%d) incorrect: it should be (%d)", name.c_str(),
            count, *size1Ref);
    return false;
  }
  if (count > 0) *size1Ref = count;

  int size = VH::count(*vec);
  if (*size2Ref > 0 && size != *size2Ref)
  {
    messerr("Second dimension of %s (%d) incorrect: it should be (%d)", name.c_str(),
            size, *size2Ref);
    return false;
  }
  if (size > 0) *size2Ref = size;
  return true;
}

bool KrigingAlgebraSimpleCase::_checkDimensionMatrix(const String& name,
                                                     const AMatrix* mat,
                                                     int* nrowsRef,
                                                     int* ncolsRef)
{
  int nrows = mat->getNRows();
  int ncols = mat->getNCols();
  if (*nrowsRef > 0 && nrows != *nrowsRef)
  {
    messerr("Number of Rows of %s (%d) incorrect: it should be (%d)",
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

VectorDouble& KrigingAlgebraSimpleCase::getEstimation()
{
  _needZstar();
  return _Zstar;
}

bool KrigingAlgebraSimpleCase::_forbiddenWhenDual() const
{
  if (!_flagDual) return true;
  messerr("This option is not available as 'Dual' is switched ON");
  return false;
}

VectorDouble KrigingAlgebraSimpleCase::getStdv()
{
  if (!_forbiddenWhenDual()) return VectorDouble();
  if (_needStdv()) return VectorDouble();
  return _Stdv.getDiagonal();
}

const MatrixSquareSymmetric* KrigingAlgebraSimpleCase::getStdvMat()
{
  if (!_forbiddenWhenDual()) return nullptr;
  if (_needStdv()) return nullptr;
  return &_Stdv;
}

double KrigingAlgebraSimpleCase::getVarianceZstar(int i)
{
  if (!_forbiddenWhenDual()) return TEST;
  if (_flagSK)
  {
    if (_needVarZSK()) return TEST;
    return _VarZSK(i, i);
  }
  if (_needVarZUK()) return TEST;
  return _VarZUK(i, i);
}

VectorDouble KrigingAlgebraSimpleCase::getVarianceZstar()
{
  if (!_forbiddenWhenDual()) return VectorDouble();
  if (_flagSK)
  {
    if (_needVarZSK()) return VectorDouble();
    return _VarZSK.getDiagonal();
  }
  if (_needVarZUK()) return VectorDouble();
  return _VarZUK.getDiagonal();
}

const MatrixSquareSymmetric* KrigingAlgebraSimpleCase::getVarianceZstarMat()
{
  if (!_forbiddenWhenDual()) return nullptr;
  if (_flagSK)
  {
    if (_needVarZSK()) return nullptr;
    return &_VarZSK;
  }
  if (_needVarZUK()) return nullptr;
  return &_VarZUK;
}

const MatrixDense* KrigingAlgebraSimpleCase::getLambda()
{
  if (!_forbiddenWhenDual()) return nullptr;
  if (_flagSK)
  {
    if (_needLambdaSK()) return nullptr;
    return _LambdaSK.get();
  }
  if (_needLambdaUK()) return nullptr;
  return &_LambdaUK;
}

const MatrixDense* KrigingAlgebraSimpleCase::getMu()
{
  if (_needMuUK()) return nullptr;
  return &_MuUK;
}

int KrigingAlgebraSimpleCase::_needInvSigma()
{
  if (!_invSigmaHasChanged) return 0;
  if (_notFindSigma()) return 1;
  if (_flagCholesky)
  {
    _cholSigma = std::make_shared<CholeskyDense>(_Sigma.get());
  }
  else
  {
    _InvSigma = std::make_shared<MatrixSquareSymmetric>();
    _InvSigma->resize(_Sigma->getNRows(), _Sigma->getNCols());
    _Sigma->invert2(*_InvSigma);
  }
  _invSigmaHasChanged = false;
  return 0;
}

double KrigingAlgebraSimpleCase::getLTerm()
{
  if (!_flagDual)
  {
    messerr("This Option requires 'Dual' programming");
    return TEST;
  }
  if (_needDual()) return 1;
  if (_notFindZ()) return 1;

  return VH::innerProduct(*_bDual, *_Z);
}

int KrigingAlgebraSimpleCase::prepare()
{
  if (!_neighUnique) return 0;

  if (_flagDual)
  {
    if (_needDual()) return 1;
  }
  else
  {

    if (_nbfl > 0)
    {
      if (_needBeta()) return 1;
    }
    else if (_needInvSigma())
      return 1;
  }
  return 0;
}
int KrigingAlgebraSimpleCase::_computeZstarWithDual()
{
  if (_needDual()) return 1;
  if (_notFindSigma0()) return 1;
  vect vZstar(_Zstar);
  _Sigma0->prodMatVecInPlace(*_bDual, vZstar, true);
  if (_nbfl > 0)
  {
    if (_notFindX0()) return 1;
    if (_needBeta()) return 1;
    _X0Beta.resize(_nbfl);
    if (_notFindX0()) return 1;
    _X0->prodMatVecInPlace(*_Beta, _X0Beta);
    VH::linearCombinationInPlace(1., _Zstar, 1., _X0Beta, _Zstar);
  }
  else
  {
    if (!_Means.empty())
      VH::linearCombinationInPlace(1., _Zstar, 1., _Means, _Zstar);
  }
  return 0;
}

int KrigingAlgebraSimpleCase::_computeZstarSK()
{
  
  if (_needLambdaSK()) return 1;
    _LambdaSK->prodVecMatInPlace(*_Z, _Zstar, false);

    // Adding Mean per Variable
  if (!_Means.empty())
    VH::linearCombinationInPlace(1., _Zstar, 1., _Means, _Zstar);
  return 0;
}

int KrigingAlgebraSimpleCase::_needZstar()
{
  if (!_Zstar.empty()) return 0;

  _Zstar.resize(_nrhs);
  if (_flagDual)
   return _computeZstarWithDual();

  if (_notFindZ()) return 1;
  if (_flagSK)
    return _computeZstarSK();

  if (_needLambdaUK()) return 1;
    _LambdaUK.prodVecMatInPlace(*_Z, _Zstar);
  return 0;
}

int KrigingAlgebraSimpleCase::_needVarZSK()
{
  if (!_VarZSK.empty()) return 0;
  if (_needLambdaSK()) return 1;
  _VarZSK.resize(_nrhs, _nrhs);
  _VarZSK.prodMatMatInPlace(_LambdaSK.get(), _Sigma0.get(), true, false);
  return 0;
}

int KrigingAlgebraSimpleCase::_needVarZUK()
{
  if (!_VarZUK.empty()) return 0;
  if (_needLambdaUK()) return 1;
  _VarZUK.resize(_nrhs, _nrhs);
  _VarZUK.prodNormMatMatInPlace(&_LambdaUK, _Sigma.get(), true);
  return 0;
}

int KrigingAlgebraSimpleCase::_needStdv()
{
  if (!_Stdv.empty()) return 0;
  if (_notFindSigma00()) return 1;
  _Stdv.resize(_nrhs, _nrhs);

  if (_flagSK)
  {
    if (_needVarZSK()) return 1;
    _Stdv.linearCombination(1., _Sigma00.get(), -1., &_VarZSK);
  }
  else
  {
    if (_needLambdaUK()) return 1;
    if (_needMuUK()) return 1;

    _Stdv.resize(_nrhs, _nrhs);
    _LambdaUKtSigma0.resize(_nrhs, _nrhs);
    _LambdaUKtSigma0.prodMatMatInPlace(&_LambdaUK, _Sigma0.get(), true);
    _MuUKtX0t.resize(_nrhs, _nrhs);
    _MuUKtX0t.prodMatMatInPlace(&_MuUK, _X0.get(), true, true);
    _Stdv.linearCombination(1, _Sigma00.get(), -1., &_LambdaUKtSigma0, +1., &_MuUKtX0t);
  }

  // Transform variance into standard deviation

  for (int icol = 0; icol < _nrhs; icol++)
  {
    vect colcur = _Stdv.getViewOnColumnModify(icol);
    
    double* value = &colcur[icol];
    if (*value < -EPSILON10)
    {
      messerr("Negative variance (%g) element %d", *value, icol);
      return 1;
    }
    *value        = sqrt(ABS(*value));
  }
  return 0;
}

bool KrigingAlgebraSimpleCase::_notFindX()
{
  return !_isPresentMatrix("X", _X.get());
}

bool KrigingAlgebraSimpleCase::_notFindX0()
{
  return !_isPresentMatrix("X0", _X0.get());
}

bool KrigingAlgebraSimpleCase::_notFindSigma()
{
  return !_isPresentMatrix("Sigma", _Sigma.get());
}

bool KrigingAlgebraSimpleCase::_notFindSigma00()
{
  return !_isPresentMatrix("Sigma00", _Sigma00.get());
}

bool KrigingAlgebraSimpleCase::_notFindSigma0()
{
  return !_isPresentMatrix("Sigma0", _Sigma0.get());
}

int KrigingAlgebraSimpleCase::_needXtInvSigmaZ()
{
  if (!_XtInvSigmaZ->empty()) return 0;
  if (_needXtInvSigma()) return 1;
  if (_notFindZ()) return 1;
  _XtInvSigmaZ->resize(_nbfl);
  constvect vZ(*_Z);
  vect vres(*_XtInvSigmaZ);
  if (_flagCholesky)
    _invSigmaX->prodMatVecInPlace(vZ, vres, true);
  else
    _XtInvSigma->prodMatVecInPlace(vZ, vres); // TODO in place
  return 0;
}

int KrigingAlgebraSimpleCase::_needXtInvSigma()
{
  if (!_XtInvSigmaHasChanged) return 0;
  if (_notFindX()) return 1;
  if (_needInvSigma()) return 1;
  if (_flagCholesky)
  {
    _invSigmaX->resize(_neq, _nbfl);
    _cholSigma->solveMatInPlace(*_X, *_invSigmaX);
  }
  else
  {
    _XtInvSigma->resize(_nbfl, _neq);
    _XtInvSigma->prodMatMatInPlace(_X.get(), _InvSigma.get(), true, false);
  }
  _XtInvSigmaHasChanged = false;
  return 0;
}

int KrigingAlgebraSimpleCase::_needSigmac()
{
  if (!_invSigmac->empty()) return 0;
  if (_needXtInvSigma()) return 1;
  _Sigmac.resize(_nbfl, _nbfl);
  if (_flagCholesky)
    _Sigmac.prodMatMatInPlace(_invSigmaX.get(), _X.get(), true);
  else
    _Sigmac.prodMatMatInPlace(_XtInvSigma.get(), _X.get());
  // Compute the inverse matrix
  _invSigmac->resize(_nbfl, _nbfl);
  if (_Sigmac.invert2(*_invSigmac)) return 1;
  return 0;
}

int KrigingAlgebraSimpleCase::_needBeta()
{
  if (!_Beta->empty()) return 0;
  if (_needSigmac()) return 1;
  if (_needXtInvSigmaZ()) return 1;
  _Beta->resize(_nbfl);
  _invSigmac->prodMatVecInPlace(*_XtInvSigmaZ, *_Beta);
  return 0;
}

int KrigingAlgebraSimpleCase::_needMuUK()
{
  if (!_MuUK.empty()) return 0;
  if (_flagSK) return 1;
  if (_notFindX0()) return 1;
  if (_needSigmac()) return 1;
  if (_needLambdaSK()) return 1;

  _MuUK.resize(_nbfl, _nrhs);
  _LambdaSKtX.resize(_nrhs, _nbfl);
  _Y0.resize(_nrhs, _nbfl);

  _LambdaSKtX.prodMatMatInPlace(_LambdaSK.get(), _X.get(), true, false);
  _Y0.linearCombination(1., _X0.get(), -1., &_LambdaSKtX);

  _MuUK.prodMatMatInPlace(_invSigmac.get(), &_Y0, false, true);
  return 0;
}

int KrigingAlgebraSimpleCase::_needLambdaSK()
{
  if (!_LambdaSK->empty()) return 0;
  if (_notFindSigma0()) return 1;
  if (_needInvSigma()) return 1;
  _LambdaSK->resize(_neq, _nrhs);
  if (_flagCholesky)
  {
    _cholSigma->solveMatInPlace(*_Sigma0, *_LambdaSK);
  }
  else
    _LambdaSK->prodMatMatInPlace(_InvSigma.get(), _Sigma0.get());
  return 0;
}

bool KrigingAlgebraSimpleCase::_notFindSampleRanks()
{
  return !_isPresentIIVector("SampleRanks", _sampleRanks.get());
}

bool KrigingAlgebraSimpleCase::_notFindZ()
{
  return !_isPresentVector("Z", _Z.get());
}

int KrigingAlgebraSimpleCase::_needLambdaUK()
{
  if (!_LambdaUK.empty()) return 0;
  _LambdaUK.resize(_neq, _nrhs);

  _needXtInvSigma();
  _needLambdaSK();
  _needMuUK();

  _invSigmaXMuUK.resize(_neq, _nrhs);
  if (_flagCholesky)
    _invSigmaXMuUK.prodMatMatInPlace(_invSigmaX.get(), &_MuUK);
  else
    _invSigmaXMuUK.prodMatMatInPlace(_XtInvSigma.get(), &_MuUK, true);
  _LambdaUK.linearCombination(1., _LambdaSK.get(), 1., &_invSigmaXMuUK);

  return 0;
}

int KrigingAlgebraSimpleCase::_needDual()
{
  if (!_dualHasChanged) return 0;
  if (_notFindZ()) return 1;
  if (_needInvSigma()) return 1;

  _bDual->resize(_neq);
  constvect vZ(*_Z);
  vect vB(*_bDual);

  if (_flagCholesky)
    _cholSigma->solve(vZ, vB);
  else
    _InvSigma->prodMatVecInPlace(vZ, vB, true);
  if (_nbfl > 0)
  {
    if (_needBeta()) return 1;
    _invSigmaXBeta->resize(_neq);
    vect vISXD(*_invSigmaXBeta);
    if (_flagCholesky)
      _invSigmaX->prodMatVecInPlace(*_Beta, vISXD);
    else
      _XtInvSigma->prodMatVecInPlace(*_Beta, vISXD, true);
    VH::linearCombinationInPlace(1., *_bDual, -1., *_invSigmaXBeta, *_bDual);
  }
  _dualHasChanged = false;
  return 0;
}

void KrigingAlgebraSimpleCase::_printMatrix(const String& name, const AMatrix* mat)
{
  if (mat == nullptr || mat->empty()) return;
  message(" - %s (%d, %d)\n", name.c_str(), mat->getNRows(), mat->getNCols());
}

void KrigingAlgebraSimpleCase::_printVector(const String& name, const VectorDouble* vec)
{
  if (vec == nullptr) return;
  if (vec->size() <= 0) return;
  message(" - %s (%d)\n", name.c_str(), (int)vec->size());
}

void KrigingAlgebraSimpleCase::printStatus() const
{
  mestitle(1, "List of arrays used in 'KrigingAlgebraSimpleCase'");
  message("\nGeneral Parameters\n");
  message("Number of Covariance Rows ('_neq') = %d\n", _neq);
  message("Number of Drift equations ('_nbfl') = %d\n", _nbfl);
  message("Number of Right_Hand sides ('_nrhs') = %d\n", _nrhs);

  if (_flagSK)
    message("Working with Known Mean(s)\n");
  else
    message("Working with Unknown Mean(s)\n");

  message("\nExternal Pointers\n");
  _printMatrix("Sigma00", _Sigma00.get());
  _printMatrix("Sigma", _Sigma.get());
  _printMatrix("Sigma0", _Sigma0.get());
  _printMatrix("X", _X.get());
  _printMatrix("X0", _X0.get());
  _printVector("Z", _Z.get());
  _printVector("Means", &_Means);

  message("\nInternal Memory (retrievable)\n");
  _printVector("Zstar", &_Zstar);
  _printVector("Beta", _Beta.get());
  _printMatrix("LambdaSK", _LambdaSK.get());
  _printMatrix("LambdaUK", &_LambdaUK);
  _printMatrix("MuUK", &_MuUK);
  _printMatrix("Stdv", &_Stdv);
  _printMatrix("VarZSK", &_VarZSK);
  _printMatrix("VarZUK", &_VarZUK);

  message("\nInternal Memory (hidden)\n");
  _printMatrix("XtInvSigma", _invSigmaX->transpose());
  _printMatrix("Y0", &_Y0);
  if (!_flagCholesky)
    _printMatrix("InvSigma", _InvSigma.get());
  _printMatrix("Sigmac", _invSigmac.get());
}

bool KrigingAlgebraSimpleCase::_isPresentMatrix(const String& name, const AMatrix* mat)
{
  if (mat != nullptr) return true;
  messerr(">>> Matrix %s is missing (required)", name.c_str());
  messerr("    (generated in KrigingAlgebraSimpleCase::_isPresentMatrix)");
  return false;
}

bool KrigingAlgebraSimpleCase::_isPresentVector(const String& name,
                                                const VectorDouble* vec)
{
  if (vec != nullptr) return true;
  messerr(">>> Vector %s is missing (required)", name.c_str());
  messerr("    (generated in KrigingAlgebraSimpleCase::_isPresentVector)");
  return false;
}

bool KrigingAlgebraSimpleCase::_isPresentIVector(const String& name,
                                                 const VectorInt* vec)
{
  if (vec != nullptr) return true;
  messerr(">>> Vector %s is missing (required)", name.c_str());
  messerr("    (generated in KrigingAlgebraSimpleCase::_isIPresentVector)");
  return false;
}

bool KrigingAlgebraSimpleCase::_isPresentIIVector(const String& name, const VectorVectorInt* vec)
{
  if (vec != nullptr) return true;
  messerr(">>> VectorVector %s is missing (required)", name.c_str());
  messerr("    (generated in KrigingAlgebraSimpleCase::_isIIPresentVector)");
  return false;
}

void KrigingAlgebraSimpleCase::dumpLHS(int nbypas) const
{
  int size = _neq;
  if (!_flagSK) size += _nbfl;
  int npass = (size - 1) / nbypas + 1;

  /* General Header */

  mestitle(0, "LHS of Kriging matrix");
  if (_Sigma != nullptr)
    message("Dimension of the Covariance Matrix  = %d\n", _Sigma->getNRows());
  if (_X != nullptr && !_flagSK)
    message("Dimension of the Drift Matrix       = %d\n", _nbfl);

  // LHS matrices
  for (int ipass = 0; ipass < npass; ipass++)
  {
    int ideb = ipass * nbypas;
    int ifin = MIN(size, ideb + nbypas);
    message("\n");

    // Header line
    tab_prints(NULL, "Rank");
    for (int j = ideb; j < ifin; j++) tab_printi(NULL, j + 1);
    message("\n");

    // LHS Matrix
    for (int i = 0; i < size; i++)
    {
      tab_printi(NULL, i + 1);
      if (i < _neq)
      {
        for (int j = ideb; j < ifin; j++)
        {
          if (j < _neq)
            tab_printg(NULL, _Sigma->getValue(i, j, false));
          else
            tab_printg(NULL, _X->getValue(i, j - _neq, false));
        }
        message("\n");
      }
      else
      {
        for (int j = ideb; j < ifin; j++)
        {
          if (j < _neq)
            tab_printg(NULL, _X->getValue(j, i - _neq, false));
          else
            tab_printg(NULL, 0.);
        }
        message("\n");
      }
    }
  }
}

void KrigingAlgebraSimpleCase::dumpRHS() const
{
  int size = _Sigma0->getNRows();
  // Note: X0 is transposed!
  if (_X0 != nullptr) size += _X0->getNCols();

  // Header line
  tab_prints(NULL, "Rank");
  for (int irhs = 0; irhs < _nrhs; irhs++) tab_printi(NULL, irhs + 1);
  message("\n");

  // RHS Matrix
  for (int i = 0; i < size; i++)
  {
    tab_printi(NULL, i + 1);
    if (i < _neq)
    {
      for (int irhs = 0; irhs < _nrhs; irhs++)
        tab_printg(NULL, _Sigma0->getValue(i, irhs, false));
    }
    else
    {
      if (_X0 != nullptr)
        for (int irhs = 0; irhs < _nrhs; irhs++)
          tab_printg(NULL, _X0->getValue(irhs, i - _neq, false));
    }
    message("\n");
  }
}

// This method cannot be const as it may compute _lambda internally upon request
void KrigingAlgebraSimpleCase::dumpWGT()
{
  MatrixDense* lambda;
  if (_flagSK)
  {
    if (_needLambdaSK()) return;
    lambda = _LambdaSK.get();
  }
  else
  {
    if (_needLambdaUK()) return;
    lambda = &_LambdaUK;
  }
  if (_notFindSampleRanks()) return;
  char string[20];

  /* Header Line */

  tab_prints(NULL, "Rank");
  tab_prints(NULL, "Data");
  for (int irhs = 0; irhs < _nrhs; irhs++)
  {
    (void)gslSPrintf(string, "Z%d*", irhs + 1);
    tab_prints(NULL, string);
  }
  message("\n");

  // Matrix lines
  VectorDouble sum(_nrhs);
  int lec = 0;
  for (int ivar = 0; ivar < _nvar; ivar++)
  {
    if (_nvar > 1) message("Using variable Z%-2d\n", ivar + 1);
    int nbyvar = (*_sampleRanks)[ivar].size();
    sum.fill(0.);

    for (int j = 0; j < nbyvar; j++)
    {
      tab_printi(NULL, lec + 1);
      double value = (*_Z)[lec];
      // Correct printout by the mean locally in case of SK
      if (_flagSK && !_Means.empty()) value += _Means[ivar];
      tab_printg(NULL, value);
      for (int irhs = 0; irhs < _nrhs; irhs++)
      {
        value = lambda->getValue(lec, irhs, false);
        tab_printg(NULL, value);
        sum[irhs] += value;
      }
      message("\n");
      lec++;
    }

    // Display sum of weights
    tab_prints(NULL, "Sum of weights", 2, EJustify::LEFT);
    for (int irhs = 0; irhs < _nrhs; irhs++) tab_printg(NULL, sum[irhs]);
    message("\n");
  }
}

void KrigingAlgebraSimpleCase::dumpAux()
{
  if (_notFindSampleRanks()) return;
  char string[20];

  // For Simple Kriging, dump the information on Means
  if (_nbfl <= 0)
  {
    if (!_Means.empty())
    {
      for (int ivar = 0; ivar < _nvar; ivar++)
        message("Mean for Variable Z%d = %lf\n", ivar + 1, _Means[ivar]);
    }
    return;
  }

  if (_needMuUK()) return;
  if (_needBeta()) return;

  // Header Line
  tab_prints(NULL, "Rank");
  for (int irhs = 0; irhs < _nrhs; irhs++)
  {
    (void)gslSPrintf(string, "Mu%d*", irhs + 1);
    tab_prints(NULL, string);
  }
  tab_prints(NULL, "Coeff");
  message("\n");

  for (int ibfl = 0; ibfl < _nbfl; ibfl++)
  {
    tab_printi(NULL, ibfl + 1);
    for (int irhs = 0; irhs < _nrhs; irhs++)
      tab_printg(NULL, _MuUK.getValue(ibfl, irhs, false));
    tab_printg(NULL, _Beta->at(ibfl));
    message("\n");
  }
}
