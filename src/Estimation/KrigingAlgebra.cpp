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
#include "Estimation/KrigingAlgebra.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/String.hpp"
#include "Basic/AStringable.hpp"

KrigingAlgebra::KrigingAlgebra(bool flagDual,
                               const VectorVectorInt* sampleRanks,
                               const VectorDouble* Z,
                               const MatrixSymmetric* Sigma,
                               const MatrixDense* X,
                               const MatrixSymmetric* Sigma00,
                               const VectorDouble* Means)
  : _Sigma00(nullptr)
  , _Sigma(nullptr)
  , _Sigma0(nullptr)
  , _X(nullptr)
  , _PriorCov(nullptr)
  , _Z(nullptr)
  , _X0(nullptr)
  , _PriorMean(nullptr)
  , _Means(nullptr)
  , _Zp(nullptr)
  , _rankColCok(nullptr)
  , _rankXvalidEqs(nullptr)
  , _rankXvalidVars(nullptr)
  , _sampleRanks(nullptr)
  , _Zstar()
  , _Beta()
  , _LambdaSK()
  , _LambdaUK()
  , _MuUK()
  , _Stdv()
  , _VarZSK()
  , _VarZUK()
  , _XtInvSigma()
  , _Y0()
  , _InvSigmaSigma0()
  , _InvSigma()
  , _Sigmac()
  , _InvPriorCov()

  , _Sigma00pp()
  , _Sigma00p()
  , _Sigma0p()
  , _X0p()
  , _Y0p()
  , _rankColVars()
  , _Z0p()
  , _Lambda0()
  , _bDual()
  , _cDual()
  , _C_RHS()
  , _X_RHS()

  , _nvar(0)
  , _neq(0)
  , _nbfl(0)
  , _nrhs(0)
  , _ncck(0)
  , _nxvalid(0)
  , _flagSK(true)
  , _flagBayes(false)
  , _flagDual(flagDual) {
  (void)setData(Z, sampleRanks, Means);
  (void)setLHS(Sigma, X);
  (void)setVariance(Sigma00);
}

KrigingAlgebra::~KrigingAlgebra() {
  _resetAll();
}

void KrigingAlgebra::_resetAll() {
  _resetLinkedToZ();
  _resetLinkedToLHS();
  _resetLinkedToRHS();
  _resetLinkedtoVar0();
  _resetLinkedToBayes();
  _resetLinkedToColCok();
  _resetLinkedToXvalid();
}

void KrigingAlgebra::setDual(bool status) {
  _resetAll();
  _flagDual = status;
}
void KrigingAlgebra::_resetLinkedToSampleRanks() {
  _deleteIndices();
}
void KrigingAlgebra::_resetLinkedToZ() {
  _deleteZ();
}
void KrigingAlgebra::_resetLinkedToLHS() {
  _deleteSigma();
  _deleteX();
  _deleteDual();
}
void KrigingAlgebra::_resetLinkedToRHS() {
  _deleteSigma0();
  _deleteX0();
}
void KrigingAlgebra::_resetLinkedtoVar0() {
  _deleteSigma00();
}
void KrigingAlgebra::_resetLinkedToBayes() {
  _deletePriorCov();
  _deletePriorMean();
}
void KrigingAlgebra::_resetLinkedToColCok() {
  _deleteZp();
  _deleteColCok();
}
void KrigingAlgebra::_resetLinkedToXvalid() {
  _deleteXvalid();
}

void KrigingAlgebra::_deleteX() {
  _deleteXtInvSigma();
  _deleteSigmac();
  // Cannot delete _X due to constness
}
void KrigingAlgebra::_deleteX0() {
  _deleteX0p();
  _deleteY0();
  _deleteStdv();
  // Cannot delete _X0 due to constness
}
void KrigingAlgebra::_deleteSigma() {
  _deleteInvSigma();
  // Cannot delete _Sigma due to constness
}
void KrigingAlgebra::_deleteSigma0() {
  _deleteStdv();
  _deleteSigma0p();
  _deleteInvSigmaSigma0();
  _deleteVarZSK();
  _deleteVarZUK();
  // Cannot delete _Sigma0 due to constness
}
void KrigingAlgebra::_deleteSigma00() {
  _deleteSigma00p();
  _deleteSigma00pp();
  _deleteLambda0();
  _deleteStdv();
  // Cannot delete _Sigma00 due to constness
}
void KrigingAlgebra::_deleteBeta() {
  _deleteZstar();
  _Beta.clear();
}
void KrigingAlgebra::_deleteInvSigma() {
  _deleteInvSigmaSigma0();
  _deleteLambdaSK();
  _deleteLambda0();
  _deleteXtInvSigma();

  _InvSigma.clear();
}
void KrigingAlgebra::_deleteLambdaSK() {
  _deleteLambdaUK();
  _deleteZstar();
  _deleteVarZSK();

  _LambdaSK.clear();
}
void KrigingAlgebra::_deleteLambdaUK() {
  _deleteVarZUK();
  _deleteStdv();
  _deleteZstar();

  _LambdaUK.clear();
}
void KrigingAlgebra::_deleteMuUK() {
  _deleteLambdaUK();
  _deleteStdv();

  _MuUK.clear();
}
void KrigingAlgebra::_deleteSigmac() {
  _deleteLambda0();
  _deleteBeta();
  _deleteMuUK();
  _deleteDual();

  _Sigmac.clear();
}
void KrigingAlgebra::_deleteZstar() {
  _Zstar.clear();
}
void KrigingAlgebra::_deleteY0() {
  _deleteZstar();
  _deleteMuUK();
  _deleteLambda0();

  _Y0.clear();
}
void KrigingAlgebra::_deleteXtInvSigma() {
  _deleteY0p();
  _deleteLambdaUK();
  _deleteSigmac();
  _deleteBeta();

  _XtInvSigma.clear();
}
void KrigingAlgebra::_deleteStdv() {
  _Stdv.clear();
}
void KrigingAlgebra::_deleteVarZSK() {
  _deleteStdv();

  _VarZSK.clear();
}
void KrigingAlgebra::_deleteVarZUK() {
  _VarZUK.clear();
}
void KrigingAlgebra::_deleteInvPriorCov() {
  _deleteSigmac();
  _deleteBeta();

  _InvPriorCov.clear();
}
void KrigingAlgebra::_deleteSigma0p() {
  _deleteY0p();
  _deleteLambdaSK();
  _deleteLambda0();
  _deleteVarZUK();

  _Sigma0p.clear();
}
void KrigingAlgebra::_deleteSigma00p() {
  _deleteLambda0();
  _deleteStdv();

  _Sigma00p.clear();
}
void KrigingAlgebra::_deleteSigma00pp() {
  _deleteLambda0();
  _deleteVarZUK();

  _Sigma00pp.clear();
}
void KrigingAlgebra::_deleteX0p() {
  _deleteY0p();

  _X0p.clear();
}
void KrigingAlgebra::_deleteY0p() {
  _deleteZstar();
  _deleteMuUK();
  _deleteLambda0();

  _Y0.clear();
}
void KrigingAlgebra::_deleteZ0p() {
  _deleteZstar();

  _Z0p.clear();
}
void KrigingAlgebra::_deleteLambda0() {
  _deleteVarZUK();
  _deleteMuUK();
  _deleteLambdaSK();

  _Lambda0.clear();
}
void KrigingAlgebra::_deleteInvSigmaSigma0() {
  _deleteY0();
  _deleteLambdaSK();

  _InvSigmaSigma0.clear();
}

void KrigingAlgebra::_deletePriorCov() {
  _deleteInvPriorCov();
  // Cannot delete _PriorCov due to constness
}
void KrigingAlgebra::_deletePriorMean() {
  _deleteBeta();
  // Cannot delete _PriorMean due to constness
}
void KrigingAlgebra::_deleteIndices() {
  _deleteZ();
}
void KrigingAlgebra::_deleteZ() {
  _deleteZstar();
  _deleteBeta();
  // Cannot delete _Z due to constness

  _deleteDual();
}
void KrigingAlgebra::_deleteZp() {
  _deleteZ0p();
  // Cannot delete _Zp due to constness
}
void KrigingAlgebra::_deleteColCok() {
  _deleteX0p();
  _deleteZ0p();
  _deleteSigma0p();
  _deleteSigma00p();
  _deleteSigma00pp();
  // Cannot delete _rankColCok due to constness
}
void KrigingAlgebra::_deleteXvalid() {
  _nxvalid = 0;
  // Cannot delete _rankXvalidEqs or _rankXvalidVars due to constness
  _C_RHS.clear();
  _X_RHS.clear();
}
void KrigingAlgebra::_deleteDual() {
  if (!_flagDual) return;
  _bDual.clear();
  _cDual.clear();
}

/**
 * @brief Method to be used when the data has changed (e.g. Moving Neighborhood)
 */
void KrigingAlgebra::resetNewData() {
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
int KrigingAlgebra::setData(const VectorDouble* Z,
                            const VectorVectorInt* indices,
                            const VectorDouble* Means) {
  _resetLinkedToZ();
  _resetLinkedToSampleRanks();

  // Argument Z
  if (Z != nullptr) {
    if (!_checkDimensionVD("Z", Z, &_neq)) return 1;
    _Z = Z;
  }

  // Argument indices
  if (indices != nullptr) {
    if (!_checkDimensionVVI("sampleRanks", indices, &_nvar, &_neq)) return 1;
    _sampleRanks = indices;
  }

  // Argument Means
  if (Means != nullptr) {
    if (!_checkDimensionVD("Means", Means, &_nrhs)) return 1;
    _Means = Means;
  }
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
int KrigingAlgebra::setLHS(const MatrixSymmetric* Sigma,
                           const MatrixDense* X) {
  _resetLinkedToLHS();

  // Argument Sigma
  if (Sigma != nullptr) {
    if (!_checkDimensionMatrix("Sigma", Sigma, &_neq, &_neq)) return 1;
    _Sigma = Sigma;
  }

  // Argument X
  if (X == nullptr || X->getNRows() <= 0 || X->getNCols() <= 0) {
    _X      = nullptr;
    _flagSK = true;
  } else {
    if (!_checkDimensionMatrix("X", X, &_neq, &_nbfl)) return 1;
    _X      = X;
    _flagSK = (_nbfl <= 0);
  }

  return 0;
}

int KrigingAlgebra::setVariance(const MatrixSymmetric* Sigma00) {
  if (Sigma00 != nullptr) {
    if (!_checkDimensionMatrix("Sigma00", Sigma00, &_nrhs, &_nrhs)) return 1;
    _Sigma00 = Sigma00;
  }
  return 0;
}

int KrigingAlgebra::setRHS(const MatrixDense* Sigma0,
                           const MatrixDense* X0) {
  _resetLinkedToRHS();

  // Argument Sigma0
  if (Sigma0 == nullptr) {
    _Sigma0 = nullptr;
  } else {
    if (!_checkDimensionMatrix("Sigma0", Sigma0, &_neq, &_nrhs)) return 1;
    _Sigma0 = Sigma0;
  }

  // Argument X0
  if (X0 == nullptr || X0->empty()) {
    _X0 = nullptr;
  } else {
    if (!_checkDimensionMatrix("X0", X0, &_nrhs, &_nbfl)) return 1;
    _X0 = X0;
  }
  return 0;
}

bool KrigingAlgebra::_checkDimensionVD(const String& name,
                                       const VectorDouble* vec,
                                       int* sizeRef) {
  int size = (int)vec->size();
  if (*sizeRef > 0 && size > 0 && size != *sizeRef) {
    messerr("Dimension of %s (%d) incorrect: it should be (%d)",
            name.c_str(), size, *sizeRef);
    return false;
  }
  if (size > 0) *sizeRef = size;
  return true;
}

bool KrigingAlgebra::_checkDimensionVI(const String& name,
                                       const VectorInt* vec,
                                       int* sizeRef) {
  int size = (int)vec->size();
  if (*sizeRef > 0 && size != *sizeRef) {
    messerr("Dimension of %s (%d) incorrect: it should be (%d)", name.c_str(), size,
            *sizeRef);
    return false;
  }
  if (size > 0) *sizeRef = size;
  return true;
}

bool KrigingAlgebra::_checkDimensionVVI(const String& name,
                                        const VectorVectorInt* vec,
                                        int* size1Ref,
                                        int* size2Ref) {
  int count = (int)vec->size();
  if (*size1Ref > 0 && count != *size1Ref) {
    messerr("First dimension of %s (%d) incorrect: it should be (%d)", name.c_str(),
            count, *size1Ref);
    return false;
  }
  if (count > 0) *size1Ref = count;

  int size = VH::count(*vec);
  if (*size2Ref > 0 && size != *size2Ref) {
    messerr("Second dimension of %s (%d) incorrect: it should be (%d)", name.c_str(),
            size, *size2Ref);
    return false;
  }
  if (size > 0) *size2Ref = size;
  return true;
}

bool KrigingAlgebra::_checkDimensionMatrix(const String& name,
                                           const AMatrix* mat,
                                           int* nrowsRef,
                                           int* ncolsRef) {
  int nrows = mat->getNRows();
  int ncols = mat->getNCols();
  if (*nrowsRef > 0 && nrows != *nrowsRef) {
    messerr("Number of Rows of %s (%d) incorrect: it should be (%d)",
            name.c_str(), nrows, *nrowsRef);
    return false;
  }
  if (*ncolsRef > 0 && ncols != *ncolsRef) {
    messerr("Number of Columns of %s (%d) incorrect: it should be (%d)",
            name.c_str(), ncols, *ncolsRef);
    return false;
  }
  if (nrows > 0) *nrowsRef = nrows;
  if (ncols > 0) *ncolsRef = ncols;
  return true;
}

/**
 * @brief Define the inforlation for Collocated Option
 *
 * @param Zp Vector of the Collocated variables (see note)
 * @param rankColCok Vector of ranks of Collocated variables (dim: _nvar)
 * @return int Error return code
 *
 * @note Argument 'rankColCok' gives the variable rank in Target File or -1
 * @note The argument 'Zp' must be corrected by the mean of the variables
 * for the use of Collocated Option in Simple Kriging
 */
int KrigingAlgebra::setColCokUnique(const VectorDouble* Zp, const VectorInt* rankColCok) {
  _resetLinkedToColCok();

  if (Zp == nullptr || rankColCok == nullptr) {
    _ncck = 0;
    return 0;
  }
  if (_flagDual) {
    messerr("Colocated Cokriging is incompatible with 'Dual'");
    return 1;
  }

  // Argument Zp
  if (!_checkDimensionVD("Zp", Zp, &_nrhs)) return 1;

  // Argument rankColCok
  if (!_checkDimensionVI("rankColCok", rankColCok, &_nvar)) return 1;

  _ncck = 0;
  _rankColVars.clear();
  for (int var = 0; var < _nvar; var++) {
    if ((*rankColCok)[var] < 0) continue;
    _rankColVars.push_back((*rankColCok)[var]);
    _ncck++;
  }

  _rankColCok = rankColCok;
  _Zp         = Zp;

  return 0;
}

/**
 * @brief Define the elements of the input Db to be cross-validated
 *
 * @param rankXvalidEqs Vector of equation ranks to be cross-validated
 * @param rankXvalidVars Vector of variable ranks to be cross-validated
 * @return int Error return code
 *
 * @remarks The argument 'rankXvalidVars' only serves in assigning the
 * mean of the correct cross-validated variable (SK only). It is optional in OK
 */
int KrigingAlgebra::setXvalidUnique(const VectorInt* rankXvalidEqs,
                                    const VectorInt* rankXvalidVars)
{
  if (rankXvalidEqs == nullptr || rankXvalidVars == nullptr) return 1;
  if (rankXvalidEqs->size() <= 0) return 1;
  _resetLinkedToXvalid();
  _nrhs           = 0;
  _rankXvalidEqs  = rankXvalidEqs;
  _rankXvalidVars = rankXvalidVars;
  _nxvalid        = (int)rankXvalidEqs->size();
  return _patchRHSForXvalidUnique();
}

int KrigingAlgebra::setBayes(const VectorDouble* PriorMean,
                             const MatrixSymmetric* PriorCov) {
  _resetLinkedToBayes();

  if (PriorMean == nullptr || PriorCov == nullptr) {
    _flagBayes = false;
    return 0;
  }
  if (_flagDual) {
    messerr("Bayesian option is incompatible with 'Dual'");
    return 1;
  }

  if (!_checkDimensionVD("PriorMean", PriorMean, &_nbfl)) return 1;
  if (!_checkDimensionMatrix("PriorCov", PriorCov, &_nbfl, &_nbfl)) return 1;
  _PriorMean = PriorMean;
  _PriorCov  = PriorCov;
  _flagBayes = true;

  return 0;
}

VectorDouble KrigingAlgebra::getEstimation() {
  if (_needZstar()) return VectorDouble();
  return _Zstar;
}

VectorDouble KrigingAlgebra::getPostMean() {
  if (_needBeta()) return VectorDouble();
  return _Beta;
}

bool KrigingAlgebra::_forbiddenWhenDual() const {
  if (!_flagDual) return true;
  messerr("This option is not available as 'Dual' is switched ON");
  return false;
}

VectorDouble KrigingAlgebra::getStdv() {
  if (!_forbiddenWhenDual()) return VectorDouble();
  if (_needStdv()) return VectorDouble();
  return _Stdv.getDiagonal();
}

const MatrixSymmetric* KrigingAlgebra::getStdvMat() {
  if (!_forbiddenWhenDual()) return nullptr;
  if (_needStdv()) return nullptr;
  return &_Stdv;
}

VectorDouble KrigingAlgebra::getVarianceZstar() {
  if (!_forbiddenWhenDual()) return VectorDouble();
  if (_flagSK) {
    if (_needVarZSK()) return VectorDouble();
    return _VarZSK.getDiagonal();
  }
  if (_needVarZUK()) return VectorDouble();
  return _VarZUK.getDiagonal();
}

const MatrixSymmetric* KrigingAlgebra::getVarianceZstarMat() {
  if (!_forbiddenWhenDual()) return nullptr;
  if (_flagSK) {
    if (_needVarZSK()) return nullptr;
    return &_VarZSK;
  }
  if (_needVarZUK()) return nullptr;
  return &_VarZUK;
}

const MatrixDense* KrigingAlgebra::getLambda() {
  if (!_forbiddenWhenDual()) return nullptr;
  if (_flagSK) {
    if (_needLambdaSK()) return nullptr;
    return &_LambdaSK;
  }
  if (_needLambdaUK()) return nullptr;
  return &_LambdaUK;
}

const MatrixDense* KrigingAlgebra::getLambda0() {
  if (_needLambda0()) return nullptr;
  return &_Lambda0;
}

const MatrixDense* KrigingAlgebra::getMu() {
  if (_needMuUK()) return nullptr;
  return &_MuUK;
}

const MatrixSymmetric* KrigingAlgebra::getPostCov() {
  // At this stage, the posterior covariance is contained in '_Sigmac'
  if (_needSigmac()) return nullptr;
  return &_Sigmac;
}

int KrigingAlgebra::_needInvSigma() {
  if (!_InvSigma.empty()) return 0;
  if (_needSigma()) return 1;
  _InvSigma = *_Sigma;

  if (_InvSigma.invert()) return 1;
  return 0;
}

int KrigingAlgebra::_needInvPriorCov() {
  if (!_InvPriorCov.empty()) return 0;
  if (_needPriorCov()) return 1;
  _InvPriorCov = *_PriorCov;
  if (_InvPriorCov.invert()) return 1;
  return 0;
}

double KrigingAlgebra::getLTerm() {
  if (!_flagDual) {
    messerr("This Option requires 'Dual' programming");
    return TEST;
  }
  if (_needDual()) return 1;
  if (_needZ()) return 1;

  return VH::innerProduct(_bDual, *_Z);
}

int KrigingAlgebra::_needZstar() {
  if (!_Zstar.empty()) return 0;
  if (_flagDual) 
  {
    // Particular Dual case
    if (_needDual()) return 1;
    if (_needSigma0()) return 1;

    _Zstar = _Sigma0->prodMatVec(_bDual, true);
    if (_nbfl > 0) {
      if (_needX0()) return 1;
      VectorDouble ext = _X0->prodMatVec(_cDual);
      VH::linearCombinationInPlace(1., _Zstar, 1., ext, _Zstar);
    } else {
      if (!_Means->empty())
        VH::linearCombinationInPlace(1., _Zstar, 1., *_Means, _Zstar);
    }
    return 0;
  }
  if (_needZ()) return 1;
  if (_flagSK || _flagBayes) {
    if (_needLambdaSK()) return 1;
    _Zstar = _LambdaSK.prodVecMat(*_Z, false);

    // Adding Mean per Variable
    if (_flagSK && !_Means->empty()) {
      VectorDouble localMeans = *_Means;
      if (_nxvalid > 0) localMeans = VH::sample(*_Means, *_rankXvalidVars);
      VH::linearCombinationInPlace(1., _Zstar, 1., localMeans, _Zstar);
    }

    if (_flagBayes) {
      if (_needY0()) return 1;
      if (_needBeta()) return 1;
      VectorDouble mean = _Y0.prodMatVec(_Beta);
      VH::linearCombinationInPlace(1., _Zstar, 1., mean, _Zstar);
    }
  } else {
    if (_needLambdaUK()) return 1;
    _Zstar = _LambdaUK.prodVecMat(*_Z, false);
  }

  // Collocated case
  if (_ncck > 0) {
    if (_needZ0p()) return 1;
    VectorDouble Zstar0 = _Lambda0.prodMatVec(_Z0p, true);
    VH::linearCombinationInPlace(1., _Zstar, 1., Zstar0, _Zstar);
  }
  return 0;
}

int KrigingAlgebra::_patchColCokVarianceZstar(MatrixSymmetric* varZK) {
  if (_needLambda0()) return 1;
  if (_needSigma0p()) return 1;
  if (_needSigma00pp()) return 1;
  MatrixSymmetric L0tCL0(_nrhs);
  //TODO check dims
  L0tCL0.prodNormMatMatInPlace(&_Lambda0, &_Sigma00pp, true);

  MatrixDense p2(_nrhs, _neq);
  p2.prodMatMatInPlace(&_Lambda0, &_Sigma0p, true, true);
  MatrixSymmetric L0tCLK(_nrhs);

  if (_flagSK) {
    if (_needLambdaSK()) return 1;
    L0tCLK.prodMatMatInPlace(&p2, &_LambdaSK);
  } else {
    if (_needLambdaUK()) return 1;
    L0tCLK.prodMatMatInPlace(&p2, &_LambdaUK);
  }

  varZK->linearCombination(1., varZK, 2., &L0tCLK, 1., &L0tCL0);
  return 0;
}

int KrigingAlgebra::_needVarZSK() {
  if (!_VarZSK.empty()) return 0;
  if (_needSigma0()) return 1;
  if (_needLambdaSK()) return 1;
  _VarZSK.resize(_nrhs, _nrhs);
  _VarZSK.prodMatMatInPlace(&_LambdaSK, _Sigma0, true, false);

  if (_ncck > 0) {
    if (_patchColCokVarianceZstar(&_VarZSK)) return 1;
  }
  return 0;
}

int KrigingAlgebra::_needVarZUK() {
  if (!_VarZUK.empty()) return 0;
  if (_needSigma0()) return 1;
  if (_needLambdaUK()) return 1;
  _VarZUK.resize(_nrhs, _nrhs);
  _VarZUK.prodNormMatMatInPlace(&_LambdaUK, _Sigma, true);

  if (_ncck > 0) {
    if (_patchColCokVarianceZstar(&_VarZUK)) return 1;
  }
  return 0;
}

int KrigingAlgebra::_needStdv() {
  if (!_Stdv.empty()) return 0;
  if (_needSigma00()) return 1;
  _Stdv.resize(_nrhs, _nrhs);

  if (_flagSK) {
    if (_needVarZSK()) return 1;
    _Stdv.linearCombination(1., _Sigma00, -1., &_VarZSK);
  } else {
    if (_needLambdaUK()) return 1;
    if (_needSigma0()) return 1;
    if (_needMuUK()) return 1;
    _Stdv = *_Sigma00;
    MatrixDense p1(_nrhs, _nrhs);
    p1.prodMatMatInPlace(&_LambdaUK, _Sigma0, true);
    MatrixDense p2(_nrhs, _nrhs);
    p2.prodMatMatInPlace(&_MuUK, _X0, true, true);
    _Stdv.linearCombination(1, &_Stdv, -1., &p1, +1., &p2);

    if (_ncck > 0) {
      if (_needSigma00p()) return 1;
      MatrixSymmetric p3(_nrhs);
      p3.prodMatMatInPlace(&_Sigma00p, &_Lambda0, true);
      _Stdv.linearCombination(1., &_Stdv, -1., &p3);
    }
  }

  // Transform variance into standard deviation

  for (int irow = 0; irow < _nrhs; irow++)
    for (int icol = 0; icol < _nrhs; icol++)
      if (irow >= icol) {
        double value = _Stdv.getValue(irow, icol);
        value        = (value > 0.) ? sqrt(value) : 0.;
        _Stdv.setValue(irow, icol, value);
      }
  return 0;
}

int KrigingAlgebra::_needX() {
  if (!_isPresentMatrix("X", _X)) return 1;
  return 0;
}

int KrigingAlgebra::_needX0() {
  if (!_isPresentMatrix("X0", _X0)) return 1;
  return 0;
}

int KrigingAlgebra::_needSigma() {
  if (!_isPresentMatrix("Sigma", _Sigma)) return 1;
  return 0;
}

int KrigingAlgebra::_needSigma00() {
  if (!_isPresentMatrix("Sigma00", _Sigma00)) return 1;
  return 0;
}

int KrigingAlgebra::_needSigma0() {
  if (!_isPresentMatrix("Sigma0", _Sigma0)) return 1;
  return 0;
}

int KrigingAlgebra::_needZ0p() {
  if (!_Z0p.empty()) return 0;
  if (_needZp()) return 1;
  if (_needColCok()) return 1;

  // Sample the active values for collocated information
  _Z0p = VH::sample(*_Zp, _rankColVars);
  return 0;
}

int KrigingAlgebra::_needXtInvSigma() {
  if (!_XtInvSigma.empty()) return 0;
  if (_needX()) return 1;
  if (_needInvSigma()) return 1;

  _XtInvSigma.resize(_nbfl, _neq);
  _XtInvSigma.prodMatMatInPlace(_X, &_InvSigma, true, false);
  return 0;
}

int KrigingAlgebra::_needSigmac() {
  if (!_Sigmac.empty()) return 0;
  if (_needX()) return 1;
  if (_needXtInvSigma()) return 1;

  _Sigmac.resize(_nbfl, _nbfl);
  _Sigmac.prodMatMatInPlace(&_XtInvSigma, _X);

  // Bayesian case
  if (_flagBayes) {
    if (_needInvPriorCov()) return 1;
    _Sigmac.linearCombination(1., &_Sigmac, 1., &_InvPriorCov);
  }

  // Compute the inverse matrix
  if (_Sigmac.invert()) return 1;
  return 0;
}

int KrigingAlgebra::_needSigma00p() {
  if (!_Sigma00p.empty()) return 0;
  if (_needSigma00()) return 1;
  if (_needColCok()) return 1;
  MatrixDense::sample(_Sigma00p, *_Sigma00, _rankColVars, VectorInt());
  return 0;
}

int KrigingAlgebra::_needSigma00pp() {
  if (!_Sigma00pp.empty()) return 0;
  if (_needSigma00()) return 1;
  if (_needColCok()) return 1;
  MatrixSymmetric::sample(_Sigma00pp, *_Sigma00, _rankColVars);
  return 0;
}

int KrigingAlgebra::_needSigma0p() {
  if (!_Sigma0p.empty()) return 0;
  if (_needSigma0()) return 1;
  if (_needColCok()) return 1;

  MatrixDense::sample(_Sigma0p, *_Sigma0, VectorInt(), _rankColVars);
  return 0;
}

int KrigingAlgebra::_needX0p() {
  if (!_X0p.empty()) return 0;
  if (_needX0()) return 1;
  if (_needColCok()) return 1;

  MatrixDense::sample(_X0p, *_X0, _rankColVars, VectorInt());
  return 0;
}

int KrigingAlgebra::_needBeta() {
  if (!_Beta.empty()) return 0;
  if (_needZ()) return 1;
  if (_needSigmac()) return 1;
  if (_needXtInvSigma()) return 1;

  VectorDouble XtInvCZ = _XtInvSigma.prodMatVec(*_Z);

  if (_flagBayes) {
    if (_needPriorMean()) return 1;
    if (_needInvPriorCov()) return 1;
    VectorDouble InvSMBayes = _InvPriorCov.prodMatVec(*_PriorMean);
    VH::linearCombinationInPlace(1., XtInvCZ, 1., InvSMBayes, XtInvCZ);
  }

  _Beta = _Sigmac.prodMatVec(XtInvCZ);
  return 0;
}

int KrigingAlgebra::_needY0() {
  if (!_Y0.empty()) return 0;
  if (_needX0()) return 1;
  if (_needInvSigmaSigma0()) return 1;

  MatrixDense LambdaSKtX(_nrhs, _nbfl);
  LambdaSKtX.prodMatMatInPlace(&_InvSigmaSigma0, _X, true, false);

  _Y0.resize(_nrhs, _nbfl);
  _Y0.linearCombination(1., _X0, -1., &LambdaSKtX);
  return 0;
}

int KrigingAlgebra::_needY0p() {
  if (!_Y0p.empty()) return 0;
  if (_needX0p()) return 1;
  if (_needSigma0p()) return 1;
  if (_needXtInvSigma()) return 1;

  _Y0p.resize(_ncck, _nbfl);
  _Y0p.prodMatMatInPlace(&_Sigma0p, &_XtInvSigma, true, true);
  _Y0p.linearCombination(1., &_X0p, -1., &_Y0p);
  return 0;
}

int KrigingAlgebra::_needMuUK() {
  if (!_MuUK.empty()) return 0;
  if (_flagSK) return 1;
  if (_needSigmac()) return 1;
  if (_needY0()) return 1;

  _MuUK.resize(_nbfl, _nrhs);

  if (_ncck > 0) {
    if (_needY0p()) return 1;
    if (_needLambda0()) return 1;

    MatrixDense LtY(_nrhs, _nbfl);
    LtY.prodMatMatInPlace(&_Lambda0, &_Y0p, true);
    LtY.linearCombination(1., &_Y0, -1., &LtY);

    _MuUK.prodMatMatInPlace(&_Sigmac, &LtY, false, true);
  } else {
    _MuUK.prodMatMatInPlace(&_Sigmac, &_Y0, false, true);
  }
  return 0;
}

int KrigingAlgebra::_needInvSigmaSigma0() {
  if (!_InvSigmaSigma0.empty()) return 0;
  if (_needSigma0()) return 1;
  if (_needInvSigma()) return 1;
  _InvSigmaSigma0.resize(_neq, _nrhs);
  _InvSigmaSigma0.prodMatMatInPlace(&_InvSigma, _Sigma0);
  return 0;
}

int KrigingAlgebra::_patchRHSForXvalidUnique() {
  _resetLinkedToRHS();
  _resetLinkedtoVar0();

  if (_needInvSigma()) return 1;
  if (_needSigma()) return 1;
  if (_needSigma00()) return 1;
  if (_needXvalid()) return 1;

  // Extract S00
  MatrixSymmetric S00;
  MatrixSymmetric::sample(S00, *_Sigma, *_rankXvalidEqs);

  // Extract alpha and invert it
  MatrixSymmetric alpha;
  MatrixSymmetric::sample(alpha, _InvSigma, *_rankXvalidEqs);
  MatrixSymmetric InvAlpha = alpha;
  InvAlpha.invert();

  // Calculate a1 term
  MatrixSymmetric omega(_nxvalid);
  omega.linearCombination(1., &S00, -1., &InvAlpha);

  if (_nbfl > 0) {
    // Extract beta
    MatrixDense beta;
    MatrixDense::sample(beta, _InvSigma, *_rankXvalidEqs, *_rankXvalidEqs, false, true);

    // Extracting delta
    MatrixSymmetric delta;
    MatrixSymmetric::sample(delta, _InvSigma, *_rankXvalidEqs, true);

    // Extract Drift matrix at target point
    MatrixDense X0;
    MatrixDense::sample(X0, *_X, *_rankXvalidEqs, VectorInt());

    // Extract Drift matrix at data point
    MatrixDense X;
    MatrixDense::sample(X, *_X, *_rankXvalidEqs, VectorInt(), true);

    // Compute epsilon (up to its sign); inv(alpha) * beta
    AMatrix* p1 = MatrixFactory::prodMatMat(&InvAlpha, &beta);
    MatrixDense epsilon(_nxvalid, _nbfl);
    epsilon.prodMatMatInPlace(p1, &X);
    delete p1;

    // Compute a3 (transpose)
    MatrixDense a3(_nxvalid, _nbfl);
    a3.linearCombination(1., &X0, 1., &epsilon);

    // Compute a2 (inverted)
    MatrixSymmetric a2(_nbfl);
    a2.prodNormMatMatInPlace(&X, &delta, true);
    MatrixSymmetric p3(_nbfl);
    p3.prodNormMatMatInPlace(&epsilon, &alpha, true);
    a2.linearCombination(1., &a2, -1., &p3);
    a2.invert();

    // Compute omega
    MatrixSymmetric p4(_nxvalid);
    p4.prodNormMatMatInPlace(&a3, &a2);
    omega.linearCombination(1., &omega, -1., &p4);

    // Patch the Right-hand side vector (Drift part)
    MatrixDense::sample(_X_RHS, *_X, *_rankXvalidEqs, VectorInt());
  }

  // Patch the Right-hand side vector (Covariance part)
  MatrixDense::sample(_C_RHS, *_Sigma, VectorInt(), *_rankXvalidEqs, false, false);
  _C_RHS.unsample(&omega, *_rankXvalidEqs, VectorInt());

  setRHS(&_C_RHS, &_X_RHS);

  setVariance(S00.clone());

  return 0;
}

int KrigingAlgebra::_needPriorCov() {
  if (!_isPresentMatrix("PriorCov", _PriorCov)) return 1;
  return 0;
}

int KrigingAlgebra::_needSampleRanks() {
  if (!_isPresentIIVector("SampleRanks", _sampleRanks)) return 1;
  return 0;
}

int KrigingAlgebra::_needZ() {
  if (!_isPresentVector("Z", _Z)) return 1;
  return 0;
}

int KrigingAlgebra::_needZp() {
  if (!_isPresentVector("Zp", _Zp)) return 1;
  return 0;
}

int KrigingAlgebra::_needColCok() {
  if (!_isPresentIVector("rankColVars", &_rankColVars)) return 1;
  return 0;
}

int KrigingAlgebra::_needXvalid() {
  if (!_isPresentIVector("rankXvalidEqs", _rankXvalidEqs)) return 1;
  return 0;
}

int KrigingAlgebra::_needPriorMean() {
  if (!_isPresentVector("PriorMean", _PriorMean)) return 1;
  return 0;
}

int KrigingAlgebra::_needLambdaSK() {
  if (!_LambdaSK.empty()) return 0;
  // In the case of Xvalidation when drift is present, cannot return the vector of SK weights

  if (_ncck > 0) {
    if (_needInvSigma()) return 1;
    if (_needSigma0p()) return 1;
    if (_needLambda0()) return 1;

    MatrixDense S(_neq, _nrhs);
    S.prodMatMatInPlace(&_Sigma0p, &_Lambda0);
    S.linearCombination(1., _Sigma0, -1., &S);
    _LambdaSK.resize(_neq, _nrhs);
    _LambdaSK.prodMatMatInPlace(&_InvSigma, &S);
  } else {
    if (_needInvSigmaSigma0()) return 1;
    _LambdaSK = _InvSigmaSigma0;
  }
  return 0;
}

int KrigingAlgebra::_needLambdaUK() {
  if (!_LambdaUK.empty()) return 0;
  _LambdaUK.resize(_neq, _nrhs);

  if (_needXtInvSigma()) return 1;
  if (_needLambdaSK()) return 1;
  if (_needMuUK()) return 1;

  MatrixDense p1(_neq, _nrhs);
  p1.prodMatMatInPlace(&_XtInvSigma, &_MuUK, true, false);
  //MatrixRectangular::sum(_LambdaSK, &p1, _LambdaUK);
  _LambdaUK.linearCombination(1., &_LambdaSK, 1., &p1);

  return 0;
}

int KrigingAlgebra::_needDual() {
  if (!_flagDual) return 1;
  if (_needZ()) return 1;
  if (_needInvSigma()) return 1;

  _bDual = _InvSigma.prodMatVec(*_Z, true);
  if (_nbfl > 0) {
    if (_needSigmac()) return 1;
    if (_needXtInvSigma()) return 1;

    VectorDouble wp = _XtInvSigma.prodMatVec(*_Z, false);
    _cDual          = _Sigmac.prodMatVec(wp, false);

    VectorDouble p1 = _XtInvSigma.prodMatVec(_cDual, true);
    VH::linearCombinationInPlace(1., _bDual, -1., p1, _bDual);
  }
  return 0;
}

void KrigingAlgebra::_printMatrix(const String& name, const AMatrix* mat) {
  if (mat == nullptr || mat->empty()) return;
  message(" - %s (%d, %d)\n", name.c_str(), mat->getNRows(), mat->getNCols());
}

void KrigingAlgebra::_printVector(const String& name, const VectorDouble* vec) {
  if (vec == nullptr) return;
  if (vec->size() <= 0) return;
  message(" - %s (%d)\n", name.c_str(), (int)vec->size());
}

void KrigingAlgebra::printStatus() const {
  mestitle(1, "List of arrays used in 'KrigingAlgebra'");
  message("\nGeneral Parameters\n");
  message("Number of Covariance Rows ('_neq') = %d\n", _neq);
  message("Number of Drift equations ('_nbfl') = %d\n", _nbfl);
  message("Number of Right_Hand sides ('_nrhs') = %d\n", _nrhs);
  if (_ncck > 0) {
    message("Number of Collocated Variables ('_ncck') = %d\n", _ncck);
    VH::dump("Rank of Collocated Variables", _rankColVars, false);
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
  _printVector("Z", _Z);
  _printVector("PriorMean", _PriorMean);
  _printVector("Means", _Means);
  _printVector("Zp", _Zp);

  message("\nInternal Memory (retrievable)\n");
  _printVector("Zstar", &_Zstar);
  _printVector("Beta", &_Beta);
  _printMatrix("LambdaSK", &_LambdaSK);
  _printMatrix("LambdaUK", &_LambdaUK);
  _printMatrix("MuUK", &_MuUK);
  _printMatrix("Stdv", &_Stdv);
  _printMatrix("VarZSK", &_VarZSK);
  _printMatrix("VarZUK", &_VarZUK);

  message("\nInternal Memory (hidden)\n");
  _printMatrix("XtInvSigma", &_XtInvSigma);
  _printMatrix("Y0", &_Y0);
  _printMatrix("InvSigmaSigma0", &_InvSigmaSigma0);
  _printMatrix("InvSigma", &_InvSigma);
  _printMatrix("Sigmac", &_Sigmac);
  _printMatrix("InvPriorCov", &_InvPriorCov);

  if (_ncck > 0) {
    message("\nInternal Memory (for Collocated Cokriging only)\n");
    _printMatrix("Sigma00pp", &_Sigma00pp);
    _printMatrix("Sigma00p", &_Sigma00p);
    _printMatrix("Sigma0p", &_Sigma0p);
    _printMatrix("X0p", &_X0p);
    _printMatrix("Y0p", &_Y0p);
    _printVector("Z0p", &_Z0p);
    _printMatrix("Lambda0", &_Lambda0);
  }
}

bool KrigingAlgebra::_isPresentMatrix(const String& name, const AMatrix* mat) {
  if (mat != nullptr) return true;
  messerr(">>> Matrix %s is missing (required)", name.c_str());
  messerr("    (generated in KrigingAlgebra::_isPresentMatrix)");
  return false;
}

bool KrigingAlgebra::_isPresentVector(const String& name,
                                      const VectorDouble* vec) {
  if (vec != nullptr) return true;
  messerr(">>> Vector %s is missing (required)", name.c_str());
  messerr("    (generated in KrigingAlgebra::_isPresentVector)");
  return false;
}

bool KrigingAlgebra::_isPresentIVector(const String& name,
                                       const VectorInt* vec) {
  if (vec != nullptr) return true;
  messerr(">>> Vector %s is missing (required)", name.c_str());
  messerr("    (generated in KrigingAlgebra::_isIPresentVector)");
  return false;
}

bool KrigingAlgebra::_isPresentIIVector(const String& name, const VectorVectorInt* vec) {
  if (vec != nullptr) return true;
  messerr(">>> VectorVector %s is missing (required)", name.c_str());
  messerr("    (generated in KrigingAlgebra::_isIIPresentVector)");
  return false;
}

int KrigingAlgebra::_needLambda0() {
  if (!_Lambda0.empty()) return 0;

  if (_ncck <= 0) return 1;
  if (_needSigma00()) return 1;
  if (_needSigma0p()) return 1;
  if (_needSigma00p()) return 1;
  if (_needSigma00pp()) return 1;
  if (_needInvSigma()) return 1;
  if (_nbfl > 0) {
    if (_needSigmac()) return 1;
    if (_needY0p()) return 1;
    if (_needY0()) return 1;
  }

  MatrixDense Sigma0ptInvSigma(_ncck, _neq);
  Sigma0ptInvSigma.prodMatMatInPlace(&_Sigma0p, &_InvSigma, true);

  // Determine the Bottom part of the ratio
  MatrixSymmetric bot = _Sigma00pp;

  MatrixSymmetric bot1(_ncck);
  bot1.prodMatMatInPlace(&Sigma0ptInvSigma, &_Sigma0p);

  MatrixDense Y0pSigmac;
  if (_nbfl > 0) {
    Y0pSigmac = MatrixDense(_ncck, _nbfl);
    Y0pSigmac.prodMatMatInPlace(&_Y0p, &_Sigmac);
  }

  MatrixSymmetric bot2;
  if (_nbfl > 0) {
    bot2 = MatrixSymmetric(_ncck);
    bot2.prodMatMatInPlace(&Y0pSigmac, &_Y0p, false, true);
  }
  bot.linearCombination(1., &bot, -1., &bot1, +1., (_nbfl > 0) ? &bot2 : nullptr);

  if (bot.invert()) return 1;

  // Determine the Top part of the ratio
  MatrixDense top = _Sigma00p;

  MatrixDense top1(_ncck, _nrhs);
  top1.prodMatMatInPlace(&Sigma0ptInvSigma, _Sigma0);

  MatrixDense top2;
  if (_nbfl > 0) {
    top2 = MatrixDense(_ncck, _nrhs);
    top2.prodMatMatInPlace(&Y0pSigmac, &_Y0, false, true);
  }
  top.linearCombination(1., &top, -1., &top1, +1., (_nbfl > 0) ? &top2 : nullptr);

  _Lambda0.resize(_ncck, _nrhs);
  _Lambda0.prodMatMatInPlace(&bot, &top);

  return 0;
}

void KrigingAlgebra::dumpLHS(int nbypas) const {
  int size = _neq;
  if (!_flagSK && !_flagBayes) size += _nbfl;
  int npass = (size - 1) / nbypas + 1;

  /* General Header */

  mestitle(0, "LHS of Kriging matrix");
  if (_Sigma != nullptr)
    message("Dimension of the Covariance Matrix  = %d\n", _Sigma->getNRows());
  if (_X != nullptr && !_flagSK && !_flagBayes)
    message("Dimension of the Drift Matrix       = %d\n", _nbfl);

  // LHS matrices
  for (int ipass = 0; ipass < npass; ipass++) {
    int ideb = ipass * nbypas;
    int ifin = MIN(size, ideb + nbypas);
    message("\n");

    // Header line
    tab_prints(NULL, "Rank");
    for (int j = ideb; j < ifin; j++) tab_printi(NULL, j + 1);
    message("\n");

    // LHS Matrix
    for (int i = 0; i < size; i++) {
      tab_printi(NULL, i + 1);
      if (i < _neq) {
        for (int j = ideb; j < ifin; j++) {
          if (j < _neq)
            tab_printg(NULL, _Sigma->getValue(i, j, false));
          else
            tab_printg(NULL, _X->getValue(i, j - _neq, false));
        }
        message("\n");
      } else {
        for (int j = ideb; j < ifin; j++) {
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

void KrigingAlgebra::dumpRHS() const {
  int size = _Sigma0->getNRows();
  // Note: X0 is transposed!
  if (_X0 != nullptr) size += _X0->getNCols();

  // Header line
  tab_prints(NULL, "Rank");
  for (int irhs = 0; irhs < _nrhs; irhs++) tab_printi(NULL, irhs + 1);
  message("\n");

  // RHS Matrix
  for (int i = 0; i < size; i++) {
    tab_printi(NULL, i + 1);
    if (i < _neq) {
      for (int irhs = 0; irhs < _nrhs; irhs++)
        tab_printg(NULL, _Sigma0->getValue(i, irhs, false));
    } else {
      if (_X0 != nullptr)
        for (int irhs = 0; irhs < _nrhs; irhs++)
          tab_printg(NULL, _X0->getValue(irhs, i - _neq, false));
    }
    message("\n");
  }
}

// This method cannot be const as it may compute _lambda internally upon request
void KrigingAlgebra::dumpWGT() {
  const MatrixDense* lambda;
  if (_flagSK || _flagBayes) {
    if (_needLambdaSK()) return;
    lambda = &_LambdaSK;
  } else {
    if (_needLambdaUK()) return;
    lambda = &_LambdaUK;
  }
  if (_needSampleRanks()) return;
  char string[20];

  /* Header Line */

  tab_prints(NULL, "Rank");
  tab_prints(NULL, "Data");
  for (int irhs = 0; irhs < _nrhs; irhs++) {
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
    int nbyvar = (int) (*_sampleRanks)[ivar].size();
    sum.fill(0.);

    for (int j = 0; j < nbyvar; j++)
    {
      tab_printi(NULL, lec + 1);
      double value = (*_Z)[lec];
      // Correct printout by the mean locally in case of SK
      if (_flagSK && !_Means->empty()) value += (*_Means)[ivar];
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

void KrigingAlgebra::dumpAux() {
  if (_needSampleRanks()) return;
  char string[20];

  // For Simple Kriging, dump the information on Means
  if (_nbfl <= 0) {
    if (!_Means->empty()) {
      for (int ivar = 0; ivar < _nvar; ivar++)
        message("Mean for Variable Z%d = %lf\n", ivar + 1, (*_Means)[ivar]);
    }
    return;
  }

  // In Bayesian case, dump the Prior and Posterior information
  if (_flagBayes) {
    VH::dump("Prior Mean", *_PriorMean, false);
    message("Prior Covariance Matrix\n");
    _PriorCov->display();

    VectorDouble postmean = getPostMean();
    VH::dump("Posterior Mean", postmean, false);
    message("Posterior Covariance Matrix\n");
    const MatrixSymmetric* postcov = getPostCov();
    postcov->display();
    return;
  }

  if (_needMuUK()) return;
  if (_needBeta()) return;

  // Header Line
  tab_prints(NULL, "Rank");
  for (int irhs = 0; irhs < _nrhs; irhs++) {
    (void)gslSPrintf(string, "Mu%d*", irhs + 1);
    tab_prints(NULL, string);
  }
  tab_prints(NULL, "Coeff");
  message("\n");

  for (int ibfl = 0; ibfl < _nbfl; ibfl++) {
    tab_printi(NULL, ibfl + 1);
    for (int irhs = 0; irhs < _nrhs; irhs++)
      tab_printg(NULL, _MuUK.getValue(ibfl, irhs, false));
    tab_printg(NULL, _Beta[ibfl]);
    message("\n");
  }
}
