/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Matrix/MatrixSymmetric.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/AMatrix.hpp"

/**
 * @brief Perform the Algebra for Kriging and CoKriging
 *
 * It requires the definition of:
 * - the vector of Data values *Z* (possibly multivariate and heterotopic)
 * - the Covariance matrix at data points *Sigma*
 * - the Drift matrix at data points *X* (UK if defined, SK otherwise)
 * - the Covariance matrix at target *Sigma00* (only for calculating variance)
 * - the Drift coefficients *Beta* (for SK)
 *
 * Note:
 * When using SK:
 * - the vector *Z* must be centered by the drift beforehand
 * - the vector *beta* corresponds to the vector of Means.
 */
class GSTLEARN_EXPORT KrigingAlgebra {
public:
  KrigingAlgebra(bool flagDual                        = false,
                 const VectorVectorInt* sampleRanks   = nullptr,
                 const VectorDouble* Z                = nullptr,
                 const MatrixSymmetric* Sigma   = nullptr,
                 const MatrixDense* X           = nullptr,
                 const MatrixSymmetric* Sigma00 = nullptr,
                 const VectorDouble* Means            = nullptr);

  void setDual(bool status);
  void resetNewData();
  int setData(const VectorDouble* Z          = nullptr,
              const VectorVectorInt* indices = nullptr,
              const VectorDouble* Means      = nullptr);
  int setLHS(const MatrixSymmetric* Sigma = nullptr,
             const MatrixDense* X         = nullptr);
  int setRHS(const MatrixDense* Sigma0 = nullptr,
             const MatrixDense* X0     = nullptr);
  int setVariance(const MatrixSymmetric* Sigma00 = nullptr);
  int setBayes(const VectorDouble* PriorMean         = nullptr,
               const MatrixSymmetric* PriorCov = nullptr);
  int setXvalidUnique(const VectorInt* rankXvalidEqs  = nullptr,
                      const VectorInt* rankXvalidVars = nullptr);
  int setColCokUnique(const VectorDouble* Zp      = nullptr,
                      const VectorInt* rankColCok = nullptr);

  void printStatus() const;
  void dumpLHS(int nbypas = 5) const;
  void dumpRHS() const;
  void dumpWGT();
  void dumpAux();

  VectorDouble getEstimation();
  VectorDouble getStdv();
  VectorDouble getVarianceZstar();
  VectorDouble getPostMean();
  const MatrixSymmetric* getStdvMat();
  const MatrixSymmetric* getVarianceZstarMat();
  const MatrixSymmetric* getPostCov();
  const MatrixDense*     getLambda();
  const MatrixDense*     getLambda0();
  const MatrixDense*     getMu();
  double getLTerm();
  bool isDual() const { return _flagDual; }

private:
  static bool _checkDimensionMatrix(const String& name, const AMatrix* mat, int* nrowsRef, int* ncolsRef);
  static bool _checkDimensionVD(const String& name, const VectorDouble* vec, int* sizeRef);
  static bool _checkDimensionVI(const String& name, const VectorInt* vec, int* sizeRef);
  static bool _checkDimensionVVI(const String& name, const VectorVectorInt* vec, int* size1Ref, int* size2Ref);

  static bool _isPresentMatrix(const String& name, const AMatrix* mat);
  static bool _isPresentVector(const String& name, const VectorDouble* vec);
  static bool _isPresentIVector(const String& name, const VectorInt* vec);
  static bool _isPresentIIVector(const String& name, const VectorVectorInt* vec);

  void _resetLinkedToSampleRanks();
  void _resetLinkedToZ();
  void _resetLinkedToLHS();
  void _resetLinkedToRHS();
  void _resetLinkedtoVar0();
  void _resetLinkedToBayes();
  void _resetLinkedToColCok();
  void _resetLinkedToXvalid();

  int _needX();
  int _needX0();
  int _needSigma();
  int _needSigma0();
  int _needSigma00();
  int _needBeta();
  int _needInvSigma();
  int _needLambdaSK();
  int _needLambdaUK();
  int _needMuUK();
  int _needSigmac();
  int _needZstar();
  int _needY0();
  int _needXtInvSigma();
  int _needStdv();
  int _needVarZSK();
  int _needVarZUK();
  int _needInvPriorCov();
  int _needSigma0p();
  int _needSigma00p();
  int _needSigma00pp();
  int _needX0p();
  int _needY0p();
  int _needZ0p();
  int _needLambda0();
  int _needInvSigmaSigma0();
  int _needPriorCov();
  int _needPriorMean();
  int _needSampleRanks();
  int _needZ();
  int _needZp();
  int _needColCok();
  int _needXvalid();
  int _needDual();

  int _patchRHSForXvalidUnique();
  int _patchColCokVarianceZstar(MatrixSymmetric* varZK);

  void _deleteX();
  void _deleteX0();
  void _deleteSigma();
  void _deleteSigma0();
  void _deleteSigma00();
  void _deleteBeta();
  void _deleteInvSigma();
  void _deleteLambdaSK();
  void _deleteLambdaUK();
  void _deleteMuUK();
  void _deleteSigmac();
  void _deleteZstar();
  void _deleteY0();
  void _deleteXtInvSigma();
  void _deleteStdv();
  void _deleteVarZSK();
  void _deleteVarZUK();
  void _deleteInvPriorCov();
  void _deleteSigma0p();
  void _deleteSigma00p();
  void _deleteSigma00pp();
  void _deleteX0p();
  void _deleteY0p();
  void _deleteZ0p();
  void _deleteLambda0();
  void _deleteInvSigmaSigma0();
  void _deleteInvSigma00vv();
  void _deletePriorCov();
  void _deletePriorMean();
  void _deleteIndices();
  void _deleteZ();
  void _deleteZp();
  void _deleteColCok();
  void _deleteXvalid();
  void _deleteDual();

  static void _printMatrix(const String& name, const AMatrix* mat);
  static void _printVector(const String& name, const VectorDouble* vec);

  bool _forbiddenWhenDual() const;
  void _resetAll();

private:
  // Following information should not be removed in destructor
  const MatrixSymmetric* _Sigma00;  // Variance at Target (Dim: _nrhs * _nrhs)
  const MatrixSymmetric* _Sigma;    // Covariance Matrix (Dim: _neq * _neq)
  const MatrixDense* _Sigma0;       // Covariance at Target (Dim: _neq * _nrhs)
  const MatrixDense* _X;            // Drift at Data (Dim: _neq * _nbfl)
  const MatrixSymmetric* _PriorCov; // Bayesian Prior Covariance (Dim: _nbfl * _nbfl)
  const VectorDouble* _Z;                 // Data [flattened] (Dim: _neq)
  const MatrixDense* _X0;           // Drift at Target (Dim: _nrhs * _nbfl)
  const VectorDouble* _PriorMean;         // Prior Bayesian Mean (Dim: _nbfl)
  const VectorDouble* _Means;             // Fixed drift coefficients
  const VectorDouble* _Zp;                // Vector of values for collocation
  const VectorInt* _rankColCok;           // Ranks of collocated variables
  const VectorInt* _rankXvalidEqs;        // Ranks of the cross-validated Equations
  const VectorInt* _rankXvalidVars;       // Ranks of the cross-validated Variables
  const VectorVectorInt* _sampleRanks;    // Vector of Vector of sampl indices per variable

  // Following elements can be retrieved by Interface functions
  VectorDouble _Zstar;            // Estimated values (Dim: _nrhs)
  VectorDouble _Beta;             // Drift coefficients (Dim: _nbfl)
  MatrixDense _LambdaSK;          // Weights for SK (Dim: _neq * _nrhs)
  MatrixDense _LambdaUK;          // Weights for UK (Dim: _neq * _nrhs)
  MatrixDense _MuUK;              // Lagrange multipliers (Dim: _nbfl * _nrhs)
  MatrixSymmetric _Stdv;          // Estimation stdv. (Dim: _nrhs * _nrhs)
  MatrixSymmetric _VarZSK;        // Estimator variance in SK (Dim: _nrhs * _nrhs)
  MatrixSymmetric _VarZUK;        // Estimator variance in UK (Dim: _nrhs * _nrhs)

  // Following elements are defined for internal storage
  MatrixDense _XtInvSigma;      // X^t * Inv{Sigma} (Dim: _nbfl * _neq);
  MatrixDense _Y0;              // X0 - LambdaSK * X^t (Dim: _nrhs * _nbfl)
  MatrixDense _InvSigmaSigma0;  // Inv{Sigma} * Sigma0 (Dim: _neq * _nrhs)
  MatrixSymmetric _InvSigma;    // Inv{Sigma} (Dim: _neq * _neq)
  MatrixSymmetric _Sigmac;      // Inv{X^t * Inv{Sigma} * X} (Dim: _nbfl * _nbfl)
  MatrixSymmetric _InvPriorCov; // Inv{PriorCov} (Dim: _nbfl * _nbfl)

  // Following elements are defined for internal storage (collocated case in UN)
  MatrixSymmetric _Sigma00pp;        // ColCok Variance T-T (Dim: _ncck * _ncck)
  MatrixDense _Sigma00p;             // ColCok Variance D-T (Dim: _ncck * _nrhs)
  MatrixDense _Sigma0p;              // Collocated Covariance (Dim: _neq * _ncck)
  MatrixDense _X0p;                  // Collocated Drift (Dim: _ncck * _nbfl)
  MatrixDense _Y0p;                  // X0p - Sigma0p^t * Inv{Sigma} * X (Dim: _ncck *_nbfl)
  VectorInt _rankColVars;            // Ranks of active collocated variables
  VectorDouble _Z0p;                 // Vector of (active) collocated values
  MatrixDense _Lambda0;              // Collocated weights (Dim: _ncck * _nrhs)

  // Following elements are defined for Dual programming
  VectorDouble _bDual; // Fake Covariance part in Dual (Dim: _neq)
  VectorDouble _cDual; // Fake Drift part in Dual (Dim: _nbfl)

  // Following elements are defined for internal storage (Cross-validation in UN)
  MatrixDense _C_RHS; // Fictitious Right-hand side (covariance part)
  MatrixDense _X_RHS; // Fictitious Right-hand side (drift part)

  // Additional parameters
  int _nvar;
  int _neq;
  int _nbfl;
  int _nrhs;
  int _ncck;    // Number of additional samples for ColCok in Unique Neighborhood
  int _nxvalid; // Number of samples in XValid in Unique Neighborhood
  bool _flagSK;
  bool _flagBayes;
  bool _flagDual;
};
