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
#pragma once

#include "gstlearn_export.hpp"

#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/MatrixRectangular.hpp"
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
class GSTLEARN_EXPORT KrigingCalcul
{
public:
  KrigingCalcul(bool flagDual                        = false,
                const VectorDouble* Z                = nullptr,
                const MatrixSquareSymmetric* Sigma   = nullptr,
                const MatrixRectangular* X           = nullptr,
                const MatrixSquareSymmetric* Sigma00 = nullptr,
                const VectorDouble* Means            = nullptr);
  KrigingCalcul(const KrigingCalcul& r)            = delete;
  KrigingCalcul& operator=(const KrigingCalcul& r) = delete;
  virtual ~KrigingCalcul();

  int setData(const VectorDouble* Z     = nullptr,
              const VectorDouble* Means = nullptr);
  int setLHS(const MatrixSquareSymmetric* Sigma = nullptr,
             const MatrixRectangular* X         = nullptr);
  int setRHS(const MatrixRectangular* Sigma0 = nullptr,
             const MatrixRectangular* X0     = nullptr);
  int setVar(const MatrixSquareSymmetric* Sigma00 = nullptr);
  int setColCokUnique(const VectorDouble* Zp      = nullptr,
                      const VectorInt* rankColCok = nullptr);
  int setBayes(const VectorDouble* PriorMean         = nullptr,
               const MatrixSquareSymmetric* PriorCov = nullptr);
  int setXvalidUnique(const VectorInt* rankXvalidEqs  = nullptr,
                      const VectorInt* rankXvalidVars = nullptr);

  void printStatus() const;

  VectorDouble getEstimation();
  VectorDouble getStdv();
  VectorDouble getVarianceZstar();
  VectorDouble getPostMean();
  const MatrixSquareSymmetric* getStdvMat();
  const MatrixSquareSymmetric* getVarianceZstarMat();
  const MatrixSquareSymmetric* getPostCov();
  const MatrixRectangular* getLambda();
  const MatrixRectangular* getLambda0();
  const MatrixRectangular* getMu();

  // Some debugging functions. Should be deleted later
  const MatrixRectangular* getX0();
  const MatrixRectangular* getX0p();
  const MatrixRectangular* getY0();
  const MatrixRectangular* getY0p();
  const MatrixRectangular* getSigma0();
  const MatrixRectangular* getSigma0p();

  void resetLinkedToZ();
  void resetLinkedToLHS();
  void resetLinkedToRHS();
  void resetLinkedtoVar0();
  void resetLinkedToBayes();
  void resetLinkedToColCok();
  void resetLinkedToXvalid();

private:
  static bool _checkDimensionMatrix(const String& name,
                                               const AMatrix* mat,
                                               int* nrowsRef,
                                               int* ncolsRef);
  static bool _checkDimensionVector(const String& name,
                                    const VectorDouble* vec,
                                    int* sizeRef);

  static bool _isPresentMatrix(const String& name, const AMatrix* mat);
  static bool _isPresentVector(const String& name, const VectorDouble* vec);
  static bool _isPresentIVector(const String& name, const VectorInt* vec);

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
  int _needZ();
  int _needZp();
  int _needColCok();
  int _needXvalid();
  int _needDual();

  int _patchRHSForXvalidUnique();
  int _patchColCokVarianceZstar(MatrixSquareSymmetric* varZK);

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
  void _deleteZ();
  void _deleteZp();
  void _deleteColCok();
  void _deleteXvalid();
  void _deleteDual();

  static void _printMatrix(const String& name, const AMatrix* mat);
  static void _printVector(const String& name, const VectorDouble* vec);

  bool _validForDual() const;
  void _resetAll();

private:
  // Following information should not be removed in destructor
  const MatrixSquareSymmetric* _Sigma00;  // Variance at Target (Dim: _nrhs * _nrhs)
  const MatrixSquareSymmetric* _Sigma;    // Covariance Matrix (Dim: _neq * _neq)
  const MatrixRectangular* _Sigma0;       // Covariance at Target (Dim: _neq * _nrhs)
  const MatrixRectangular* _X;            // Drift at Data (Dim: _neq * _nbfl)
  const MatrixRectangular* _X0;           // Drift at Target (Dim: _nrhs * _nbfl)
  const MatrixSquareSymmetric* _PriorCov; // Bayesian Prior Covariance (Dim: _nbfl * _nbfl)
  const VectorDouble* _Z;                 // Data [flattened] (Dim: _neq)
  const VectorDouble* _PriorMean;         // Prior Bayesian Mean (Dim: _nbfl)
  const VectorDouble* _Means;             // Fixed drift coefficients
  const VectorDouble* _Zp;                // Vector of values for collocation
  const VectorInt* _rankColCok;           // Ranks of collocated variables
  const VectorInt* _rankXvalidEqs;        // Ranks of the cross-validated Equations
  const VectorInt* _rankXvalidVars;       // Ranks of the cross-validated Variables

  // Following elements can be retrieved by Interface functions  
  VectorDouble _Zstar;                  // Estimated values (Dim: _nrhs)
  VectorDouble _Beta;                   // Drift coefficients (Dim: _nbfl)
  MatrixRectangular* _LambdaSK;         // Weights for SK (Dim: _neq * _nrhs)
  MatrixRectangular* _LambdaUK;         // Weights for UK (Dim: _neq * _nrhs)
  MatrixRectangular* _MuUK;             // Lagrange multipliers (Dim: _nbfl * _nrhs)
  MatrixSquareSymmetric* _Stdv;         // Estimation stdv. (Dim: _nrhs * _nrhs)
  MatrixSquareSymmetric* _VarZSK;       // Estimator variance in SK (Dim: _nrhs * _nrhs)
  MatrixSquareSymmetric* _VarZUK;       // Estimator variance in UK (Dim: _nrhs * _nrhs)

  // Following elements are defined for internal storage
  MatrixRectangular* _XtInvSigma;       // X^t * Inv{Sigma} (Dim: _nbfl * _neq);
  MatrixRectangular* _Y0;               // X0 - LambdaSK * X^t (Dim: _nrhs * _nbfl)
  MatrixRectangular* _InvSigmaSigma0;   // Inv{Sigma} * Sigma0 (Dim: _neq * _nrhs)
  MatrixSquareSymmetric* _InvSigma;     // Inv{Sigma} (Dim: _neq * _neq)
  MatrixSquareSymmetric* _Sigmac;       // Inv{X^t * Inv{Sigma} * X} (Dim: _nbfl * _nbfl)
  MatrixSquareSymmetric* _InvPriorCov;  // Inv{PriorCov} (Dim: _nbfl * _nbfl)

  // Following elements are defined for internal storage (collocated case in UN)
  MatrixSquareSymmetric* _Sigma00pp; // ColCok Variance T-T (Dim: _ncck * _ncck)
  MatrixRectangular* _Sigma00p;      // ColCok Variance D-T (Dim: _ncck * _nrhs)
  MatrixRectangular* _Sigma0p;       // Collocated Covariance (Dim: _neq * _ncck)
  MatrixRectangular* _X0p;           // Collocated Drift (Dim: _ncck * _nbfl)
  MatrixRectangular* _Y0p;           // X0p - Sigma0p^t * Inv{Sigma} * X (Dim: _ncck *_nbfl)
  VectorDouble _Z0p;                 // Vector of (active) collocated values
  MatrixRectangular* _Lambda0;       // Collocated weights (Dim: _ncck * _nrhs)

  // Following elements are defined for Dual programming
  VectorDouble _bDual;               // Fake Covariance part in Dual (Dim: _neq)
  VectorDouble _cDual;               // Fake Drift part in Dual (Dim: _nbfl) 
   
  // Following elements are defined for internal storage (Cross-validation in UN)
  MatrixRectangular* _C_RHS;         // Fictitious Right-hand side (covariance part)
  MatrixRectangular* _X_RHS;         // Fictitious Right-hand side (drift part)
  
  // Additional parameters
  int _neq;
  int _nbfl;
  int _nrhs;
  int _ncck;
  int _nxvalid;
  bool _flagSK;
  bool _flagBayes;
  bool _flagDual;
};
