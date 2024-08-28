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
  KrigingCalcul(const VectorDouble& Z                = VectorDouble(),
                const MatrixSquareSymmetric* Sigma   = nullptr,
                const MatrixRectangular* X           = nullptr,
                const MatrixSquareSymmetric* Sigma00 = nullptr,
                const VectorDouble& Beta             = VectorDouble());
  KrigingCalcul(const KrigingCalcul &r) = delete;
  KrigingCalcul& operator=(const KrigingCalcul &r) = delete;
  virtual ~KrigingCalcul();

  int setSigma00(const MatrixSquareSymmetric* Sigma00);
  int setSigma(const MatrixSquareSymmetric* Sigma);
  int setSigma0(const MatrixRectangular* Sigma0);
  int setX(const MatrixRectangular* X);
  int setX0(const MatrixRectangular* X0);
  int setZ(const VectorDouble& Z);
  int setBeta(const VectorDouble& beta);
  int setColCok(const VectorInt& rankColCok);
  int setBayes(const VectorDouble& PriorMean, const MatrixSquareSymmetric* PriorCov);

  void printStatus() const;

  VectorDouble getEstimation();
  VectorDouble getStdv();
  VectorDouble getVarianceZstar();
  VectorDouble getPostMean();
  const MatrixSquareSymmetric* getStdvMat();
  const MatrixSquareSymmetric* getVarianceZstarMat();
  const MatrixSquareSymmetric* getPostCov();
  const MatrixRectangular*     getLambdaSK();
  const MatrixRectangular*     getLambdaUK();
  const MatrixRectangular*     getMu();

private:
  static bool _checkDimensionMatrix(const String& name,
                                    const AMatrix* mat,
                                    int* nrowsRef,
                                    int* ncolsRef);
  static bool _checkDimensionVector(const String& name,
                                    const VectorDouble& vec,
                                    int* sizeRef);

  static bool _isPresentMatrix(const String& name, const AMatrix* mat);
  static bool _isPresentVector(const String& name, const VectorDouble& vec);
  static bool _isPresentIVector(const String& name, const VectorInt& vec);

  int _needX();
  int _needX0();
  int _needSigma();
  int _needSigma0();
  int _needBeta();
  int _needInvSigma();
  int _needLambdaSK();
  int _needLambdaUK();
  int _needMu();
  int _needSigmac();
  int _needZstar();
  int _needVarSK();
  int _needY0();
  int _needXtInvSigma();
  int _needStdv();
  int _needVarZSK();
  int _needVarZUK();
  int _needInvPriorCov();
  int _needInvCCK();
  int _needSigma0p();
  int _needSigma0pInvCCK();
  int _needX0p();
  int _needY0p();
  int _needX0ptInvCCK();
  int _needLambda0();

  static void _printMatrix(const String& name, const AMatrix* mat);
  static void _printVector(const String& name, const VectorDouble& vec);

private:
  // Following pointers should not be removed in destructor
  const MatrixSquareSymmetric* _Sigma00;  // Variance at Target (Dim: _nrhs * _nrhs)
  const MatrixSquareSymmetric* _Sigma;    // Covariance Matrix (Dim: _neq * _neq)
  const MatrixRectangular* _Sigma0;       // Covariance at Target (Dim: _neq * _nrhs)
  const MatrixRectangular* _X;            // Drift at Data (Dim: _neq * _nbfl)
  const MatrixRectangular* _X0;           // Drift at Target (Dim: _nrhs * _nbfl)
  const MatrixSquareSymmetric* _PriorCov; // Bayesian Prior Covariance (Dim: _nbfl * _nbfl)

  // Following elements can be retrieved by Interface functions  
  VectorDouble _Z;                      // Data [flattened] (Dim: _neq)
  VectorDouble _PriorMean;              // Prior Bayesian Mean (Dim: _nbfl)
  VectorDouble _Zstar;                  // Estimated values (Dim: _nrhs)
  VectorDouble _Beta;                   // Drift coefficients (Dim: _nbfl)
  MatrixRectangular* _LambdaSK;         // Weights for SK (Dim: _neq * _nrhs)
  MatrixRectangular* _LambdaUK;         // Weights for UK (Dim: _neq * _nrhs)
  MatrixRectangular* _Mu;               // Lagrange multipliers (Dim: _nbfl * _nrhs)
  MatrixSquareSymmetric* _Stdv;         // Estimation stdv. (Dim: _nrhs * _nrhs)
  MatrixSquareSymmetric* _VarSK;        // Estimation variance in SK (Dim: _nrhs * _nrhs)
  MatrixSquareSymmetric* _VarZSK;       // Estimator variance in SK (Dim: _nrhs * _nrhs)
  MatrixSquareSymmetric* _VarZUK;       // Estimator variance in UK (Dim: _nrhs * _nrhs)

  // Following elements are defined for internal storage
  MatrixRectangular* _XtInvSigma;       // X^t * InvSigma (Dim: _nbfl * _neq);
  MatrixRectangular* _Y0;               // X0 - LambdaSK * X^t (Dim: _nrhs * _nbfl)
  MatrixSquareSymmetric* _InvSigma;     // (Sigma)^{-1} (Dim: _neq * _neq)
  MatrixSquareSymmetric* _Sigmac;       // (X^t * Sigma^{-1} * X)^{-1} (Dim: _nbfl * _nbfl)
  MatrixSquareSymmetric* _InvPriorCov;  // (PriorCov)^{-1} (Dim: _nbfl * _nbfl)

  // Following elements are defined for internal storage (collocated case)
  MatrixRectangular* _Sigma0p;       // Collocated Covariance (Dim: _neq * _ncck)
  MatrixRectangular* _X0p;           // Collocated Drift (Dim: _ncck * _nbfl)
  MatrixRectangular* _Y0p;           // X0p - Sigma0p^t * InvSigma * X (Dim: _ncck * _nbfl)
  MatrixSquareSymmetric* _InvCCK;    // (C^00_kk)^{-1} (Dim: _ncck * _ncck)
  MatrixRectangular* _X0ptInvCCK;    // X0p^t * InvCCK (Dim: _nbfl * _ncck)
  MatrixRectangular* _Sigma0pInvCCK; // Sigmap0 * InvCCK (Dim: _neq * _ncck)

  // Additional parameters
  int _neq;
  int _nbfl;
  int _nrhs;
  int _ncck;
  bool _flagSK;
  bool _flagBayes;
  VectorInt _rankColCok;
};
