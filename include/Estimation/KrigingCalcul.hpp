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

class GSTLEARN_EXPORT KrigingCalcul
{
public:
  KrigingCalcul(const MatrixSquareSymmetric* C00 = nullptr,
                const MatrixSquareSymmetric* C = nullptr,
                const MatrixRectangular* X = nullptr);
  KrigingCalcul(const KrigingCalcul &r) = delete;
  KrigingCalcul& operator=(const KrigingCalcul &r) = delete;
  virtual ~KrigingCalcul();

  int setC(const MatrixSquareSymmetric* C);
  int setX(const MatrixRectangular* X);
  int setC0(const MatrixRectangular* C0);
  int setX0(const MatrixRectangular* X0);
  int setZ(const VectorDouble& Z);
  int setBeta(const VectorDouble& beta);
  int setBayes(const VectorDouble& mBayes, const MatrixSquareSymmetric* SBayes);

private:
  static bool _matchDimensions(const AMatrix* mat, int nrowsRef, int ncolsRef);
  int _computeXi();
  int _computeLambdaSK();
  int _computeLambdaUK();
  int _computeZstar(bool flagSK);
  int _computeVarianceZ(bool flagSK);
  int _computeVariance(bool flagSK);
  int _computeXtCm1();
  int _computeLambdaSKtX();
  int _computeX0mLambdaSKtX();
  int _computeVarianceBeta();
  int _computeBeta();
  int _computeCm1();

  int _needCm1();
  int _needLambdaSK();
  int _needLambdaUK();
  int _needLambdaSKtX();
  int _needVarBeta();
  int _needX0mLambdaSKtX();
  int _needXtCm1();
  int _needVarZ(bool flagSK);

private:
  const MatrixSquareSymmetric* _C00;    // Not to be deleted (Dim: _nrhs * _nrhs)
  const MatrixSquareSymmetric* _C;      // Not to be deleted (Dim: _neq * _neq)
  const MatrixRectangular* _X;          // Not to be deleted (Dim: _neq * _nbfl)
  const MatrixRectangular* _C0;         // Not to be deleted (Dim: _neq * _nrhs)
  const MatrixRectangular* _X0;         // Not to be deleted (Dim: _nbfl * _nrhs)
  const MatrixSquareSymmetric* _SBayes; // Not to be deleted (Dim: _nbfl * _nbfl)
  VectorDouble _Z;                      // Vector of data (Dim: _neq)
  VectorDouble _mBayes;                 // Vector of prior mean (Dim: _nbfl)
  VectorDouble _Xi;                     // Vector of drift at data (Dim: _neq)
  VectorDouble _Beta;                   // Vector of drift coefficients (Dim: _nbfl)
  VectorDouble _Zstar;                  // Vector of estimated values (Dim: _nrhs)
  MatrixSquareSymmetric* _Cm1;          // Inverse of Cov matrix (Dim: _neq * _neq)
  MatrixRectangular* _lambdaSK;         // Vector of Weights for SK (Dim: _neq * _nrhs)
  MatrixRectangular* _lambdaUK;         // Vector of Weights for UK (Dim: _neq * _nrhs)
  MatrixRectangular* _XtCm1;            // Product: Xt * Cm1 (Dim: _nbfl * _neq);
  MatrixRectangular* _LambdaSKtX;       // Product: Xt * LambdaSK (Dim: _nbfl * _nrhs)
  MatrixRectangular* _X0mLambdaSKtX;    // X0 - LambdaSK * Xt (Dim: _nbfl * _nrhs)
  MatrixSquareSymmetric* _var;          // Matrix of estimation var. (Dim: _nrhs * _nrhs)
  MatrixSquareSymmetric* _varZ;         // Matrix of estimator var. (Dim: _nrhs * _nrhs)
  MatrixSquareSymmetric* _varBeta;      // Variance of Beta (Dim: _nbfl * _nbfl)
  int _neq;
  int _nbfl;
  int _nrhs;
};
