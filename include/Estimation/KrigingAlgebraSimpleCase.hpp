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

#include "Basic/VectorNumT.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "gstlearn_export.hpp"
#include "Basic/OptCustom.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/AMatrix.hpp"
#include <memory>

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
class GSTLEARN_EXPORT KrigingAlgebraSimpleCase
{
public:
  KrigingAlgebraSimpleCase(bool flagDual                        = false,
                           const VectorVectorInt* sampleRanks   = nullptr,
                           const VectorDouble* Z                = nullptr,
                           const VectorDouble& Means            = VectorDouble(),
                           int flagchol                         = false,
                           bool neighUnique                     = OptCustom::query("unique", 1));
  KrigingAlgebraSimpleCase(KrigingAlgebraSimpleCase& r);
  KrigingAlgebraSimpleCase& operator=(const KrigingAlgebraSimpleCase& r) = delete;
  virtual ~KrigingAlgebraSimpleCase();
  int prepare();
  void setDual(bool status);
  void setNeighUnique(bool nu = false) { _neighUnique = nu; }
  void resetNewData();
  int setData(const VectorDouble* Z          = nullptr,
              const VectorVectorInt* indices = nullptr,
              const VectorDouble& Means      = VectorDouble());
  int setLHS(const MatrixSymmetric* Sigma = nullptr,
             const MatrixDense* X         = nullptr);
  int setRHS(MatrixDense* Sigma0 = nullptr,
             MatrixDense* X0     = nullptr);
  int setVariance(const MatrixSymmetric* Sigma00 = nullptr);

  void printStatus() const;
  void dumpLHS(int nbypas = 5) const;
  void dumpRHS() const;
  void dumpWGT();
  void dumpAux();

  VectorDouble& getEstimation();
  VectorDouble getStdv();
  double getVarianceZstar(int i);
  VectorDouble getVarianceZstar();
  const MatrixSymmetric* getStdvMat();
  const MatrixSymmetric* getVarianceZstarMat();
  const MatrixDense* getLambda();
  const MatrixDense* getLambda0();
  const MatrixDense* getMu();
  double getLTerm();
  int updateRHS();

  bool isDual() const { return _flagDual; }
  VectorDouble* getZ() { return _Z.get(); }
  MatrixSquareSymmetric* getSigma() { return _Sigma.get(); }
  MatrixSquareSymmetric* getSigma00() { return _Sigma00.get(); }
  MatrixRectangular* getX() { return _X.get(); }
  MatrixDense* getSigma0() { return _Sigma0.get(); }
  void updateSampleRanks();
  MatrixDense* getX0() { return _X0.get(); }
  VectorVectorInt* getSampleRanks() { return _sampleRanks.get(); }
  VectorInt* getNbgh() { return _nbgh.get(); }
  void setMeans(const VectorDouble& means);

private:
  void _copyPtrForUniqueNeigh(KrigingAlgebraSimpleCase& r);
  void _copyContentForMovingNeigh(const KrigingAlgebraSimpleCase& r);
  void _copyOtherContent(const KrigingAlgebraSimpleCase& r);
  void _copyMatsAndVecs(const KrigingAlgebraSimpleCase& r);
  void _copyFlags(const KrigingAlgebraSimpleCase& r);
  void _copyModelQuantities(const KrigingAlgebraSimpleCase& r);

  static bool _checkDimensionMatrix(const String& name, const AMatrix* mat, int* nrowsRef, int* ncolsRef);
  static bool _checkDimensionVD(const String& name, const VectorDouble* vec, int* sizeRef);
  static bool _checkDimensionVI(const String& name, const VectorInt* vec, int* sizeRef);
  static bool _checkDimensionVVI(const String& name, const VectorVectorInt* vec, int* size1Ref, int* size2Ref);

  static bool _isPresentMatrix(const String& name, const AMatrix* mat);
  static bool _isPresentVector(const String& name, const VectorDouble* vec);
  static bool _isPresentIVector(const String& name, const VectorInt* vec);
  static bool _isPresentIIVector(const String& name, const VectorVectorInt* vec);

  void _resetLinkedToLHS();
  void _resetLinkedToRHS();

  bool _notFindZ();
  bool _notFindX();
  bool _notFindX0();
  bool _notFindSigma();
  bool _notFindSigma0();
  bool _notFindSigma00();

  bool _notFindSampleRanks();

  void _resetLinkedToZ();
  void _resetLinkedToX();
  void _resetLinkedToX0();
  void _resetLinkedToSigma();
  void _resetLinkedToSigma0();
  void _resetLinkedToSigma00();

  int _needInvSigma();
  int _needXtInvSigma();
  int _needXtInvSigmaZ();
  int _needSigmac();
  int _needBeta();
  int _needDual();

  int _needLambdaSK();
  int _needLambdaUK();
  int _needMuUK();

  int _needZstar();
  int _needStdv();
  int _needVarZSK();
  int _needVarZUK();

  int _computeZstarWithDual();
  int _computeZstarSK();

  void _deleteBeta();
  void _deleteInvSigma();
  void _deleteLambdaSK();
  void _deleteLambdaUK();
  void _deleteMuUK();
  void _deleteSigmac();
  void _deleteXtInvSigma();
  void _deleteXtInvSigmaZ();
  void _deleteStdv();
  void _deleteVarZSK();
  void _deleteVarZUK();
  void _deleteDual();

  static void _printMatrix(const String& name, const AMatrix* mat);
  static void _printVector(const String& name, const VectorDouble* vec);

  bool _forbiddenWhenDual() const;
  void _resetAll();

private:
  // Quantities to be defined by the user
  std::shared_ptr<VectorDouble> _Z;              // Data [flattened] (Dim: _neq)
  std::shared_ptr<VectorVectorInt> _sampleRanks; // Vector of Vector of sampl indices per variable
  std::shared_ptr<VectorInt> _nbgh;
  std::shared_ptr<MatrixDense> _X;           // Drift at Data (Dim: _neq * _nbfl)
  std::shared_ptr<MatrixSymmetric> _Sigma;   // Covariance Matrix (Dim: _neq * _neq)
  std::shared_ptr<MatrixSymmetric> _Sigma00; // Variance at Target (Dim: _nrhs * _nrhs)
  std::shared_ptr<MatrixDense> _Sigma0;      // Covariance at Target (Dim: _neq * _nrhs)
  std::shared_ptr<MatrixDense> _X0;          // Drift at Target (Dim: _nrhs * _nbfl)

  VectorDouble _Means; // Fixed drift coefficients

  std::shared_ptr<MatrixSymmetric> _InvSigma; // Inv{Sigma} (Dim: _neq * _neq)
  std::shared_ptr<CholeskyDense> _cholSigma;
  std::shared_ptr<MatrixDense> _XtInvSigma;    // X^t * Inv{Sigma} (Dim: _nbfl * _neq);
  std::shared_ptr<MatrixDense> _invSigmaX;     // Inv{Sigma} X (Dim: _neq * _nbfl);
  std::shared_ptr<VectorDouble> _XtInvSigmaZ;        // X^t * Inv{Sigma} Z (Dim: _nbfl * _nvar);
  std::shared_ptr<MatrixSymmetric> _invSigmac; // Inv{X^t * Inv{Sigma} * X} (Dim: _nbfl * _nbfl)
  std::shared_ptr<VectorDouble> _Beta;               // Drift coefficients (Dim: _nbfl)
  std::shared_ptr<MatrixDense> _LambdaSK;      // Weights for SK (Dim: _neq * _nrhs)
                                                     // Following elements are defined for Dual programming
  std::shared_ptr<VectorDouble> _bDual;              // Fake Covariance part in Dual (Dim: _neq)
  std::shared_ptr<VectorDouble> _invSigmaXBeta;

  VectorDouble _Zstar; // Estimated values (Dim: _nrhs)

  MatrixDense _LambdaSKtX;
  MatrixDense _LambdaUK;   // Weights for UK (Dim: _neq * _nrhs)
  MatrixDense _MuUK;       // Lagrange multipliers (Dim: _nbfl * _nrhs)
  MatrixSymmetric _Stdv;   // Estimation stdv. (Dim: _nrhs * _nrhs)
  MatrixSymmetric _VarZSK; // Estimator variance in SK (Dim: _nrhs * _nrhs)
  MatrixSymmetric _VarZUK; // Estimator variance in UK (Dim: _nrhs * _nrhs)
  MatrixSymmetric _Sigmac;
  // Following elements are defined for internal storage

  MatrixDense _Y0; // X0 - LambdaSK * X^t (Dim: _nrhs * _nbfl)

  MatrixDense _LambdaUKtSigma0;
  MatrixDense _MuUKtX0t;
  MatrixDense _invSigmaXMuUK;
  VectorDouble _X0Beta;

  // Additional parameters
  int _nvar;
  int _neq;
  int _nbfl;
  int _nrhs;

  //  Flags
  bool _flagSK;
  bool _flagDual;
  bool _neighUnique;
  bool _flagCholesky;
  bool _dualHasChanged;
  bool _invSigmaHasChanged;
  bool _XtInvSigmaHasChanged;
};
