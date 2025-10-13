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

#include "Enum/EKrigOpt.hpp"
#include "Estimation/KrigOpt.hpp"
#include "Estimation/KrigingAlgebra.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixSquare.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Model/ModelGeneric.hpp"
#include "Neigh/ANeigh.hpp"
#include "Space/SpacePoint.hpp"
#include "Space/SpaceRN.hpp"

namespace gstlrn
{

class Db;
class DbGrid;
class Model;
class ModelGeneric;
class ANeigh;
class CovCalcMode;
class ECalcMember;
class NeighImage;
class AAnam;
class ACov;
class KrigingAlgebra;
class KrigOpt;

class GSTLEARN_EXPORT KrigingSystem
{
public:
  KrigingSystem(Db* dbin,
                Db* dbout,
                const ModelGeneric* model,
                ANeigh* neigh,
                const KrigOpt& krigopt = KrigOpt());
  KrigingSystem(const KrigingSystem& m)            = delete;
  KrigingSystem& operator=(const KrigingSystem& m) = delete;
  virtual ~KrigingSystem();

  Id resetData();
  Id setKrigOpt(const KrigOpt& krigopt);
  Id setKrigOptCalcul(const EKrigOpt& calcul,
                      const VectorInt& ndiscs = VectorInt(),
                      bool flag_per_cell      = false);
  Id setKrigOptXValid(bool flag_xvalid,
                      bool flag_kfold,
                      bool optionXValidEstim = false,
                      bool optionXValidStdev = false,
                      bool optionXValidVarZ  = false);
  Id setKrigOptBayes(bool flag_bayes,
                     const VectorDouble& prior_mean,
                     const MatrixSymmetric& prior_cov);
  Id setKrigOptDataWeights(Id iptrWeights, bool flagSet = true);
  Id setKrigOptFlagSimu(bool flagSimu, Id nbsimu = 0, Id rankPGS = -1);
  Id setKrigOptFlagGlobal(bool flag_global);
  Id setKrigOptFlagLTerm(bool flag_lterm);
  Id setKrigOptAnamophosis(AAnam* anam);
  Id setKrigOptFactorKriging(bool flag_factor_kriging);

  // The subsequent methods do not require isReady() validation
  Id updKrigOptEstim(Id iptrEst, Id iptrStd, Id iptrVarZ, bool forceNoDual = false);
  Id updKrigOptIclass(Id index_class, Id nclasses);
  Id updKrigOptNeighOnly(Id iptrNeigh);
  bool isReady();
  Id estimate(Id iech_out);
  void conclusion();

  // Methods used to return algebraic internal information
  Id getNDim() const { return (_model != nullptr) ? static_cast<Id>(_model->getNDim()) : 0; }
  Id getNVar() const { return (_model != nullptr) ? _model->getNVar() : 0; }
  Id getNech() const { return static_cast<Id>(_nbgh.size()); }
  Id getCovSize() const { return (!_Sigma.empty()) ? _Sigma.getNRows() : 0; }
  Id getDriftSize() const { return (!_X.empty()) ? _X.getNCols() : 0; }
  Id getNrhs() const { return (!_Sigma0.empty()) ? _Sigma0.getNCols() : 0; }

  VectorInt getSampleNbgh() const { return _nbgh; }
  VectorVectorDouble getSampleCoordinates() const;
  VectorDouble getSampleData() const { return _Z; };
  MatrixSymmetric getLHS() const { return _Sigma; }
  MatrixDense getLHSF() const { return _Sigma0; }
  MatrixDense getRHS() const { return _Sigma0; }
  MatrixDense getRHSF() const { return _X0; }
  MatrixSquare getVariance() const { return _Sigma00; }
  MatrixDense getWeights() const;
  MatrixDense getMu() const;
  double getLTerm() const { return _algebra.getLTerm(); }

private:
  Id _getNVar() const;
  Id _getNVarCL() const;
  Id _getNbfl() const;
  Id _getNeq() const;
  Id _getNFeq() const;

  void _resetMemoryGeneral();
  bool _isAuthorized() const;

  void _rhsDump();
  void _wgtDump();
  void _estimateCalcul(Id status);
  void _simulateCalcul(Id status);
  void _neighCalcul(Id status, const VectorDouble& tab);
  void _estimateVarZ(Id status);
  void _estimateStdv(Id status);
  void _estimateEstim(Id status);
  void _dumpKrigingResults(Id status);
  void _dumpSimulationResults(Id status);
  bool _isCorrect();
  bool _preparNoStat();

  Id _bayesPreCalculations();
  void _bayesPreSimulate();
  void _transformGaussianToRaw();

  void _setInternalShortCutVariablesGeneral();
  void _setInternalShortCutVariablesModel();
  Id _setInternalShortCutVariablesNeigh();

  Model* _castInOldModel();
  VectorInt _xvalidUniqueIndices() const;
  Id _updateForColCokMoving();

  // Deprecated function
  double _continuousMultiplier(Id rank1, Id rank2, double eps = EPSILON4);

private:
  Db* _dbin;
  Db* _dbout;
  ModelGeneric* _model;
  ANeigh* _neigh;
  const AAnam* _anam;
  bool _isReady;

  // Pointers used when plugging KrigingAlgebra (not to be deleted)
  // Note that 'algebra' is mutable not to destroy constness when calling getLambda.
  mutable KrigingAlgebra _algebra;
  mutable KrigOpt _krigopt;
  VectorVectorInt _sampleRanks; // Vector of vector of sample indices
  MatrixSymmetric _Sigma00;     // Covariance part for variance
  MatrixSymmetric _Sigma;       // Covariance part for LHS
  MatrixDense _X;               // Drift part for LHS
  MatrixDense _Sigma0;          // Covariance part for RHS
  MatrixDense _X0;              // Drift par for RHS
  VectorDouble _Z;              // Vector of Data
  VectorDouble _means;          // Means of the variables (used to center variables)
  VectorDouble _meansTarget;    // Means for target (possible using matLC)

  /// UID for storage
  Id _iptrEst;
  Id _iptrStd;
  Id _iptrVarZ;
  bool _flagEst;
  bool _flagStd;
  bool _flagVarZ;
  bool _flagDataChanged;

  /// Option for Weights at Data locations
  Id _iptrWeights;
  bool _flagWeights;
  bool _flagSet;

  /// Option for Simulation
  bool _flagSimu;
  Id _nbsimu;
  Id _rankPGS;

  /// Option for Cross_validation
  bool _xvalidEstim;
  bool _xvalidStdev;
  bool _xvalidVarZ;

  /// Option for Colocation
  VectorDouble _valuesColcok;

  /// Option for Bayesian
  bool _flagBayes;
  VectorDouble _priorMean;
  MatrixSymmetric _priorCov;
  VectorDouble _postMean;
  MatrixSymmetric _postCov;
  MatrixDense _postSimu;
  MatrixSymmetric _varCorrec;

  /// Option for (Disjunctive) Kriging of Factor
  bool _flagFactorKriging;
  Id _nclasses;
  Id _factorClass;

  /// Option for asking for Z * A-1 * Z
  bool _flagLTerm;

  /// Option for Gaussian Kriging
  bool _flagAnam;

  /// Option for Neighboring test
  bool _flagNeighOnly;
  Id _iptrNeigh;

  /// Local variables
  Id _iechOut;
  Id _ndim;
  Id _nvar;
  Id _nvarCL;
  Id _nech;
  Id _nfeq;
  Id _neq;

  /// Working arrays
  mutable VectorInt _nbgh;
  mutable VectorInt _dbinUidToBeDeleted;
  mutable VectorInt _dboutUidToBeDeleted;

  /// Some Space Point allocated once for all
  mutable ASpaceSharedPtr _space;
  mutable SpacePoint _p0;
  mutable SpacePoint _p1;
  mutable SpacePoint _p2;
  mutable SpacePoint _p0_memo;

  /// Some local flags defined in order to speed up the process
  mutable bool _flagVerr;
  mutable bool _flagNoStat;

  /// Some local variables for XValid
  mutable VectorInt _rankXvalidEqs;
  mutable VectorInt _rankXvalidVars;
};
} // namespace gstlrn