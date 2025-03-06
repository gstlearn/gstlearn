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

#include "Estimation/KrigingAlgebra.hpp"
#include "Estimation/KrigOpt.hpp"
#include "Model/ModelGeneric.hpp"
#include "Space/SpaceRN.hpp"
#include "Space/SpacePoint.hpp"
#include "Neigh/ANeigh.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Enum/EKrigOpt.hpp"

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
                ANeigh* neigh);
  KrigingSystem(const KrigingSystem &m) = delete;
  KrigingSystem& operator=(const KrigingSystem &m) = delete;
  virtual ~KrigingSystem();

  int resetData();
  int setKrigOptCalcul(const EKrigOpt& calcul,
                       const VectorInt& ndiscs = VectorInt(),
                       bool flag_per_cell      = false);
  int setKrigOptXValid(bool flag_xvalid,
                       bool flag_kfold,
                       bool optionXValidEstim = false,
                       bool optionXValidStdev = false,
                       bool optionXValidVarZ  = false);
  int setKrigOptColCok(const VectorInt& rank_colcok);
  int setKrigOptBayes(bool flag_bayes,
                      const VectorDouble& prior_mean,
                      const MatrixSquareSymmetric& prior_cov);
  int setKrigOptDataWeights(int iptrWeights, bool flagSet = true);
  int setKrigOptMatLC(const MatrixRectangular* matLC);
  int setKrigOptFlagSimu(bool flagSimu, int nbsimu = 0, int rankPGS = -1);
  int setKrigOptDGM(bool flag_dgm, double eps = EPSILON6);
  int setKrigOptFlagGlobal(bool flag_global);
  int setKrigOptFlagLTerm(bool flag_lterm);
  int setKrigOptAnamophosis(AAnam* anam);
  int setKrigOptFactorKriging(bool flag_factor_kriging);

  // The subsequent methods do not require isReady() validation
  int  updKrigOptEstim(int iptrEst, int iptrStd, int iptrVarZ, bool forceNoDual = false);
  int  updKrigOptIclass(int index_class, int nclasses);
  int  updKrigOptNeighOnly(int iptrNeigh);
  bool isReady();
  int  estimate(int iech_out);
  void conclusion();

  // Methods used to return algebraic internal information
  int getNDim() const { return (_model != nullptr) ? _model->getNDim() : 0; }
  int getNVar() const { return (_model != nullptr) ? _model->getNVar() : 0; }
  int getNech() const { return (int)_nbgh.size(); }
  int getCovSize() const { return (!_Sigma.empty()) ? _Sigma.getNRows() : 0; }
  int getDriftSize() const { return (!_X.empty()) ? _X.getNCols() : 0; }
  int getNrhs() const { return (!_Sigma0.empty()) ? _Sigma0.getNCols() : 0; }

  VectorInt             getSampleNbgh() const { return _nbgh; }
  VectorVectorDouble    getSampleCoordinates() const;
  VectorDouble          getSampleData() const { return _Z; };
  MatrixSquareSymmetric getLHS() const { return _Sigma; }
  MatrixRectangular     getLHSF() const { return _Sigma0; }
  MatrixRectangular     getRHS() const { return _Sigma0; }
  MatrixRectangular     getRHSF() const { return _X0; }
  MatrixSquareGeneral   getVariance() const { return _Sigma00; }
  MatrixRectangular     getWeights() const;
  MatrixRectangular     getMu() const;
  double                getLTerm() const { return _algebra.getLTerm(); }

private:
  int    _getNVar() const;
  int    _getNVarCL() const;
  int    _getNbfl() const;
  int    _getNeq() const;
  int    _getNFeq() const;

  void _resetMemoryGeneral();
  bool _isAuthorized() const;

  void _dumpOptions() const;
  void _rhsDump();
  void _wgtDump();
  void _estimateCalcul(int status);
  void _simulateCalcul(int status);
  void _neighCalcul(int status, const VectorDouble& tab);
  void _estimateVarZ(int status);
  void _estimateStdv(int status);
  void _estimateEstim(int status);
  void _dumpKrigingResults(int status);
  void _dumpSimulationResults(int status);
  bool _isCorrect();
  bool _preparNoStat();

  int    _bayesPreCalculations();
  void   _bayesPreSimulate();
  void   _transformGaussianToRaw();

  void   _setInternalShortCutVariablesGeneral();
  void   _setInternalShortCutVariablesModel();
  int    _setInternalShortCutVariablesNeigh();

  Model* _castInOldModel();
  VectorInt _xvalidUniqueIndices() const;
  int  _updateForColCokMoving();

  // Deprecated function
  double _continuousMultiplier(int rank1, int rank2, double eps = EPSILON4);

private:
  Db* _dbin;
  Db* _dbout;
  ModelGeneric* _model;
  ANeigh*       _neigh;
  const AAnam*  _anam;
  bool          _isReady;

  // Pointers used when plugging KrigingAlgebra (not to be deleted)
  // Note that 'algebra' is mutable not to destroy constness when calling getLambda.
  mutable KrigingAlgebra _algebra;
  mutable KrigOpt        _krigopt;
  VectorVectorInt        _sampleRanks; // Vector of vector of sample indices
  MatrixSquareSymmetric  _Sigma00; // Covariance part for variance
  MatrixSquareSymmetric  _Sigma;   // Covariance part for LHS
  MatrixRectangular      _X;       // Drift part for LHS
  MatrixRectangular      _Sigma0;  // Covariance part for RHS
  MatrixRectangular      _X0;      // Drift par for RHS
  VectorDouble           _Z;       // Vector of Data
  VectorDouble _means;            // Means of the variables (used to center variables)
  VectorDouble _meansTarget;      // Means for target (possible using matLC)

  /// UID for storage
  int  _iptrEst;
  int  _iptrStd;
  int  _iptrVarZ;
  bool _flagEst;
  bool _flagStd;
  bool _flagVarZ;
  bool _flagDataChanged;

  /// Option for Calculation
  EKrigOpt _calcul;

  /// Option for Weights at Data locations
  int  _iptrWeights;
  bool _flagWeights;
  bool _flagSet;

  /// Option for Simulation
  bool _flagSimu;
  int  _nbsimu;
  int  _rankPGS;

  /// Option for Block estimation
  VectorInt _ndiscs;

  /// Option for Cross_validation
  bool _xvalidEstim;
  bool _xvalidStdev;
  bool _xvalidVarZ;

  /// Option for Colocation
  VectorInt _rankColCok;
  VectorDouble _valuesColCok;

  /// Option for Bayesian
  bool _flagBayes;
  VectorDouble          _priorMean; 
  MatrixSquareSymmetric _priorCov;  
  VectorDouble          _postMean;
  MatrixSquareSymmetric _postCov;
  MatrixRectangular     _postSimu; 
  MatrixSquareSymmetric _varCorrec;

  /// Option for Discrete Gaussian case
  bool _flagDGM;

  /// Option for (Disjunctive) Kriging of Factor
  bool _flagFactorKriging;
  int  _nclasses;
  int  _factorClass;

  /// Option for Estimating the Linear Combination of Variables
  const MatrixRectangular* _matLC;

  /// Option for asking for Z * A-1 * Z
  bool _flagLTerm;

  /// Option for Gaussian Kriging
  bool _flagAnam;

  /// Option for Neighboring test
  bool _flagNeighOnly;
  int _iptrNeigh;

  /// Local variables
  int _iechOut;
  int _ndim;
  int _nvar;
  int _nvarCL;
  int _nech;
  int _nfeq;
  int _neq;
  int _nred;

  /// Working arrays
  mutable VectorInt    _nbgh;
  mutable VectorInt    _dbinUidToBeDeleted;
  mutable VectorInt    _dboutUidToBeDeleted;

  /// Some Space Point allocated once for all
  mutable ASpaceSharedPtr    _space;
  mutable SpacePoint _p0;
  mutable SpacePoint _p1;
  mutable SpacePoint _p2;
  mutable SpacePoint _p0_memo;

  /// Some local flags defined in order to speed up the process
  mutable bool _flagNoMatLC;
  mutable bool _flagVerr;
  mutable bool _flagNoStat;
};
