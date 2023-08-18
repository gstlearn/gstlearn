/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Space/SpaceRN.hpp"
#include "Space/SpacePoint.hpp"
#include "Neigh/ANeigh.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Enum/EKrigOpt.hpp"

class Db;
class DbGrid;
class Model;
class ANeigh;
class CovCalcMode;
class ECalcMember;
class NeighImage;
class AAnam;

class GSTLEARN_EXPORT KrigingSystem
{
public:
  KrigingSystem(Db* dbin,
                Db* dbout,
                const Model* model,
                ANeigh* neigh);
  KrigingSystem(const KrigingSystem &m) = delete;
  KrigingSystem& operator=(const KrigingSystem &m) = delete;
  virtual ~KrigingSystem();

  int  setKrigOptCalcul(const EKrigOpt& calcul,
                        const VectorInt& ndiscs = VectorInt(),
                        bool flag_per_cell = false);
  int  setKrigOptXValid(bool flag_xvalid,
                        bool flag_kfold,
                        bool optionXValidEstim = false,
                        bool optionXValidStdev = false,
                        bool optionXValidVarZ = false);
  int  setKrigOptColCok(const VectorInt& rank_colcok);
  int  setKrigOptBayes(bool flag_bayes,
                       const VectorDouble& prior_mean,
                       const VectorDouble& prior_cov,
                       int seed = 414371);
  int  setKrigOptImage(int seed = 133271);
  int  setKrigOptDataWeights(int iptrWeights, bool flagSet = true);
  int  setKrigOptMatCL(const VectorVectorDouble& matCL);
  int  setKrigoptCode(bool flag_code);
  int  setKrigOptFlagSimu(bool flagSimu, int nbsimu = 0, int rankPGS = -1);
  int  setKrigOptSaveWeights(bool flag_save);
  int  setKrigOptDGM(bool flag_dgm, double eps = EPSILON6);
  int  setKrigOptFlagGlobal(bool flag_global);
  int  setKrigOptFlagLTerm(bool flag_lterm);
  int  setKrigOptAnamophosis(AAnam* anam);
  int  setKrigOptFactorKriging(bool flag_factor_kriging);

  // The subsequent methods do not require isReady() validation
  int  updKrigOptEstim(int iptrEst, int iptrStd, int iptrVarZ);
  int  updKrigOptIclass(int index_class, int nclasses);
  int  updKrigOptNeighOnly(int iptrNeigh);

  void setFlagCheckAddress(bool flagCheckAddress) { _flagCheckAddress = flagCheckAddress; }

  bool isReady();
  int  estimate(int iech_out);
  void conclusion();

  int  getNDim() const;
  int  getNech() const;
  int  getNeq()  const;
  int  getNRed() const { return _nred; }
  VectorInt    getSampleIndices() const { return _nbgh; }
  VectorVectorDouble getSampleCoordinates() const;
  VectorDouble getSampleData() const;
  VectorDouble getZam() const { return _zam; }
  VectorDouble getLHS() const { return _lhs; }
  VectorDouble getLHSInv() const { return _lhsinv; }
  VectorDouble getRHSC() const { return _rhs; }
  VectorDouble getRHSC(int ivar) const;
  VectorDouble getWeights() const { return _wgt; }
  VectorDouble getVariance() const { return _var0.getValues(); }
  double getLTerm() const { return _lterm; }

private:
  int    _getNVar() const;
  int    _getNVarCL() const;
  int    _getNbfl() const;
  int    _getNFeq() const;
  int    _getNFex() const;
  int    _getNDisc() const;
  void   _setFlag(int iech, int ivar, int value);
  int    _getFlag(int iech, int ivar);
  double _getIdim(int loc_rank, int idim) const;
  double _getFext(int rank, int ibfl) const;
  double _getIvar(int rank, int ivar) const;
  double _getVerr(int rank, int ivar) const;
  double _getMean(int ivarCL) const;
  double _getCoefDrift(int ivar, int il, int ib) const;
  int    _getFLAG(int iech,int ivar) const;
  double _getCOVTAB(int ivar,int jvar) const;
  void   _addCOVTAB(int ivar,int jvar,double value);
  void   _prodCOVTAB(double value);
  double _getRHS(int iech, int ivar, int jvCL) const;
  void   _setRHS(int iech, int ivar, int jvCL, double value, bool isForDrift = false);
  double _getRHSC(int i, int jvCL) const;
  double _getWGTC(int i,int jvCL) const;
  double _getLHS(int iech, int ivar, int jech, int jvar) const;
  double _getLHSINV(int iech, int ivar, int jech, int jvar) const;
  void   _setLHS(int iech, int ivar, int jech, int jvar, double value, bool isForDrift = false);
  void   _addLHS(int iech, int ivar, int jech, int jvar, double value);
  double _getLHSC(int i, int j) const;
  double _getDISC1(int idisc, int idim) const;
  double _getZAM(int i) const;
  double _getZEXT(int i) const;
  void   _setZEXT(int i, double value) const;
  VectorDouble _getDISC1Vec(int idisc) const;
  VectorVectorDouble _getDISC1s() const;
  double _getDISC2(int idisc,int idim) const;
  VectorDouble _getDISC2Vec(int idisc) const;
  VectorVectorDouble _getDISC2s() const;
  double _getVAR0(int ivCL, int jvCL) const;
  void   _setVAR0(int ivCL, int jvCL, double value);

  const double* _getRHSCAdd(int i = 0, int jvCL = 0) const;
  const double* _getWGTCAdd(int i = 0, int jvCL = 0) const;
  const double* _getZamAdd(int i = 0) const;
  const double* _getZextAdd(int i = 0) const;

  void _resetMemoryGeneral();
  void _resetMemoryPerNeigh();
  void _flagDefine();
  void _zextInit();
  void _lhsInit();
  void _covUpdate(const ECalcMember &member, int iech1, int iech2);
  void _covtabInit();
  void _covtabCalcul(int iech1,
                     int iech2,
                     const CovCalcMode* mode,
                     bool flagSameData = false);
  void _covCvvCalcul(const CovCalcMode* mode);
  int  _drftabCalcul(const ECalcMember &member, int iech);
  bool _isAuthorized();
  double _continuousMultiplier(int rank1,int rank2, double eps = EPSILON4);
  void _lhsCalcul();
  void _lhsIsoToHetero();
  void _lhsDump(int nbypas = 5);
  int  _rhsCalcul();
  void _rhsCalculPoint();
  void _rhsCalculBlock();
  void _rhsCalculDrift();
  void _rhsCalculDGM();
  void _rhsStore(int iech);
  void _rhsIsoToHetero();
  void _rhsDump();
  void _wgtCalcul();
  void _wgtDump(int status);
  VectorInt _getRelativePosition();
  int  _lhsInvert();
  void _dualCalcul();
  int  _prepar();
  void _estimateCalcul(int status);
  void _estimateCalculImage(int status);
  void _estimateCalculXvalidUnique(int status);
  void _simulateCalcul(int status);
  void _neighCalcul(int status, const VectorDouble& tab);
  double _estimateVarZ(int ivarCL, int jvarCL);
  double _variance(int ivarCL, int jvarCL);
  void _variance0();
  void _krigingDump(int status);
  void _simulateDump(int status);
  void _saveWeights(int status);
  void _blockDiscretize();
  bool _isCorrect();

  void   _checkAddress(const String& title,
                       const String& theme,
                       int ival,
                       int nval) const;
  bool   _prepareForImage(const NeighImage* neighI);
  bool   _prepareForImageKriging(Db* dbaux, const NeighImage* neighI);
  int    _bayesPreCalculations();
  void   _bayesPreSimulate();
  void   _bayesCorrectVariance();
  void   _transformGaussianToRaw();
  int    _getFlagAddress(int iech0, int ivar0);
  bool   _isMatCLempty() const;

  void   _setLocalModel(Model* model);
  void   _setInternalShortCutVariablesGeneral();
  void   _setInternalShortCutVariablesModel();
  int    _setInternalShortCutVariablesNeigh();

private:
  // Aggregated classes
  Db*                  _dbin;
  Db*                  _dbout;
  Model*               _modelInit; // Copy of the input model
  ANeigh*              _neigh;
  const AAnam*         _anam;
  bool                 _isReady;

  // Pointer to the Model currently used (must not be freed)
  Model*               _model;
  bool                 _optimEnabled;

  // Calculation modes
  CovCalcMode          _calcModeLHS;
  CovCalcMode          _calcModeRHS;
  CovCalcMode          _calcModeVAR;

  // Options

  /// UID for storage
  int  _iptrEst;
  int  _iptrStd;
  int  _iptrVarZ;
  bool _flagEst;
  bool _flagStd;
  bool _flagVarZ;
  bool _flagGlobal;
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

  /// Options complement for neighborhood
  bool _flagCode;  // True when kriging by Code (Profile)

  /// Option for Block estimation
  bool _flagPerCell;
  int  _ndiscNumber;
  VectorInt _ndiscs;
  VectorVectorDouble _disc1; // Dimension: ndiscNumber, ndim
  VectorVectorDouble _disc2; // Dimension: ndiscNumber, ndim

  /// Option for Cross_validation
  bool _xvalidEstim;
  bool _xvalidStdev;
  bool _xvalidVarZ;

  /// Option for Colocation
  VectorInt _rankColCok;

  /// Option for Bayesian
  bool _flagBayes;
  int  _seedForBayes;
  VectorDouble _priorMean; // Dimension NF
  VectorDouble _priorCov;  // Dimension NF * NF
  VectorDouble _postMean;
  VectorDouble _postCov;
  VectorVectorDouble _postSimu;
  VectorDouble _varCorrec;
  Model* _modelSimple;

  /// Option for Discrete Gaussian Model
  bool   _flagDGM;

  /// Option for (Disjunctive) Kriging of Factor
  bool _flagFactorKriging;
  int _nclasses;

  /// Option for Estimating the Linear Combination of Variables
  VectorVectorDouble _matCL;

  /// Option for asking for Z * A-1 * Z
  bool   _flagLTerm;
  double _lterm;

  /// Option for Gaussian Kriging
  bool   _flagAnam;

  /// Option for Estimation based on Image
  int     _seedForImage;
  DbGrid* _dbaux;

  /// Option for saving the Weights using Keypair mechanism
  bool _flagKeypairWeights;

  /// Option for Neighboring test
  bool _flagNeighOnly;
  int  _iptrNeigh;

  /// Local variables
  int _iechOut;
  int _ndim;
  int _nvar;
  int _nvarCL;
  int _nech;
  int _nbfl;
  int _nfeq;
  int _nfex;
  int _neq;
  int _nred;
  bool _flagIsotopic;

  /// Working arrays
  mutable bool _flagCheckAddress;
  mutable VectorInt    _nbgh;
  mutable VectorInt    _flag;
  mutable MatrixSquareGeneral _covtab;
  mutable MatrixSquareGeneral _covref;
  mutable VectorDouble _drftab;
  mutable VectorDouble _lhs;
  mutable VectorDouble _lhsinv;
  mutable VectorDouble _rhs;
  mutable VectorDouble _wgt;
  mutable VectorDouble _zam;
  mutable VectorDouble _zext;
  mutable MatrixSquareGeneral _var0;
  mutable VectorInt    _dbinUidToBeDeleted;
  mutable VectorInt    _dboutUidToBeDeleted;

  /// Some Space Point allocated once for all
  mutable SpaceRN    _space;
  mutable SpacePoint _p0;
  mutable SpacePoint _p1;
  mutable SpacePoint _p2;
  mutable SpacePoint _p0_memo;
};
