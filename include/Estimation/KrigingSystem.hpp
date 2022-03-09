/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Neigh/NeighWork.hpp"
#include "Basic/Vector.hpp"
#include "Enum/EKrigOpt.hpp"

class Db;
class DbGrid;
class Model;
class ANeighParam;
class CovCalcMode;
class ECalcMember;
class NeighImage;

class GSTLEARN_EXPORT KrigingSystem
{
public:
  KrigingSystem(Db* dbin,
                Db* dbout,
                Model* model,
                ANeighParam* neighParam);
  KrigingSystem(const KrigingSystem &m) = delete;
  KrigingSystem& operator=(const KrigingSystem &m) = delete;
  virtual ~KrigingSystem();

  int  setKrigOptEstim(int iptrEst, int iptrStd, int iptrVarZ);
  int  setKrigOptCalcul(const EKrigOpt& calcul, const VectorInt& ndiscs);
  int  setKrigOptXValid(bool optionXValidEstim = false,
                        bool optionXValidStdev = false);
  int  setKrigOptColCok(const VectorInt& rank_colcok);
  int  setKrigOptBayes(bool flag_bayes);
  int  setKrigOptMatCL(const VectorVectorDouble& matCL);
  int  setKrigoptCode(bool flag_code);
  void setKrigOptFlagSimu(bool flagSimu);
  void setKrigOptSaveWeights(bool flag_save);
  int  setKrigOptImageSmooth(bool flag_smooth, int type = 1, double range = 0.);
  void setKrigOptCheckAddress(bool flagCheckAddress);

  bool isReady();
  int  estimate(int iech_out);

private:
  int  _getNVar() const;
  int  _getNVarCL() const;
  int  _getNbfl() const;
  int  _getNFeq() const;
  int  _getNech() const;
  int  _getNDim() const;
  int  _getNFex() const;
  int  _getNeq()  const;
  int  _getNDisc() const;
  void _setFlag(int iech, int ivar, int value);
  int  _getFlag(int iech, int ivar);
  double _getIdim(int loc_rank, int idim, int iech_out = 0) const;
  double _getFext(int rank, int ibfl, int iech_out = 0) const;
  double _getIvar(int rank, int ivar, int iech_out = 0) const;
  double _getVerr(int rank, int ivar, int iech_out = 0) const;
  double _getMean(int ivarCL);
  void _resetMemoryPerNeigh();
  void _resetMemoryGeneral();
  void _flagDefine();
  void _covtabInit(bool flag_init);
  void _covtabCalcul(bool flag_init,
                     const CovCalcMode& mode,
                     int iech1,
                     int iech2,
                     VectorDouble d1);
  void _drftabCalcul(const ECalcMember &member, const Db* db, int iech);
  bool _isAuthorized();
  double _continuousMultiplier(int rank1,int rank2, double eps = EPSILON4);
  void _lhsCalcul();
  void _lhsIsoToHetero();
  void _lhsDump(int nbypas = 5);
  int  _rhsCalcul(int rankRandom = -1);
  void _rhsIsoToHetero();
  void _rhsDump();
  void _wgtCalcul();
  void _wgtDump(int status);
  VectorInt _getRelativePosition();
  int  _lhsInvert();
  void _dual(bool flagLterm, double *lterm);
  int  _prepar();
  void _estimateCalcul(int status);
  void _estimateCalculImage(int status);
  void _estimateCalculSmoothImage(int status);
  void _estimateCalculXvalidUnique(int status);
  double _estimateVarZ(int ivarCL, int jvarCL);
  double _variance(int ivarCL, int jvarCL, const double* varb = nullptr);
  void _variance0();
  void _krigingDump(int status);
  void _saveWeights(int status);
  void _blockDiscretize();
  bool _isCorrect();

  int    _IND(int iech, int ivar,int nech);
  int    _getFLAG(int iech,int ivar);
  double _getCOVTAB(int ivar,int jvar);
  void   _setCOVTAB(int ivar,int jvar,double value);
  void   _addCOVTAB(int ivar,int jvar,double value);
  void   _prodCOVTAB(int ivar,int jvar,double value);
  double _getRHS(int iech, int ivar, int jvCL);
  void   _setRHS(int iech, int ivar, int jvCL, double value, bool isForDrift = false);
  double _getRHSC(int i, int jvCL);
  double _getWGTC(int i,int jvCL);
  double _getLHS(int iech, int ivar, int jech, int jvar);
  double _getLHSINV(int iech, int ivar, int jech, int jvar);
  void   _setLHS(int iech, int ivar, int jech, int jvar, double value, bool isForDrift = false);
  void   _addLHS(int iech, int ivar, int jech, int jvar, double value);
  void   _prodLHS(int iech, int ivar, int jech, int jvar, double value);
  double _getDISC1(int idisc, int idim);
  void   _setDISC1(int idisc, int idim, double value);
  double _getDISC2(int idisc,int idim);
  void   _setDISC2(int idisc,int idim, double value);
  double _getVAR0(int ivCL, int jvCL);
  void   _setVAR0(int ivCL, int jvCL, double value);
  void   _checkAddress(const String& title,const String& theme,int ival,int nval);
  bool   _prepareForImage(const NeighImage* neighI);
  bool   _prepareForImageKriging(Db* dbaux);

private:
  // Aggregated classes
  Db*                  _dbin;
  Db*                  _dbout;
  Model*               _model;
  ANeighParam*         _neighParam;
  bool _isReady;

  // Options

  /// UID for storage
  int  _iptrEst;
  int  _iptrStd;
  int  _iptrVarZ;
  bool _flagEst;
  bool _flagStd;
  bool _flagVarZ;
  bool _flagSimu;

  /// Option for Calculation
  EKrigOpt _calcul;

  /// Options complement for neighborhood
  bool _flagCode;  // True when kriging by Code

  /// Option for Block estimation
  int _discreteMode;  // 1 : constant; 2 : per Target cell
  int _ndiscNumber;
  VectorInt    _ndiscs;
  VectorDouble _disc1;
  VectorDouble _disc2;

  /// Option for Cross_validation
  bool _xvalidEstim;
  bool _xvalidStdev;

  /// Option for Colocation
  VectorInt _rankColCok;

  /// Option for Bayesian
  bool _flagBayes;
  VectorDouble _rmean;

  /// Option for Disjunctive Kriging
  bool _flagDGM;
  double _supportCoeff;

  // Option for Estimating the Linear Combination of Variables
  VectorVectorDouble _matCL;

  // Option for Estimation based on Image
  DbGrid* _dbaux;
  bool _flagSmooth;
  int  _smoothType;
  double _smoothRange;

  // Option for saving the Weights using Keypair mechanism
  bool _flagSaveWeights;

  // Local variables
  int _iechOut;
  int _nred;

  // Local arrays
  mutable bool _flagCheckAddress;
  mutable NeighWork    _nbghWork;
  mutable VectorInt    _nbgh;
  mutable VectorInt    _flag;
  mutable VectorDouble _covtab;
  mutable VectorDouble _drftab;
  mutable VectorDouble _lhs;
  mutable VectorDouble _lhsinv;
  mutable VectorDouble _rhs;
  mutable VectorDouble _wgt;
  mutable VectorDouble _zam;
  mutable VectorDouble _var0;
  mutable VectorInt    _dbinUidToBeDeleted;
  mutable VectorInt    _dboutUidToBeDeleted;
};
