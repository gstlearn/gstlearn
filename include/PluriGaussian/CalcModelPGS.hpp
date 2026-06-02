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

#include "Basic/VectorNumT.hpp"
#include "Calculators/ACalculator.hpp"
#include "Db/Db.hpp"
#include "LithoRule/PropDef.hpp"
#include "LithoRule/Rule.hpp"
#include "LithoRule/RuleProp.hpp"
#include "Model/Model.hpp"
#include "PluriGaussian/CorPGS.hpp"
#include "PluriGaussian/DiscretePGS.hpp"
#include "PluriGaussian/TracePGS.hpp"
#include "Variogram/Vario.hpp"
#include "Variogram/VarioOrder.hpp"
#include "Variogram/VarioParam.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
  class PileRelem;

  class GSTLEARN_EXPORT CalcModelPGS: public ACalculator
  {
  public:
    CalcModelPGS(
      Db* db = nullptr,
      const VarioParam* varioparam = nullptr,
      const RuleProp* ruleprop = nullptr);
    CalcModelPGS(const CalcModelPGS& r) = delete;
    CalcModelPGS& operator=(const CalcModelPGS& r) = delete;
    virtual ~CalcModelPGS();

    bool _check() override;
    bool _preprocess() override;
    bool _run() override;

    Id getMemInt(Id ipair) const { return _memint[ipair]; }

    Id getIFirst() const { return _ifirst; }

    Id getILast() const { return _ilast; }

    Id getIpasCur() const { return _ipascur; }

    Id getIdirCur() const { return _idircur; }

    Id getIgrfCur() const { return _igrfcur; }

    Id getNFacies() const { return _nfacies; }

    bool getFlagStat() const { return _flagStat; }

    bool getUseDiscrete() const { return _useDiscrete; }

    const Db* getDb() const { return _db; }

    const DiscretePGS* getDiscretePGS() const { return _discretepgs; }

    double getStatThresh(Id ifac, Id igrf, Id rank)
    {
      return _statThresh[2 * (_nfacies * (igrf) + (ifac)) + (rank)];
    }

    VectorDouble& getStatProba() { return _statProba; }

    double getStatProba(Id i, Id j) { return _statProba[_nfacies * i + j]; }

    const VarioOrder& getVarioOrder() const { return _vorder; }

    Vario* getVario() const { return _vario; }

    Vario* getVarioInd() const { return _varioind; }

    const Rule* getRule() const { return _rule; }

    void setMemInt(Id ipair, double memint) { _memint[ipair] = memint; }

    void setFlagRho(bool flagRho) { _flagRho = flagRho; }

    void setStatThresh(Id ifac, Id igrf, Id rank, double thresh)
    {
      _statThresh[2 * (_nfacies * (igrf) + (ifac)) + (rank)] = thresh;
    }

    void setRunType(Id runType) { _runType = runType; }

    void setOptCorrel(Id optCorrel) { _optCorrel = optCorrel; }

    void setStatProba(Id i, Id j, double proba)
    {
      _statProba[_nfacies * i + j] = proba;
    }

    void setUseDb(bool useDb) { _useDb = useDb; }

    void setUseDiscrete(bool useDiscrete) { _useDiscrete = useDiscrete; }

    void setRho(double rho);

    void setModel1(const Model* model1) { _model1 = model1; }

    void setModel2(const Model* model2) { _model2 = model2; }

    void setNgrfMax(Id ngrfmax) { _ngrfmax = ngrfmax; }

    double getProbaInd(
      double correl,
      double low[2],
      double up[2],
      Id iconf,
      Id maxpts = 8000);

    double varcalcCorrelatedGRF(Id idir);
    void varcalcUncorrelatedGRF(Id idir);

    std::vector<Rule> getSortedRules(double eps = EPSILON6) const;

  private:
    double _ruleCalcul(const VectorInt& string);
    bool _checkTestDiscret() const;
    VectorDouble _evaluate(
      PileRelem* relem,
      VectorInt& fgrf,
      VectorInt& fcmp,
      Id* nscore,
      Id* r_opt,
      VectorDouble& scores);
    void _variogramDefineVars();
    void _setBounds(
      bool flag_one,
      Id ifac,
      Id iech,
      double t1min,
      double t1max,
      double t2min,
      double t2max);
    Id _calculateThreshStat();
    void _copySwhh();
    void _variogramPatchC00(Id idir, double rho);
    Id _varcalcFromVarioStat();
    Id _varioPgsVariable(Id mode, Id flag_one, Id flag_prop);
    Id _varioPgsStat();
    Id _varioPgsNoStat();
    Id _discardPoint(Id iech);
    Id _variogramGeometryPgsCalcul(Id idir);
    void _variogramGeometryPgsCorrect(Id idir);
    Id _variogramGeometryPgsFinal();
    bool _isLagUnused(Id idir, Id ilag);
    Id _getCount(Id ifac1, Id ifac2);
    double _optimOneLagPgs(
      Id new_val,
      double maxiter = 100,
      double delta0 = 1,
      double tolsort = EPSILON3);
    bool _preparePropdefAndRule();
    Id _variopgsCalculNoRho(Id flag_geometry);
    Id _variopgsCalculRho();
    void _makeSomeLagsInactive();
    void _makeAllLagsActive();
    void _defineBounds(
      Id iech1,
      Id iech2,
      Id ifac1,
      Id ifac2,
      double* low,
      double* up,
      double* ploc);
    double _calculGlobalScore(
      bool flag_deriv,
      bool flag_reset,
      const VectorDouble& params,
      VectorDouble& Grad,
      MatrixSymmetric& Hess,
      MatrixSymmetric& JJ);
    Id _varioIndicStat();
    Id _varioIndicNoStat();
    Id _ruleAuto();
    double _getProba(
      bool flagIndependent,
      double* low,
      double* up,
      const VectorInt& iconf,
      double* cov,
      Id maxpts = 8000);
    double _getValue(
      bool flagIndependent,
      Id iech1,
      Id iech2,
      Id ifac1,
      Id ifac2,
      const VectorInt& iconf,
      double* cov);
    bool _calculCovmatrix(
      VectorDouble& d0,
      const VectorDouble& d1,
      VectorInt& iconf,
      double* cov);
    void _updateVarianceStat();
    void _updateVarianceNoStat();
    double _calculGlobalScoreStat(
      Id flag_deriv,
      Id flag_reset,
      MatrixSymmetric& correl,
      VectorDouble& Grad,
      MatrixSymmetric& Hess,
      MatrixSymmetric& JJ,
      Id maxpts2 = 4000,
      Id maxpts4 = 10000);

    double _calculGlobalScoreNoStat(
      Id flag_deriv,
      Id flag_reset,
      MatrixSymmetric& correl,
      VectorDouble& Grad,
      MatrixSymmetric& Hess,
      MatrixSymmetric& JJ,
      Id maxpts2 = 4000,
      Id maxpts4 = 10000);

    static double _st_rkl(
      double x,
      double y,
      VectorDouble& lower,
      VectorDouble& upper,
      MatrixSymmetric& correl,
      MatrixSquare& covar,
      MatrixSquare& temp,
      Id maxpts);
    static double _st_ikl(
      Id index1,
      Id index2,
      VectorDouble& lower,
      VectorDouble& upper,
      MatrixSymmetric& correl,
      Id maxpts);
    static double _st_d2_dkldij(
      VectorDouble& lower,
      VectorDouble& upper,
      MatrixSymmetric& correl);
    static double _st_d2_dkldkl(
      Id index1,
      Id index2,
      VectorDouble& lower,
      VectorDouble& upper,
      MatrixSymmetric& correl,
      Id maxpts);
    static double _st_d2_dkldkj(
      Id index1,
      Id index2,
      VectorDouble& lower,
      VectorDouble& upper,
      MatrixSymmetric& correl);
    static double _st_nkl(
      VectorDouble& u,
      double lower,
      double upper,
      VectorDouble& invvari,
      Id index2,
      double meanj,
      double varj,
      double stdj);

  private:
    Db* _db; // External pointer: not to be deleted
    const VarioParam* _varioParam; // External pointer: not to be deleted
    const RuleProp* _ruleProp; // External pointer: not to be deleted
    const Model* _model1;
    const Model* _model2;

    mutable const Rule* _rule; // Internal pointer
    const Db* _dbprop; // Internal pointer
    PropDef _propdef;

    Model* _model;
    Vario* _vario; // Used for variogram_pgs calculation
    Vario* _varioind;
    DiscretePGS* _discretepgs;

    Id _runType; // 1 for variopgs; 2 for modelpgs; 3 for rulepgs
    bool _useDb;
    bool _useDiscrete;
    bool _flagStat;
    bool _flagFacies;
    bool _flagDist;
    bool _flagRho;
    Id _optCorrel;
    Id _ngrfmax;

    Id _igrfcur;
    Id _idircur;
    Id _ipascur;
    Id _ndim;
    Id _ngrf;
    Id _nfacies;
    Id _ifirst;
    Id _ilast;

    VectorDouble _props;
    VectorDouble _memint;
    VectorDouble _statProba;
    VectorDouble _statThresh;
    VectorDouble _ruleAutoScores;
    VectorVectorInt _ruleAutoRules;

    CorPGS _corpgs;
    TracePGS _tracepgs;
    VarioOrder _vorder;
  };

  GSTLEARN_EXPORT Id getRuleAutoNRULE();
  GSTLEARN_EXPORT Id getRuleAutoNCOLOR();
  GSTLEARN_EXPORT Id getRuleAutoNGRF();
  GSTLEARN_EXPORT void setRuleAutoNRULE(Id nrule);
  GSTLEARN_EXPORT void setRuleAutoNCOLOR(Id ncolor);
  GSTLEARN_EXPORT void setRuleAutoNGRF(Id ngrf);

  GSTLEARN_EXPORT Vario* variogram_pgs(
    Db* db,
    const VarioParam* varioparam,
    const RuleProp* ruleprop,
    bool use_discrete = false,
    Id flag_rho = false,
    Id opt_correl = 2);

  GSTLEARN_EXPORT Vario* model_pgs(
    Db* db,
    const VarioParam* varioparam,
    const RuleProp* ruleprop,
    const Model* model1,
    const Model* model2 = nullptr,
    bool use_discrete = false);

} // namespace gstlrn
