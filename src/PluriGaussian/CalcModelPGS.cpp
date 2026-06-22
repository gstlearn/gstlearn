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
#include "PluriGaussian/CalcModelPGS.hpp"

#include "Basic/Law.hpp"
#include "Basic/MathFunc.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/VectorHelper.hpp"
#include "Core/Keypair.hpp"
#include "LithoRule/RuleProp.hpp"
#include "LithoRule/RuleShift.hpp"
#include "Matrix/EigenVectors.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "PluriGaussian/CorPGS.hpp"
#include "PluriGaussian/PileRelem.hpp"
#include "PluriGaussian/PileSplit.hpp"
#include "PluriGaussian/TracePGS.hpp"
#include "Stats/Classical.hpp"
#include "Variogram/VarioOrder.hpp"
#include "geoslib_old_f.h"
#include <boost/mpl/bool_fwd.hpp>

namespace gstlrn
{
#define COVS(ivar, jvar) (covs[(ivar) * _nfacies + (jvar)])

  thread_local Id NRULE = 0;
  thread_local Id NCOLOR = 0;
  thread_local Id NGRF = 0;

  static double _epsGoldenSearch = 5 * EPSILON2;
  static double _releps = 0.;
  static double _abseps = EPSILON6;
  static double _deltaparam = EPSILON6;
  static Id _maxpts4 = 10000;

  Id getRuleAutoNRULE()
  {
    return NRULE;
  }

  Id getRuleAutoNCOLOR()
  {
    return NCOLOR;
  }

  Id getRuleAutoNGRF()
  {
    return NGRF;
  }

  void setRuleAutoNRULE(Id nrule)
  {
    NRULE = nrule;
  }

  void setRuleAutoNCOLOR(Id ncolor)
  {
    NCOLOR = ncolor;
  }

  void setRuleAutoNGRF(Id ngrf)
  {
    NGRF = ngrf;
  }

  static double _goldenSearchFuncStat(double correl, void* user_data)
  {
    auto* la = static_cast<CalcModelPGS*>(user_data);

    Id ilag = la->getIpasCur();
    Id idir = la->getIdirCur();
    Id igrf = la->getIgrfCur();
    Id nfacies = la->getNFacies();
    const auto* varioind = la->getVarioInd();

    Id iconf0 = 0;
    if (la->getUseDiscrete())
    {
      double cround;
      iconf0 = la->getDiscretePGS()->getCovrank(correl, &cround);
    }

    double sum = 0.;
    double low[2], up[2];
    for (Id ifac1 = 0; ifac1 < nfacies; ifac1++)
      for (Id ifac2 = 0; ifac2 < nfacies; ifac2++)
      {
        low[0] = la->getStatThresh(ifac1, igrf, 0);
        up[0] = la->getStatThresh(ifac1, igrf, 1);
        low[1] = la->getStatThresh(ifac2, igrf, 0);
        up[1] = la->getStatThresh(ifac2, igrf, 1);
        double proba = la->getProbaInd(correl, low, up, iconf0);

        double logp = (proba <= 0.) ? MINIMUM_BIG : log(proba);
        Id iad = varioind->getAddressForGg(idir, ifac1, ifac2, ilag, 1);
        double sw = varioind->getSwByIndex(idir, iad);
        double gg = varioind->getGgByIndex(idir, iad);
        iad = varioind->getAddressForGg(idir, ifac1, ifac2, ilag, -1);
        gg += varioind->getGgByIndex(idir, iad);
        sum -= logp * gg * sw / 2.;
      }

    return (0.5 * sum);
  }

  static double _goldenSearchFuncNoStat(double correl, void* user_data)
  {
    auto* la = static_cast<CalcModelPGS*>(user_data);
    const auto* db = la->getDb();

    Id iconf0 = 0;
    if (la->getUseDiscrete())
    {
      double cround;
      iconf0 = la->getDiscretePGS()->getCovrank(correl, &cround);
    }

    /* Reset the pre-calculation array (only if flag_stat) */
    if (la->getFlagStat())
      for (Id i = 0; i < static_cast<Id>(la->getStatProba().size()); i++)
        la->getStatProba()[i] = TEST;

    double sum = 0.;
    for (Id ipair = la->getIFirst(); ipair < la->getILast(); ipair++)
    {
      Id i1;
      Id i2;
      double dist;
      la->getVarioOrder().getIndices(ipair, &i1, &i2, &dist);
      double w1 = db->getWeight(i1);
      double w2 = db->getWeight(i2);
      Id ifac1 = -1;
      Id ifac2 = -1;
      double proba = TEST;
      if (la->getFlagStat())
      {

        /* In the stationary case, search in the lookup table first */

        ifac1 = static_cast<Id>(db->getZVariable(i1, 0)) - 1;
        if (ifac1 < 0 || ifac1 >= la->getNFacies()) continue;
        ifac2 = static_cast<Id>(db->getZVariable(i2, 0)) - 1;
        if (ifac2 < 0 || ifac2 >= la->getNFacies()) continue;
        proba = la->getStatProba(ifac1, ifac2);
      }

      if (FFFF(proba))
      {
        double low[2], up[2];
        if (!la->getUseDiscrete())
        {
          low[0] = db->getLocVariable(ELoc::L, i1, la->getIgrfCur());
          up[0] = db->getLocVariable(ELoc::U, i1, la->getIgrfCur());
          low[1] = db->getLocVariable(ELoc::L, i2, la->getIgrfCur());
          up[1] = db->getLocVariable(ELoc::U, i2, la->getIgrfCur());
        }
        else
        {
          low[0] = db->getLocVariable(ELoc::RKLOW, i1, la->getIgrfCur());
          up[0] = db->getLocVariable(ELoc::RKUP, i1, la->getIgrfCur());
          low[1] = db->getLocVariable(ELoc::RKLOW, i2, la->getIgrfCur());
          up[1] = db->getLocVariable(ELoc::RKUP, i2, la->getIgrfCur());
        }

        proba = la->getProbaInd(correl, low, up, iconf0);

        if (la->getFlagStat())
        {
          la->setStatProba(ifac1, ifac2, proba);
          la->setStatProba(ifac2, ifac1, proba);
        }
      }
      double logp = (proba <= 0.) ? MINIMUM_BIG : log(proba);
      sum -= w1 * w2 * logp;
    }
    return (0.5 * sum);
  }

  static double _geoldenSearchFuncRho(double rho, void* user_data)
  {
    auto* la = static_cast<CalcModelPGS*>(user_data);
    double sum = 0;
    Id ndir = la->getVario()->getNDir();
    la->setRho(rho);

    /* Evaluation of the global cost-function */

    for (Id idir = 0; idir < ndir; idir++)
      sum += la->varcalcCorrelatedGRF(idir);

    if (OptDbg::query(EDbg::CONVERGE))
      message(
        "Value of the evaluating function = %lf - rho value %lf\n", sum, rho);
    return (sum);
  }

  CalcModelPGS::CalcModelPGS(
    Db* db,
    const VarioParam* varioparam,
    const RuleProp* ruleprop)
    : ACalculator()
    , _db(db)
    , _varioParam(varioparam)
    , _ruleProp(ruleprop)
    , _model1(nullptr)
    , _model2(nullptr)
    , _rule(nullptr)
    , _dbprop(nullptr)
    , _propdef()
    , _model(nullptr)
    , _vario(nullptr)
    , _varioind(nullptr)
    , _discretepgs(nullptr)
    , _runType(0)
    , _useDb(false)
    , _useDiscrete(false)
    , _flagStat(true)
    , _flagFacies(false)
    , _flagDist(false)
    , _flagRho(false)
    , _optCorrel(0)
    , _ngrfmax(0)
    , _igrfcur(0)
    , _idircur(0)
    , _ipascur(0)
    , _ndim(0)
    , _ngrf(0)
    , _nfacies(0)
    , _ifirst(0)
    , _ilast(0)
    , _props()
    , _memint()
    , _statProba()
    , _statThresh()
    , _corpgs()
    , _tracepgs()
    , _vorder(_flagDist, 0)
  {
  }

  CalcModelPGS::~CalcModelPGS()
  {
    delete _model;
    delete _varioind;
    delete _discretepgs;
  }

  bool CalcModelPGS::_check()
  {
    // Structures that msut always be defined
    if (_useDb && _db == nullptr)
    {
      messerr("The Db must be provided");
      return false;
    }
    if (_varioParam == nullptr)
    {
      messerr("The VarioParam must be provided");
      return false;
    }

    // Extracting information from RuleProp class (compulsory)
    _flagStat = _ruleProp->isFlagStat();
    _dbprop = _ruleProp->getDbprop();

    // Extract information from Rule (obviously not defined for ruleAuto)
    if (_runType != 3)
    {
      _rule = _ruleProp->getRule();
      if (_rule == nullptr)
      {
        messerr("The Rule must be defined in the RuleProp");
        return false;
      }
      _ngrf = _rule->getNGRF();
      if (_rule->getModeRule() == ERule::SHIFT) _ngrf++;
    }
    else
    {
      _ngrf = _ngrfmax;
    }
    if (_ngrf <= 0)
    {
      messerr("The number of GRF is not valid");
      return false;
    }

    // Check the number of Facies
    _nfacies = _ruleProp->getPropCst().size();
    if (_db != nullptr)
    {
      if (_nfacies > 0 && _nfacies != _db->getNFacies())
      {
        messerr(
          "Number of Facies in 'propcst' (%d) should match Number of Facies "
          "in 'db' (%d)",
          _nfacies, _db->getNFacies());
        return false;
      }
      _nfacies = _db->getNFacies();
    }

    if (_rule != nullptr)
    {
      if (_nfacies > 0 && _nfacies != _rule->getNFacies())
      {
        messerr(
          "Number of Facies in 'propcst' (%d) should match Number of Facies "
          "in 'rule' (%d)",
          _nfacies, _rule->getNFacies());
        return false;
      }
      _nfacies = _rule->getNFacies();
    }
    if (_nfacies <= 0)
    {
      messerr("No Facies class has been found");
      return false;
    }

    // Dimension _statThresh and _statProba
    if (_flagStat)
    {
      _statProba.resize(_nfacies * _nfacies, 0.);
      _statThresh.resize(
        _nfacies * 2 * 2, 0.); // Do not use ngrf, use 2 instead
    }
    else
    {
      if (_db == nullptr)
      {
        messerr("The Db must be provided in the non-stationary case");
        return false;
      }
    }

    // Check based on Db structure
    if (_db != nullptr)
    {
      if (_db->getNLoc(ELoc::Z) != 1)
      {
        messerr(
          "The number of variables (%d) must be equal to 1",
          _db->getNLoc(ELoc::Z));
        return false;
      }
    }

    // Checks based on VarioParam structure
    if (_varioParam != nullptr)
    {
      if (_varioParam->getNDir() <= 0)
      {
        messerr(
          "The variogram must contain at least one calculation Direction");
        return false;
      }
      if (_varioParam->getNDir() <= 0)
      {
        messerr(
          "The variogram must contain at least one calculation Direction");
        return false;
      }
    }

    // Create the combined model (only when _runType == 2)
    if (_runType == 2)
    {
      _model = model_rule_combine(_model1, _model2, _rule);
      if (_model == nullptr)
      {
        messerr("The Model(s) must be defined");
        return false;
      }
      if (_model->getNVar() != _ngrf)
      {
        messerr("The number of GRF is not equal to the number of variables");
        messerr("defined in the combined Model");
        return false;
      }
    }

    // Check on Space dimension
    _ndim = 0;
    if (_db != nullptr)
    {
      if (_ndim > 0 && _ndim != _db->getNDim())
      {
        messerr("Inconsistent Space Dimension:");
        messerr("- Current value = %d", _ndim);
        messerr("- Db value = %d", _db->getNDim());
        return false;
      }
      _ndim = _db->getNDim();
    }

    if (_varioParam != nullptr)
    {
      if (_ndim > 0 && _ndim != _varioParam->getNDim())
      {
        messerr("Inconsistent Space Dimension:");
        messerr("- Current value = %d", _ndim);
        messerr("- Variogram value = %d", _varioParam->getNDim());
        return false;
      }
      _ndim = _varioParam->getNDim();
    }
    if (_dbprop != nullptr)
    {
      if (_ndim > 0 && _ndim != _dbprop->getNDim())
      {
        messerr("Inconsistent Space Dimension:");
        messerr("- Current value = %d", _ndim);
        messerr("- Dbprop value = %d", _dbprop->getNDim());
        return false;
      }
      _ndim = _dbprop->getNDim();
    }
    if (_model != nullptr)
    {
      if (_ndim > 0 && _ndim != static_cast<Id>(_model->getNDim()))
      {
        messerr("Inconsistent Space Dimension:");
        messerr("- Current value = %d", _ndim);
        messerr("- Model value = %d", _model->getNDim());
        return false;
      }
      _ndim = _model->getNDim();
    }

    if (!_checkTestDiscret()) return false;

    if (_runType != 2 && _rule != nullptr && _rule->getModeRule() != ERule::STD)
    {
      messerr("This function is only programmed for standard rule");
      return (1);
    }

    return true;
  }

  bool CalcModelPGS::_preprocess()
  {
    Id nvar = (_runType == 1) ? _ngrf : _nfacies;
    const ECalcVario calculType = (_runType == 1 || _runType == 3)
                                  ? ECalcVario::COVARIANCE_NC
                                  : ECalcVario::VARIOGRAM;

    _props.clear();
    if (_flagStat)
    {
      const VectorDouble& propcst = _ruleProp->getPropCst();
      if (!propcst.empty())
      {
        if (static_cast<Id>(propcst.size()) != _nfacies)
        {
          messerr(
            "Number of proportions in 'propcst' (%d) should match Number of "
            "Facies in 'rule' (%d)",
            static_cast<Id>(propcst.size()), _rule->getNFacies());
          return false;
        }
        _props = propcst;
      }
      else
      {
        // Calculate the number of Facies in 'Db'
        _props = dbStatisticsFacies(_db);
        if (static_cast<Id>(_props.size()) != _nfacies)
        {
          messerr(
            "Number of Facies in 'db' (%d) should match Number of facies in "
            "'rule' (%d)",
            static_cast<Id>(_props.size()), _rule->getNFacies());
          return false;
        }
      }

      // Calculate the variogram of Indicators
      _varioind = new Vario(*_varioParam);
      if (_varioind->computeIndic(_db, calculType)) return false;
    }

    /* Pre-calculation of integrals: Define the structure */

    if (_useDiscrete) _discretepgs = new DiscretePGS(1, 200, 100, -1., 1.);

    // Initiate the output variogram

    _vario = Vario::create(*_varioParam);
    _vario->setDb(_db);
    _vario->setNVar(nvar);
    return _vario->prepare(calculType) == 0;
  }

  bool CalcModelPGS::_run()
  {

    // Dispatch

    switch (_runType)
    {
      case 1:
        if (!_preparePropdefAndRule()) return 1;
        if (_flagStat)
        {
          if (_varioPgsStat()) return false;
        }
        else
        {
          if (_varioPgsNoStat()) return false;
        }
        break;

      case 2:
        if (!_preparePropdefAndRule()) return 1;
        if (_flagStat)
        {
          if (_varioIndicStat()) return false;
        }
        else
        {
          if (_varioIndicNoStat()) return false;
        }
        break;

      case 3:

        if (_ruleAuto()) return false;
        break;

      default: messerr("Unknown run type: %d", _runType); return false;
    }
    return true;
  }

  double CalcModelPGS::_ruleCalcul(const VectorInt& string)
  {
    _rule = new Rule(string);
    _ngrf = _rule->getNGRF();
    _vario->setNVar(_ngrf);
    _vario->internalVariableResize();
    _vario->internalMemoryResize();
    _tracepgs.define(_ngrf, _corpgs.getNpar(), false, false);

    if (_flagStat)
    {
      setRho(0.);
      (void)_calculateThreshStat();
      _varcalcFromVarioStat();
    }
    else
    {
      (void)_varioPgsVariable(0, 1, 0);
      setRho(0.);
      for (Id idir = 0, ndir = _vario->getNDir(); idir < ndir; idir++)
      {
        _idircur = idir;
        _variogramPatchC00(idir, _rule->getRho());
        varcalcUncorrelatedGRF(idir);
      }
    }

    /* Deallocation of the Rule */

    delete _rule;
    _rule = nullptr;

    return _tracepgs.extractTrace();
  }

  VectorDouble CalcModelPGS::_evaluate(
    PileRelem* relem,
    VectorInt& fgrf,
    VectorInt& fcmp,
    Id* nscore,
    Id* r_opt,
    VectorDouble& scores)
  {
    bool flag_check = static_cast<Id>(get_keypone("Multi_Score_Check", 0.));
    Id nmax = static_cast<Id>(pow(2., static_cast<double>(NGRF)));
    Id nrule = relem->getNrule();

    /* Core allocation */

    *nscore = nrule;
    scores.resize(nrule);
    *r_opt = 0;
    for (Id ir = 0; ir < nrule; ir++)
    {

      // Check if the same flag has already been found

      Id indice = -1;
      Id igrf_opt = -1;
      for (Id igrf_cas = 1; igrf_cas < nmax && igrf_opt < 0; igrf_cas++)
      {
        indice = relem->sameScore(ir, igrf_cas, fgrf, fcmp);
        if (indice >= 0) igrf_opt = igrf_cas;
      }

      // Set the score

      if (indice >= 0)
      {
        scores[ir] = scores[indice];
      }
      else
      {
        scores[ir] = _ruleCalcul(relem->getRrule(ir));
        _propdef.reset();
      }

      // When Multi_Score_Check, calculate the score even if already defined
      // and compare both results

      if (flag_check && indice >= 0)
      {
        double score_ref = _ruleCalcul(relem->getRrule(ir));
        if (ABS(scores[ir] - score_ref) > 1.e-10 * score_ref)
        {
          messerr("Warning: Difference between score stored and re-evaluated:");
          messerr("- as already stored = %lf", scores[ir]);
          messerr("- as re-evaluated   = %lf", score_ref);
        }
      }

      // Optional printout

      if (getVerbose())
      {
        if (indice < 0)
          st_rule_print2(
            ir, NRULE, relem->getRrules(), relem->getRfipos(), true, indice,
            igrf_opt, scores[ir]);
      }

      // Ranking the Minimum score

      if (scores[ir] < scores[*r_opt]) *r_opt = ir;
    }

    return (scores);
  }

  Id CalcModelPGS::_ruleAuto()
  {
    Id nscore, r_opt;

    _nfacies = _db->getNFacies();
    _ngrf = _ngrfmax;
    Id nrule = 2 * _nfacies - 1;
    Id flag_rho = 0;
    Id opt_correl = 0;
    bool flag_correl = false;
    setRuleAutoNRULE(nrule);
    setRuleAutoNCOLOR(_nfacies);
    setRuleAutoNGRF(_ngrf);

    VectorInt facies = VH::sequence(_nfacies, 1);
    if (!_preparePropdefAndRule()) return 1;

    /* Allocation */

    if (_flagStat)
    {
      _corpgs.define(0, 0, 0);
      _tracepgs.define(_ngrf, _corpgs.getNpar(), false, false);
    }
    else
    {
      _corpgs.define(opt_correl, flag_rho, 0.);
      _tracepgs.define(_ngrf, _corpgs.getNpar(), flag_rho, flag_correl);

      // Prepare the geometry

      for (Id idir = 0, ndir = _vario->getNDir(); idir < ndir; idir++)
      {
        _idircur = idir;
        if (_variogramGeometryPgsCalcul(idir)) return 1;
        _variogramGeometryPgsCorrect(idir);
      }
      if (_variogramGeometryPgsFinal()) return 1;

      // The thresholds are added lately in order to allow calculation of
      // geometry (without checking the threshold interval (not defined yet)
      if (_varioPgsVariable(1, 1, 0)) return 1;
    }

    /* Elaborate the whole tree of possible Lithotype Rules */

    std::unique_ptr<PileRelem> Pile_Relem = PileRelem::alloc();
    Pile_Relem->define(facies);
    Pile_Relem->subdivide(1, 1);
    Pile_Relem->explore(getVerbose() && OptDbg::query(EDbg::CONVERGE));

    // Evaluate all possibilities

    VectorInt fcmp(_nfacies);
    VectorInt fgrf(1 + _ngrf);
    if (getVerbose())
    {
      mestitle(1, "List of Rules and corresponding scores:");
      if (_flagStat)
        message("Stationary case");
      else
        message("Non-stationary case");
      if (_useDiscrete)
        message(" (Discrete Integration)\n");
      else
        message("\n");
    }
    VectorDouble scores;
    _evaluate(Pile_Relem.get(), fgrf, fcmp, &nscore, &r_opt, scores);

    /* Get the resulting optimal Rule */

    if (getVerbose())
    {
      mestitle(1, "Optimal Lithotype Rule:");
      st_rule_print2(
        r_opt, nrule, Pile_Relem->getRrules(), Pile_Relem->getRfipos(), false,
        -1, -1, TEST);
    }

    // Generate the list of rule contents (by increasing scores)
    auto order = VH::orderRanks(scores);
    _ruleAutoRules.clear();
    _ruleAutoScores.clear();
    for (Id i = 0; i < nscore; i++)
    {
      Id j = order[i];
      _ruleAutoRules.push_back(Pile_Relem->getRrule(j));
      _ruleAutoScores.push_back(scores[j]);
    }

    // Clean the geometry (non-stationary case)
    if (!_flagStat) _vorder.clear();
    (void)_varioPgsVariable(-1, 1, 0);
    return 0;
  }

  std::vector<Rule> CalcModelPGS::getSortedRules(double eps) const
  {
    // Get the decoration from the Rule if provided in input
    VectorString facnames;
    VectorInt faccols;
    VectorInt facvalues;
    if (_ruleProp->getNRule() > 0)
    {
      const auto* ruleRef = _ruleProp->getRule();
      facnames = ruleRef->getFaciesNames();
      faccols = ruleRef->getFaciesColors();
      facvalues = ruleRef->getFaciesValues();
    }

    std::vector<Rule> sortedRules;
    double currentScore = TEST;

    for (Id i = 0, n = _ruleAutoRules.size(); i < n; i++)
    {
      double score = _ruleAutoScores[i];

      // Eliminate the rule if its score is too close to the score of the previous one
      if (ABS(score - currentScore) < eps * currentScore) continue;

      currentScore = score;

      auto rule = Rule(_ruleAutoRules[i]);

      rule.setScore(_ruleAutoScores[i]);

      if (_ruleProp->getNRule() > 0)
      {
        rule.setFaciesNames(facnames);
        rule.setFaciesColors(faccols);
        rule.setFaciesValues(facvalues);
      }

      sortedRules.push_back(rule);
    }
    return sortedRules;
  }

  void CalcModelPGS::_variogramDefineVars()
  {
    for (Id igrf = 0; igrf < _ngrf; igrf++)
      for (Id jgrf = 0; jgrf < _ngrf; jgrf++)
      {
        if (igrf == jgrf)
          _vario->setVar(1., igrf, jgrf);
        else
          _vario->setVar(_rule->getRho(), igrf, jgrf);
      }
  }

  void CalcModelPGS::_setBounds(
    bool flag_one,
    Id ifac,
    Id iech,
    double t1min,
    double t1max,
    double t2min,
    double t2max)
  {
    Id jfac;
    if (!_useDiscrete)
    {
      jfac = (flag_one) ? 0 : ifac;
      _db->setBound(iech, jfac, t1min, t1max);
      if (_ngrf > 1)
      {
        jfac = (flag_one) ? 1 : _nfacies + ifac;
        _db->setBound(iech, jfac, t2min, t2max);
      }
    }
    else
    {
      jfac = (flag_one) ? 0 : ifac;
      _db->setInterval(
        iech, jfac, static_cast<double>(_discretepgs->getRankFromProba(t1min)),
        static_cast<double>(_discretepgs->getRankFromProba(t1max)));
      if (_ngrf > 1)
      {
        jfac = (flag_one) ? 1 : _nfacies + ifac;
        _db->setInterval(
          iech, jfac,
          static_cast<double>(_discretepgs->getRankFromProba(t2min)),
          static_cast<double>(_discretepgs->getRankFromProba(t2max)));
      }
    }
  }

  void CalcModelPGS::setRho(double rho)
  {
    double t1min, t1max, t2min, t2max;

    _corpgs.setRho(rho);
    _rule->setRho(rho);
    double rho2 = rho * rho;
    _variogramDefineVars();

    /* Define the Thresholds */

    if (_flagStat)
    {
      _rule->setProportions(_props);
    }
    else
    {
      for (Id iech = 0; iech < _db->getNSample(); iech++)
      {
        if (!_db->isActive(iech)) continue;
        Id ifac = static_cast<Id>(_db->getZVariable(iech, 0));
        if (_propdef.ruleThreshDefine(
              _db, _rule, ifac, iech, 0, 0, 0, &t1min, &t1max, &t2min, &t2max))
          return;
        _setBounds(true, ifac, iech, t1min, t1max, t2min, t2max);
      }
    }

    /* Update the modif matrix if necessary */

    if (_corpgs.getOptCorrel() == 2)
    {
      _corpgs.setModif(0, 1, rho);
      _corpgs.setModif(0, 2, rho);
      _corpgs.setModif(0, 3, rho2);
      _corpgs.setModif(1, 3, 1 - rho2);
    }
  }

  Id CalcModelPGS::_calculateThreshStat()
  {
    double t1min, t1max, t2min, t2max;

    for (Id ifac = 0; ifac < _nfacies; ifac++)
    {
      if (_propdef.ruleThreshDefine(
            _db, _rule, ifac + 1, 0, 0, 0, 0, &t1min, &t1max, &t2min, &t2max))
        return (1);
      if (!_useDiscrete)
      {
        setStatThresh(ifac, 0, 0, t1min);
        setStatThresh(ifac, 0, 1, t1max);
        setStatThresh(ifac, 1, 0, t2min);
        setStatThresh(ifac, 1, 1, t2max);
      }
      else
      {
        setStatThresh(
          ifac, 0, 0,
          static_cast<double>(_discretepgs->getRankFromProba(t1min)));
        setStatThresh(
          ifac, 0, 1,
          static_cast<double>(_discretepgs->getRankFromProba(t1max)));
        if (_ngrf > 1)
        {
          setStatThresh(
            ifac, 1, 0,
            static_cast<double>(_discretepgs->getRankFromProba(t2min)));
          setStatThresh(
            ifac, 1, 1,
            static_cast<double>(_discretepgs->getRankFromProba(t2max)));
        }
        else
        {
          setStatThresh(
            ifac, 1, 0,
            static_cast<double>(_discretepgs->getRankFromProba(-10.)));
          setStatThresh(
            ifac, 1, 1,
            static_cast<double>(_discretepgs->getRankFromProba(+10.)));
        }
      }
    }

    // Verification
    double total = 0.;
    double correl = _corpgs.getRho();
    double low[2], up[2];

    Id iconf0 = 0;
    if (_useDiscrete)
    {
      double cround;
      iconf0 = _discretepgs->getCovrank(correl, &cround);
    }

    for (Id ifac = 0; ifac < _nfacies; ifac++)
    {
      low[0] = getStatThresh(ifac, 0, 0);
      up[0] = getStatThresh(ifac, 0, 1);
      low[1] = getStatThresh(ifac, 1, 0);
      up[1] = getStatThresh(ifac, 1, 1);
      double proba = getProbaInd(correl, low, up, iconf0);
      total += proba;
    }
    if (ABS(total - 1.) > EPSILON3)
      messerr("The sum of Probabilities (%lf) is not close to 1.", total);
    return (0);
  }

  /****************************************************************************/
  /*!
   **  Duplicate the information from the input variogram to the output variogram
   **  This operation considers (Sw, Hh) and replicates this information for
   **  all directions.
   **  Per direction, the information of the first simple input variogram is
   **  replicated to all simple and cross-variograms outputs
   **
   ** \return  Error return code
   **
   *****************************************************************************/
  void CalcModelPGS::_copySwhh()
  {
    Id nvar = _vario->getNVar();
    for (Id idir = 0, ndir = _vario->getNDir(); idir < ndir; idir++)
    {
      for (Id ilag = 0, nlag = _vario->getNLagTotal(idir); ilag < nlag; ilag++)
      {
        for (Id ivar = 0; ivar < nvar; ivar++)
          for (Id jvar = 0; jvar <= ivar; jvar++)
          {
            Id iadlag = _vario->getAddressForGg(idir, ivar, jvar, ilag);
            _vario->setSwByIndex(
              idir, iadlag, _varioind->getSwByIndex(idir, ilag));
            _vario->setHhByIndex(
              idir, iadlag, _varioind->getHhByIndex(idir, ilag));
          }
      }
    }
  }

  void CalcModelPGS::_variogramPatchC00(Id idir, double rho)
  {
    Id nech = (_db == nullptr) ? 0 : _db->getNSample(true);
    _vario->patchCenter(idir, nech, rho);
  }

  /****************************************************************************/
  /*!
   **  Performing the variogram calculations (stationary case)
   **
   ** \return  Error return code
   **
   *****************************************************************************/
  Id CalcModelPGS::_varcalcFromVarioStat()
  {
    Id iad;
    double result, testval, varloc, niter;

    /* Initializations */

    setRho(0.);

    /* Loop on the directions */

    for (Id idir = 0, ndir = _vario->getNDir(); idir < ndir; idir++)
    {
      _idircur = idir;

      /* Set the value of C(0) */

      _variogramPatchC00(idir, 0.);

      /* Loop on the lags */

      for (Id ilag = 0, nlag = _vario->getNLag(idir); ilag < nlag; ilag++)
      {
        mes_process("Inverting Variogram Lag", _vario->getNLag(idir), ilag);
        _ipascur = ilag;
        _tracepgs.addRow();

        /* Loop on the GRFs */

        for (Id igrf = 0; igrf < _ngrf; igrf++)
        {
          _igrfcur = igrf;
          result = golden_search(
            _goldenSearchFuncStat, static_cast<void*>(this), _epsGoldenSearch,
            -1., 1., &testval, &niter);
          _tracepgs.update(idir + 1, ilag + 1, 2 * igrf, 1, &testval);
          _tracepgs.update(idir + 1, ilag + 1, 2 * igrf + 1, 1, &niter);

          for (Id jgrf = 0; jgrf <= igrf; jgrf++)
          {
            varloc = (igrf == jgrf) ? result : 0.;
            iad = _vario->getAddressForGg(idir, igrf, jgrf, ilag, 1);
            _vario->setGgByIndex(idir, iad, varloc);
            iad = _vario->getAddressForGg(idir, igrf, jgrf, ilag, -1);
            _vario->setGgByIndex(idir, iad, varloc);

            if (OptDbg::query(EDbg::CONVERGE))
              message(
                "Lag:%d - Grf:%d - Variogram(%d) = %lf\n", ilag, igrf, iad,
                varloc);
          }
        }
      }
    }
    return (0);
  }

  /****************************************************************************/
  /*!
   **  Manage local variables for variopgs calculation (Non-stationary case)
   **
   ** \return  Error return code
   **
   ** \param[in]  mode          Type of usage
   **                           1 for allocation
   **                           0 for valuation (rule dependent)
   **                          -1 for deallocation
   ** \param[in]  flag_one      1 for considering only the Facies at data point
   **                           0 for considering all facies
   ** \param[in]  flag_prop     1 for allocating variable for proportions
   **
   *****************************************************************************/
  Id CalcModelPGS::_varioPgsVariable(Id mode, Id flag_one, Id flag_prop)
  {
    Id ifac, jfac, iptr, nloop;
    double t1min, t1max, t2min, t2max;
    static bool is_prop_defined;

    // Dispatch

    if (_db == nullptr) return 0;
    Id number = (flag_one) ? _ngrf : _ngrf * _nfacies;

    switch (mode)
    {
      case 1:

        /* The proportions (if not already present in correct number) */

        is_prop_defined = false;
        if (flag_prop && _db->getNLoc(ELoc::P) != _nfacies)
        {
          iptr = _db->addColumnsByConstant(_nfacies, 0., String(), ELoc::P);
          if (iptr < 0) return (1);
          is_prop_defined = true;
        }

        /* The bounds */

        if (!_useDiscrete)
        {
          iptr = _db->addColumnsByConstant(number, 0., "Lower", ELoc::L);
          if (iptr < 0) return (1);

          iptr = _db->addColumnsByConstant(number, 0., "Upper", ELoc::U);
          if (iptr < 0) return (1);
        }
        else
        {
          iptr =
            _db->addColumnsByConstant(number, 0., "Lower Rank", ELoc::RKLOW);
          if (iptr < 0) return (1);

          iptr =
            _db->addColumnsByConstant(number, 0., "Upper Rank", ELoc::RKUP);
          if (iptr < 0) return (1);
        }
        break;

      case 0:

        /* Evaluate the bounds */
        /* Use dummy rho value in order to avoid discarding pairs in geometry */

        nloop = (flag_one) ? 1 : _nfacies;
        for (Id iech = 0, nech = _db->getNSample(); iech < nech; iech++)
        {
          if (!_db->isActive(iech)) continue;

          for (Id i = 0; i < nloop; i++)
          {
            ifac = (flag_one) ? static_cast<Id>(_db->getZVariable(iech, 0)) : i;
            jfac = (flag_one) ? ifac : ifac + 1;
            if (_propdef.ruleThreshDefine(
                  _db, _rule, jfac, iech, 0, 0, 0, &t1min, &t1max, &t2min,
                  &t2max))
              return (1);

            /* Define the proportions */

            if (flag_prop)
              _db->setLocVariable(ELoc::P, iech, ifac, _propdef._propmem[ifac]);

            /* Define the bounds */

            _setBounds(flag_one, ifac, iech, t1min, t1max, t2min, t2max);
          }
        }
        break;

      case -1:

        /* Deallocation */

        if (flag_prop && is_prop_defined)
        {
          _db->deleteColumnsByLocator(ELoc::P);
        }
        if (!_useDiscrete)
        {
          _db->deleteColumnsByLocator(ELoc::L);
          _db->deleteColumnsByLocator(ELoc::U);
        }
        else
        {
          _db->deleteColumnsByLocator(ELoc::RKLOW);
          _db->deleteColumnsByLocator(ELoc::RKUP);
        }
        break;
    }
    return (0);
  }

  bool CalcModelPGS::_preparePropdefAndRule()
  {
    if (_flagStat)
    {
      if (_propdef.define(
            true, _flagStat, {_ngrf, 0}, {_nfacies, 0}, NULL, NULL, _props))
        return false;
      if (_rule != nullptr
          && _rule->particularities(NULL, NULL, NULL, 1, _flagStat))
        return false;
    }
    else
    {
      if (_propdef.define(
            true, _flagStat, {_ngrf, 0}, {_nfacies, 0}, _db, _dbprop, _props))
        return false;
      if (_rule != nullptr
          && _rule->particularities(_db, _dbprop, NULL, 1, _flagStat))
        return false;
    }

    _propdef.defineRuleMethod(EProcessOper::COPY);
    return true;
  }

  void CalcModelPGS::_defineBounds(
    Id iech1,
    Id iech2,
    Id ifac1,
    Id ifac2,
    double* low,
    double* up,
    double* ploc)
  {
    if (_flagStat)
    {
      ploc[0] = _propdef.getPropFix(ifac1);
      ploc[1] = _propdef.getPropFix(ifac2);
      low[0] = getStatThresh(ifac1, 0, 0);
      up[0] = getStatThresh(ifac1, 0, 1);
      low[1] = getStatThresh(ifac2, 0, 0);
      up[1] = getStatThresh(ifac2, 0, 1);
      if (_ngrf > 1)
      {
        low[2] = getStatThresh(ifac1, 1, 0);
        up[2] = getStatThresh(ifac1, 1, 1);
        low[3] = getStatThresh(ifac2, 1, 0);
        up[3] = getStatThresh(ifac2, 1, 1);
      }
    }
    else
    {
      ploc[0] = _db->getLocVariable(ELoc::P, iech1, ifac1);
      ploc[1] = _db->getLocVariable(ELoc::P, iech2, ifac2);

      if (!_useDiscrete)
      {
        low[0] = _db->getLocVariable(ELoc::L, iech1, ifac1);
        up[0] = _db->getLocVariable(ELoc::U, iech1, ifac1);
        low[1] = _db->getLocVariable(ELoc::L, iech2, ifac2);
        up[1] = _db->getLocVariable(ELoc::U, iech2, ifac2);
        if (_ngrf > 1)
        {
          low[2] = _db->getLocVariable(ELoc::L, iech1, _nfacies + ifac1);
          up[2] = _db->getLocVariable(ELoc::U, iech1, _nfacies + ifac1);
          low[3] = _db->getLocVariable(ELoc::L, iech2, _nfacies + ifac2);
          up[3] = _db->getLocVariable(ELoc::U, iech2, _nfacies + ifac2);
        }
      }
      else
      {
        low[0] = _db->getLocVariable(ELoc::RKLOW, iech1, ifac1);
        up[0] = _db->getLocVariable(ELoc::RKUP, iech1, ifac1);
        low[1] = _db->getLocVariable(ELoc::RKLOW, iech2, ifac2);
        up[1] = _db->getLocVariable(ELoc::RKUP, iech2, ifac2);
        if (_ngrf > 1)
        {
          low[2] = _db->getLocVariable(ELoc::RKLOW, iech1, _nfacies + ifac1);
          up[2] = _db->getLocVariable(ELoc::RKUP, iech1, _nfacies + ifac1);
          low[3] = _db->getLocVariable(ELoc::RKLOW, iech2, _nfacies + ifac2);
          up[3] = _db->getLocVariable(ELoc::RKUP, iech2, _nfacies + ifac2);
        }
      }
    }
  }

  double CalcModelPGS::_getProba(
    bool flagIndependent,
    double* low,
    double* up,
    const VectorInt& iconf,
    double* cov,
    Id maxpts)
  {
    double proba = TEST;

    if (_ngrf == 1)
    {
      proba = getProbaInd(cov[0], low, up, iconf[0]);
    }
    else
    {
      if (flagIndependent)
      {
        proba = getProbaInd(cov[0], low, up, iconf[0])
              * getProbaInd(cov[5], &low[2], &up[2], iconf[1]);
      }
      else
      {
        // Case when the two GRFs are correlated (presence of RHO)

        if (!_useDiscrete)
        {
          Id ier;
          double err;

          Id infin[4];
          infin[0] = mvndst_infin(low[0], up[0]);
          infin[1] = mvndst_infin(low[1], up[1]);
          infin[2] = mvndst_infin(low[2], up[2]);
          infin[3] = mvndst_infin(low[3], up[3]);
          mvndst(
            4, low, up, infin, cov, maxpts, _abseps, _releps, &err, &proba,
            &ier);
        }
        else
        {
          my_throw(
            "Discrete calculation for correlated GRFs is not performed yet");
        }
      }
    }

    return (proba);
  }

  double CalcModelPGS::_getValue(
    bool flagIndependent,
    Id iech1,
    Id iech2,
    Id ifac1,
    Id ifac2,
    const VectorInt& iconf,
    double* cov)
  {
    double value, ploc[2], low[4], up[4];
    if (_vario->getCalcul() == ECalcVario::VARIOGRAM)
    {
      if (ifac1 == ifac2)
      {
        _defineBounds(iech1, iech2, ifac1, ifac2, low, up, ploc);
        double g1 = _getProba(flagIndependent, low, up, iconf, cov);
        value = (ploc[0] + ploc[1]) * 0.5 - g1;
      }
      else
      {
        _defineBounds(iech1, iech2, ifac1, ifac2, low, up, ploc);
        double g1 = _getProba(flagIndependent, low, up, iconf, cov);
        _defineBounds(iech2, iech1, ifac1, ifac2, low, up, ploc);
        double g2 = _getProba(flagIndependent, low, up, iconf, cov);
        value = -0.5 * (g1 + g2);
      }
    }
    else
    {
      _defineBounds(iech1, iech2, ifac1, ifac2, low, up, ploc);
      value = _getProba(flagIndependent, low, up, iconf, cov);
    }
    return value;
  }

  bool CalcModelPGS::_calculCovmatrix(
    VectorDouble& d0,
    const VectorDouble& d1,
    VectorInt& iconf,
    double* cov)
  {

    Id nvar = _model->getNVar();
    MatrixSquare cov0(nvar);
    MatrixSquare covh(nvar);

    /* Calculate the covariance for the zero distance */
    _model->evaluateMatInPlace(nullptr, d0, cov0);

    /* Calculate the covariance for the given shift */
    _model->evaluateMatInPlace(nullptr, d1, covh);

    if (_rule->getModeRule() == ERule::STD)
    {
      cov[0] = covh.getValue(0, 0); /* C11(h)  */
      if (_ngrf > 1)
      {
        cov[1] = cov0.getValue(1, 0); /* C21(0)  */
        cov[2] = covh.getValue(0, 1); /* C21(-h) */
        cov[3] = covh.getValue(1, 0); /* C21(h)  */
        cov[4] = cov0.getValue(0, 1); /* C21(0)  */
        cov[5] = covh.getValue(1, 1); /* C22(h)  */
      }
    }
    else if (_rule->getModeRule() == ERule::SHIFT)
    {
      auto* ruleshift = static_cast<RuleShift*>(const_cast<Rule*>(_rule));
      cov[0] = covh.getValue(0, 0); /* C11(h)  */
      cov[5] =
        (nvar == 1) ? covh.getValue(0, 0) : covh.getValue(1, 1); /* C22(h)  */

      for (Id i = 0; i < static_cast<Id>(_model->getNDim()); i++)
        d0[i] = ruleshift->getShift(i);

      _model->evaluateMatInPlace(nullptr, d0, covh);
      cov[1] =
        (nvar == 1) ? covh.getValue(0, 0) : covh.getValue(1, 0); /* C21(s)  */
      cov[4] =
        (nvar == 1) ? covh.getValue(0, 0) : covh.getValue(1, 0); /* C21(s)  */

      for (Id i = 0; i < _ndim; i++) d0[i] = d1[i] - ruleshift->getShift(i);
      _model->evaluateMatInPlace(nullptr, d0, covh);
      cov[2] =
        (nvar == 1) ? covh.getValue(0, 0) : covh.getValue(1, 0); /* C21(h-s) */

      for (Id i = 0; i < _ndim; i++) d0[i] = d1[i] + ruleshift->getShift(i);
      _model->evaluateMatInPlace(nullptr, d0, covh);
      cov[3] =
        (nvar == 1) ? covh.getValue(0, 0) : covh.getValue(1, 0); /* C21(h+s)  */
    }
    else
      messageAbort("This rule is not expected in st_calcul_covmatrix");

    /* Check if the two GRFs are spatially independent */

    bool flagIndependent = true;
    if (_ngrf > 1)
    {
      for (Id i = 1; i <= 4; i++)
        if (ABS(cov[i]) > 1.e-8) flagIndependent = false;
    }

    /* In TEST_DISCRET case, identify the ranks of the discretized covariance */

    if (_useDiscrete)
    {
      double cround;
      iconf[0] = _discretepgs->getCovrank(cov[0], &cround);
      if (_ngrf > 1) iconf[1] = _discretepgs->getCovrank(cov[5], &cround);
    }
    return flagIndependent;
  }

  void CalcModelPGS::_updateVarianceStat()
  {
    for (Id ivar = 0; ivar < _nfacies; ivar++)
      for (Id jvar = 0; jvar < _nfacies; jvar++)
      {
        double pivar = _propdef.getPropFix(ivar);
        double pjvar = _propdef.getPropFix(jvar);
        if (ivar == jvar)
          _vario->setVar(pivar * (1. - pivar), ivar, jvar);
        else
          _vario->setVar(-pivar * pjvar, ivar, jvar);
        if (!_vario->getFlagAsym()) continue;

        for (Id idir = 0; idir < _vario->getNDir(); idir++)
        {
          Id iad = _vario->getAddressCenterForGg(idir, ivar, jvar);
          _vario->setSwByIndex(idir, iad, 1);
          _vario->setHhByIndex(idir, iad, 0);

          switch (_vario->getCalcul().toEnum())
          {
            case ECalcVario::E_VARIOGRAM: break;

            case ECalcVario::E_COVARIANCE:
              _vario->setGgByIndex(idir, iad, _vario->getVar(ivar, jvar));
              break;

            case ECalcVario::E_COVARIANCE_NC:
              _vario->setGgByIndex(idir, iad, (ivar == jvar) ? pivar : 0.);
              break;
            default: break;
          }
        }
      }
  }

  void CalcModelPGS::_updateVarianceNoStat()
  {
    VectorDouble mean(_nfacies, 0.);
    VectorDouble covs(_nfacies * _nfacies, 0.);

    /* Loop on the samples */

    Id number = 0;
    for (Id iech = 0; iech < _db->getNSample(); iech++)
    {
      if (!_db->isActive(iech)) continue;
      for (Id ivar = 0; ivar < _nfacies; ivar++)
      {
        double p1 = _db->getLocVariable(ELoc::P, iech, ivar);
        mean[ivar] += p1;

        for (Id jvar = 0; jvar < _nfacies; jvar++)
        {
          double p2 = _db->getLocVariable(ELoc::P, iech, jvar);
          COVS(ivar, jvar) += p1 * p2;
        }
      }
      number++;
    }

    /* Normalization */

    for (Id ivar = 0; ivar < _nfacies; ivar++)
      mean[ivar] /= static_cast<double>(number);
    for (Id ivar = 0; ivar < _nfacies; ivar++)
      for (Id jvar = 0; jvar < _nfacies; jvar++)
        COVS(ivar, jvar) = COVS(ivar, jvar) / static_cast<double>(number);

    /* Store the results */

    for (Id ivar = 0; ivar < _nfacies; ivar++)
      for (Id jvar = 0; jvar < _nfacies; jvar++)
      {
        _vario->setVar(COVS(ivar, jvar), ivar, jvar);

        if (!_vario->getFlagAsym()) continue;

        // Set the C00 term
        for (Id idir = 0; idir < _vario->getNDir(); idir++)
        {
          Id iad = _vario->getAddressCenterForGg(idir, ivar, jvar);
          _vario->setSwByIndex(idir, iad, _db->getNSample());
          _vario->setHhByIndex(idir, iad, 0);

          switch (_vario->getCalcul().toEnum())
          {
            case ECalcVario::E_VARIOGRAM: break;

            case ECalcVario::E_COVARIANCE:
              _vario->setGgByIndex(idir, iad, _vario->getVar(ivar, jvar));
              break;

            case ECalcVario::E_COVARIANCE_NC:
              _vario->setGgByIndex(idir, iad, (ivar == jvar) ? mean[ivar] : 0.);
              break;
            default: break;
          }
        }
      }
  }

  Id CalcModelPGS::_varioIndicNoStat()
  {
    double dist, cov[6];
    Id iech, jech, iad;
    bool flagIndependent;
    VectorInt iconf(2);
    VectorDouble d0(_ndim);
    VectorDouble d1(_ndim);

    /* Loop on the directions */

    if (_varioPgsVariable(1, 0, 1)) return 1;

    if (_varioPgsVariable(0, 0, 1)) return 1;

    for (Id idir = 0, ndir = _vario->getNDir(); idir < ndir; idir++)
    {
      /* Establish the geometry */

      if (_variogramGeometryPgsCalcul(idir)) return (1);
      if (_variogramGeometryPgsFinal()) return (1);

      /* Loop on the lags */

      for (Id ilag = 0, nlag = _vario->getNLag(idir); ilag < nlag; ilag++)
      {
        _vorder.getBounds(idir, ilag, &_ifirst, &_ilast);
        if (_ifirst >= _ilast) continue;

        /* Loop on the pairs of the lag */

        for (Id ipair = _ifirst; ipair < _ilast; ipair++)
        {
          _vorder.getIndices(ipair, &iech, &jech, &dist);

          /* Calculate the distance vector */

          dist = distance_intra(_db, iech, jech, d1.data());
          flagIndependent = _calculCovmatrix(d0, d1, iconf, cov);

          /* Loops on the facies */

          for (Id ifac = 0; ifac < _nfacies; ifac++)
            for (Id jfac = 0; jfac <= ifac; jfac++)
            {
              if (_vario->getFlagAsym())
              {
                iad = _vario->getAddressForGg(idir, ifac, jfac, ilag, 1);
                _vario->updateGgByIndex(
                  idir, iad,
                  _getValue(
                    flagIndependent, iech, jech, ifac, jfac, iconf, cov));
                iad = _vario->getAddressForGg(idir, ifac, jfac, ilag, -1);
                _vario->updateGgByIndex(
                  idir, iad,
                  _getValue(
                    flagIndependent, jech, iech, ifac, jfac, iconf, cov));
              }
              else
              {
                iad = _vario->getAddressForGg(idir, ifac, jfac, ilag, ITEST);
                _vario->updateGgByIndex(
                  idir, iad,
                  _getValue(
                    flagIndependent, iech, jech, ifac, jfac, iconf, cov));
              }
            }
        }
      }

      _updateVarianceNoStat();

      (void)_varioPgsVariable(-1, 0, 1);

      /* Clear the geometry */

      _vorder.clear();

      /* Scale the variogram */

      _vario->finalScaleByWeights(idir);
    }
    return (0);
  }

  Id CalcModelPGS::_varioIndicStat()
  {
    double cov[6];
    Id iad;
    bool flagIndependent;
    VectorInt iconf(2);
    VectorDouble d0(_ndim);
    VectorDouble d1(_ndim);

    /* Initializations */

    for (Id i = 0; i < 6; i++) cov[i] = 0.;
    if (_calculateThreshStat()) return 1;
    _copySwhh();

    /* Loop on the directions */

    for (Id idir = 0, ndir = _vario->getNDir(); idir < ndir; idir++)
    {

      /* Loop on the lags */

      for (Id ilag = 0, nlag = _vario->getNLag(idir); ilag < nlag; ilag++)
      {

        /* Calculate the distance vector */

        Id jpas = _vario->getAddressForGg(idir, 0, 0, ilag, 1);
        for (Id idim = 0; idim < _ndim; idim++)
          d1[idim] =
            _vario->getHhByIndex(idir, jpas) * _vario->getCodir(idir, idim);
        flagIndependent = _calculCovmatrix(d0, d1, iconf, cov);

        /* Loops on the facies */

        for (Id ifac = 0; ifac < _nfacies; ifac++)
          for (Id jfac = 0; jfac <= ifac; jfac++)
          {
            if (_vario->getFlagAsym())
            {
              iad = _vario->getAddressForGg(idir, ifac, jfac, ilag, 1);
              _vario->setGgByIndex(
                idir, iad,
                _getValue(flagIndependent, 0, 0, ifac, jfac, iconf, cov));
              iad = _vario->getAddressForGg(idir, ifac, jfac, ilag, -1);
              _vario->setGgByIndex(
                idir, iad,
                _getValue(flagIndependent, 0, 0, jfac, ifac, iconf, cov));
            }
            else
            {
              iad = _vario->getAddressForGg(idir, ifac, jfac, ilag, ITEST);
              _vario->setGgByIndex(
                idir, iad,
                _getValue(flagIndependent, 0, 0, ifac, jfac, iconf, cov));
            }
          }
      }
    }

    _updateVarianceStat();
    return (0);
  }

  /****************************************************************************/
  /*!
   **  Calculate the gaussian variograms in the stationary case
   **
   *****************************************************************************/
  Id CalcModelPGS::_varioPgsStat()
  {
    _flagFacies = 1;
    _flagDist = 0;
    _corpgs.define(0, 0, _rule->getRho());
    _tracepgs.define(_ngrf, _corpgs.getNpar(), false, false);

    setRho(0.);
    if (_calculateThreshStat()) return 1;
    _copySwhh();

    // Infer the variogram of PGS
    _varcalcFromVarioStat();

    return 0;
  }

  /****************************************************************************/
  /*!
   **  Discard a data if:
   **  - its facies is unknown
   **  - its thresholds are so close that it leads to a zero probability
   **
   ** \return  Error return code
   **
   ** \param[in]  iech        Rank of the sample
   **
   *****************************************************************************/
  Id CalcModelPGS::_discardPoint(Id iech)
  {
    double low, up;

    /* This function is bypassed in the stationary case */

    if (_flagStat) return (0);

    /* The following checks must not be performed if not on facies */

    if (!_flagFacies) return (0);

    /* Check on the facies */

    Id ifac = static_cast<Id>(_db->getZVariable(iech, 0));
    if (ifac < 1 || ifac > _nfacies) return (1);

    /* Check on the thresholds */

    if (!_useDiscrete)
    {
      if (_db->getNInterval() <= 0) return (0);
      low = _db->getLocVariable(ELoc::L, iech, _igrfcur);
      up = _db->getLocVariable(ELoc::U, iech, _igrfcur);
    }
    else
    {
      if (get_LOCATOR_NITEM(_db, ELoc::RKLOW) <= 0
          && get_LOCATOR_NITEM(_db, ELoc::RKUP) <= 0)
        return (0);
      low = _db->getLocVariable(ELoc::RKLOW, iech, _igrfcur);
      up = _db->getLocVariable(ELoc::RKUP, iech, _igrfcur);
    }
    if (up <= low) return (1);

    return (0);
  }

  /****************************************************************************/
  /*!
   **  Determine the Geometry of all pairs
   **
   ** \return  Error return code
   **
   ** \param[in]  idir        Rank of the direction
   **
   *****************************************************************************/
  Id CalcModelPGS::_variogramGeometryPgsCalcul(Id idir)
  {
    Id iad;
    SpaceTarget T1(_vario->getSpace());
    SpaceTarget T2(_vario->getSpace());

    /* Retrieve information */

    Id nech = _db->getNSample();
    Id nvar = _vario->getNVar();
    double maxdist = _vario->getMaximumDistance(idir);
    const DirParam& dirparam = _vario->getDirParam(idir);

    // Local variables to speed up calculations
    bool hasSel = _db->hasLocVariable(ELoc::SEL);
    bool hasWeight = _db->hasLocVariable(ELoc::W);
    double dist = 0.;

    /* Sort the data */

    VectorInt rindex = _db->getSortArray();

    /* Loop on the first point */

    for (Id iiech = 0; iiech < nech - 1; iiech++)
    {
      Id iech = rindex[iiech];
      if (hasSel && !_db->isActive(iech)) continue;
      if (hasWeight && FFFF(_db->getWeight(iech))) continue;
      if (_discardPoint(iech)) continue;
      _db->getSampleAsSTInPlace(iech, T1);
      mes_process("Calculating Variogram Geometry", nech, iech);

      for (Id jjech = iiech + 1; jjech < nech; jjech++)
      {
        Id jech = rindex[jjech];
        if (_db->getDistance1D(iech, jech) > maxdist) break;
        if (hasSel && !_db->isActive(jech)) continue;
        if (hasWeight && FFFF(_db->getWeight(jech))) continue;
        if (_discardPoint(jech)) continue;
        _db->getSampleAsSTInPlace(jech, T2);

        // Reject the point as soon as one BiTargetChecker is not correct
        if (!_vario->keepPair(idir, T1, T2, &dist)) continue;

        /* Get the rank of the lag */

        auto ilag = dirparam.getLagRank(dist);
        if (isNA(ilag)) continue;

        /* Add the sample (only positive lags are of interest) */

        if (ilag < 0) ilag = -ilag;
        if (_vorder.add(iech, jech, NULL, NULL, ilag, idir, dist)) return 1;
        dist = ABS(dist);

        /* Update the distance and weight for all GRFs */

        for (Id ivar = 0; ivar < nvar; ivar++)
          for (Id jvar = 0; jvar <= ivar; jvar++)
          {
            if (_vario->getFlagAsym())
            {
              iad = _vario->getAddressForGg(idir, ivar, jvar, ilag, 1);
              _vario->setGgByIndex(idir, iad, 0.);
              _vario->setHhByIndex(
                idir, iad, _vario->getHhByIndex(idir, iad) - dist);
              _vario->setSwByIndex(
                idir, iad, _vario->getSwByIndex(idir, iad) + 1);
              iad = _vario->getAddressForGg(idir, ivar, jvar, ilag, -1);
              _vario->setGgByIndex(idir, iad, 0.);
              _vario->setHhByIndex(
                idir, iad, _vario->getHhByIndex(idir, iad) + dist);
              _vario->setSwByIndex(
                idir, iad, _vario->getSwByIndex(idir, iad) + 1);
            }
            else
            {
              iad = _vario->getAddressForGg(idir, ivar, jvar, ilag);
              _vario->setGgByIndex(idir, iad, 0.);
              _vario->setHhByIndex(
                idir, iad, _vario->getHhByIndex(idir, iad) + dist);
              _vario->setSwByIndex(
                idir, iad, _vario->getSwByIndex(idir, iad) + 1);
            }
          }
      }
    }
    return 0;
  }

  /****************************************************************************/
  /*!
   **  Correct the experimental variogram for GRFs
   **
   ** \param[in]  idir        Rank of the direction
   **
   *****************************************************************************/
  void CalcModelPGS::_variogramGeometryPgsCorrect(Id idir)
  {
    Id iad;

    for (Id ilag = 0; ilag < _vario->getNLag(idir); ilag++)
      for (Id igrf = 0; igrf < _ngrf; igrf++)
        for (Id jgrf = 0; jgrf <= igrf; jgrf++)
        {
          iad = _vario->getAddressForGg(idir, igrf, jgrf, ilag, 1);
          _vario->setGgByIndex(idir, iad, _corpgs.paramExpand(igrf, jgrf, 1));
          if (_vario->getSwByIndex(idir, iad) > 0.)
            _vario->setHhByIndex(
              idir, iad,
              _vario->getHhByIndex(idir, iad)
                / _vario->getSwByIndex(idir, iad));
          iad = _vario->getAddressForGg(idir, igrf, jgrf, ilag, -1);
          _vario->setGgByIndex(idir, iad, _corpgs.paramExpand(igrf, jgrf, -1));
          if (_vario->getSwByIndex(idir, iad) > 0.)
            _vario->setHhByIndex(
              idir, iad,
              _vario->getHhByIndex(idir, iad)
                / _vario->getSwByIndex(idir, iad));
        }
  }

  /****************************************************************************/
  /*!
   **  Compress the Geometry of all pairs
   **
   ** \return  Error return code
   **
   *****************************************************************************/
  Id CalcModelPGS::_variogramGeometryPgsFinal()
  {
    Id npair = _vorder.final();
    if (_vorder.empty()) return (1);
    if (npair > 0 && !_flagStat) _memint.resize(npair);
    return (0);
  }

  bool CalcModelPGS::_isLagUnused(Id idir, Id ilag)
  {
    return (
      _vario->getSwByIndex(idir, _vario->getNLag(idir) + ilag + 1) <= 0
      || _vario->getUtilizeByIndex(idir, _vario->getNLag(idir) + ilag + 1)
           == 0);
  }

  /****************************************************************************/
  /*!
   ** PRUPOSE: Internal function
   **
   *****************************************************************************/
  double CalcModelPGS::_st_rkl(
    double x,
    double y,
    VectorDouble& lower,
    VectorDouble& upper,
    MatrixSymmetric& correl,
    MatrixSquare& covar,
    MatrixSquare& temp,
    Id maxpts)
  {
    double v2, error;
    Id inform;

    VectorDouble cste{0., 0.};
    VectorDouble vec{x, y};
    VectorDouble mean = AMatrix::product(temp, vec);
    double v1 = law_df_bigaussian(vec, cste, correl);
    mvndst2n(
      lower.data(), upper.data(), mean.data(), covar.getValues().data(), maxpts,
      _abseps, _releps, &error, &v2, &inform);
    return (v1 * v2);
  }

  /****************************************************************************/
  /*!
   ** \return  Calculate
   ** \return    d/(d xk) * d/(d xl)
   ** \return    int_lower1^upper1
   ** \return    int_lower2^upper2
   ** \return    int_lower3^upper3
   ** \return    int_lower4^upper4
   ** \return    gaussian density(x1,x2,x3,x4,rho) dx1 dx2 dx3 dx4
   **
   ** \param[in]  index1       First derivative index
   ** \param[in]  index2       Second derivative index
   ** \param[in]  lower        Array of lower bounds
   ** \param[in]  upper        Array of upper bounds
   ** \param[in]  correl       Correlation matrix (Dimension = 4*4)
   ** \param[in]  maxpts       Maximum number of evaluations
   **
   *****************************************************************************/
  double CalcModelPGS::_st_ikl(
    Id index1,
    Id index2,
    VectorDouble& lower,
    VectorDouble& upper,
    MatrixSymmetric& correl,
    Id maxpts)
  {
    VectorInt index = {index1, index2};

    // Build submatrices
    VectorDouble low = VH::reduce(lower, index);
    VectorDouble upp = VH::reduce(upper, index);
    auto* corr1 = dynamic_cast<MatrixSymmetric*>(
      MatrixFactory::createReduce(&correl, index, index, true, true));
    auto* corrc = dynamic_cast<MatrixSymmetric*>(
      MatrixFactory::createReduce(&correl, index, index, false, true));
    auto* corr2 = dynamic_cast<MatrixSymmetric*>(
      MatrixFactory::createReduce(&correl, index, index, false, false));
    MatrixSymmetric inv_corr1(*corr1);
    if (inv_corr1.invert()) messageAbort("st_ikl #1");
    auto* temp =
      dynamic_cast<MatrixSquare*>(MatrixFactory::prodMatMat(corrc, &inv_corr1));

    // Derive covar
    MatrixSquare covar(2);
    for (Id i = 0; i < 2; i++)
      for (Id j = 0; j < 2; j++)
      {
        double value = 0;
        for (Id k = 0; k < 2; k++)
          value += temp->getValue(k, i) * corrc->getValue(k, j);
        covar.setValue(i, j, corr2->getValue(i, j) - value);
      }

    double x, y;
    double S = 0.;
    x = upper[index1];
    if (IS_GAUSS_DEF(x))
    {
      y = upper[index2];
      if (IS_GAUSS_DEF(y))
        S += _st_rkl(x, y, low, upp, *corr1, covar, *temp, maxpts);
      y = lower[index2];
      if (IS_GAUSS_DEF(y))
        S -= _st_rkl(x, y, low, upp, *corr1, covar, *temp, maxpts);
    }
    x = lower[index1];
    if (IS_GAUSS_DEF(x))
    {
      y = lower[index2];
      if (IS_GAUSS_DEF(y))
        S += _st_rkl(x, y, low, upp, *corr1, covar, *temp, maxpts);
      y = upper[index2];
      if (IS_GAUSS_DEF(y))
        S -= _st_rkl(x, y, low, upp, *corr1, covar, *temp, maxpts);
    }
    return (S / 2.);
  }

  /****************************************************************************/
  /*!
   **  Calculate the gradients for a pair of facies and return it
   **
   ** \return  Calculate d²S/dC12dC21 and d²S/dC1dC2
   **
   *****************************************************************************/
  double CalcModelPGS::_st_d2_dkldij(
    VectorDouble& lower,
    VectorDouble& upper,
    MatrixSymmetric& correl)
  {
    Id grid[4];
    VectorDouble u(4);

    double S = 0.;
    for (Id i4 = 0; i4 < 2; i4++)
      for (Id i3 = 0; i3 < 2; i3++)
        for (Id i2 = 0; i2 < 2; i2++)
          for (Id i1 = 0; i1 < 2; i1++)
          {
            grid[0] = i1;
            grid[1] = i2;
            grid[2] = i3;
            grid[3] = i4;

            bool flag_out = false;
            for (Id i = 0; i < 4 && flag_out == 0; i++)
            {
              u[i] = (grid[i]) ? upper[i] : lower[i];
              flag_out = ISNOT_GAUSS_DEF(u[i]);
            }
            if (!flag_out)
              S += pow(-1., i1 + i2 + i3 + i4) * law_df_quadgaussian(u, correl);
          }
    return (S / 2.);
  }

  double CalcModelPGS::_st_nkl(
    VectorDouble& u,
    double lower,
    double upper,
    VectorDouble& invvari,
    Id index2,
    double meanj,
    double varj,
    double stdj)
  {
    double dfupp = law_dnorm(upper, meanj, stdj);
    double dflow = law_dnorm(lower, meanj, stdj);
    double cdfupp = law_cdf_gaussian((upper - meanj) / stdj);
    double cdflow = law_cdf_gaussian((lower - meanj) / stdj);
    double invval = invvari[index2];
    VectorDouble invpart = VH::reduceOne(invvari, index2);
    double total = VH::innerProductCV(invpart, u);
    double S = (dfupp - dflow) * varj * invval
             - (cdfupp - cdflow) * (invval * meanj + total);
    return (S);
  }

  /****************************************************************************/
  /*!
   ** \return  Calculate the first-order derivatives
   **
   ** \param[in]  index1      First derivative index
   ** \param[in]  index2      Second derivative index
   ** \param[in]  lower       Array of lower bounds (Dimension = 4)
   ** \param[in]  upper       Array of upper bounds (Dimension = 4)
   ** \param[in]  correl      Correlation matrix (Dimension = 4*4)
   **
   *****************************************************************************/
  double CalcModelPGS::_st_d2_dkldkj(
    Id index1,
    Id index2,
    VectorDouble& lower,
    VectorDouble& upper,
    MatrixSymmetric& correl)
  {
    auto* varcori = dynamic_cast<MatrixSymmetric*>(
      MatrixFactory::createReduceOne(&correl, index2, index2, false, false));
    MatrixSymmetric invvarcor(correl);
    if (invvarcor.invert()) messageAbort("st_d2_dkldkj #1");
    VectorDouble invvarcori = invvarcor.getRow(index1);
    auto* corr1 = dynamic_cast<MatrixSymmetric*>(
      MatrixFactory::createReduceOne(&correl, index2, index2, false, false));
    VectorDouble crosscor = correl.getRow(index2);
    crosscor = VH::reduceOne(crosscor, index2);
    double corr2 = correl.getValue(index2, index2);

    MatrixSymmetric invcorr1(*corr1);
    if (invcorr1.invert()) messageAbort("st_d2_dkldkj #2");

    VectorDouble temp = AMatrix::product(invcorr1, crosscor);
    double covar = invcorr1.normVec(crosscor);
    double sdcovar = sqrt(corr2 - covar);

    VectorDouble lowi = VH::reduceOne(lower, index2);
    VectorDouble uppi = VH::reduceOne(upper, index2);
    Id lowj = lower[index2];
    Id uppj = upper[index2];

    double S = 0.;
    VectorDouble u(3);
    Id grid[3];

    for (Id i3 = 0; i3 < 2; i3++)
      for (Id i2 = 0; i2 < 2; i2++)
        for (Id i1 = 0; i1 < 2; i1++)
        {
          grid[0] = i1;
          grid[1] = i2;
          grid[2] = i3;

          bool flag_out = false;
          for (Id i = 0; i < 3 && !flag_out; i++)
          {
            u[i] = (grid[i]) ? uppi[i] : lowi[i];
            flag_out = ISNOT_GAUSS_DEF(u[i]);
          }
          if (flag_out) continue;
          double mu = VH::innerProductCV(temp, u);
          double random = law_df_multigaussian(u, *varcori);

          S += pow(-1., 3 - i1 + i2 + i3) * random
             * _st_nkl(u, lowj, uppj, invvarcori, index2, mu, covar, sdcovar);
        }
    return (S / 2);
  }

  /****************************************************************************/
  /*!
   ** \return  Calculate second-order derivatives
   **
   ** \param[in]  index1       First derivative index
   ** \param[in]  index2       Second derivative index
   ** \param[in]  lower        Array of lower bounds (Dimension = 4)
   ** \param[in]  upper        Array of upper bounds (Dimension = 4)
   ** \param[in]  correl       Correlation matrix (Dimension = 4*4)
   ** \param[in]  maxpts       Maximum number of evaluations
   **
   *****************************************************************************/
  double CalcModelPGS::_st_d2_dkldkl(
    Id index1,
    Id index2,
    VectorDouble& lower,
    VectorDouble& upper,
    MatrixSymmetric& correl,
    Id maxpts)
  {
    MatrixSymmetric corri = correl;
    corri.updValue(index1, index2, EOperator::ADD, _deltaparam);
    double v1 = _st_ikl(index1, index2, lower, upper, corri, maxpts);

    corri = correl;
    corri.updValue(index1, index2, EOperator::SUBTRACT, _deltaparam);
    double v2 = _st_ikl(index1, index2, lower, upper, corri, maxpts);

    double S = (v1 - v2) / (2. * _deltaparam);
    return (S / 2.);
  }

  /****************************************************************************/
  /*!
   **  Count the number of pairs with the target facies
   **
   ** \param[in]  ifac1        First target facies (starting from 1)
   ** \param[in]  ifac2        Second target facies (starting from 1)
   **
   *****************************************************************************/
  Id CalcModelPGS::_getCount(Id ifac1, Id ifac2)
  {
    Id i1, i2;
    double dist;

    /* Initializations */

    Id number = 0;
    for (Id ipair = _ifirst; ipair < _ilast; ipair++)
    {
      _vorder.getIndices(ipair, &i1, &i2, &dist);
      if (ifac1 != _db->getZVariable(i1, 0)) continue;
      if (ifac2 != _db->getZVariable(i2, 0)) continue;
      double w1 = _db->getWeight(i1);
      double w2 = _db->getWeight(i2);
      number += static_cast<Id>(w1 * w2);
    }
    return (number);
  }

  double CalcModelPGS::_calculGlobalScoreStat(
    Id flag_deriv,
    Id flag_reset,
    MatrixSymmetric& correl,
    VectorDouble& Grad,
    MatrixSymmetric& Hess,
    MatrixSymmetric& JJ,
    Id maxpts2,
    Id maxpts4)
  {
    double s, rj2, erval, ggval;
    Id inform;
    MatrixSymmetric hess(4);
    MatrixSymmetric gradgrad(4);
    VectorDouble grad(4);
    VectorDouble lower(4);
    VectorDouble upper(4);

    double S = 0.;
    grad.fill(0.);
    hess.fill(0.);

    for (Id ifac1 = 0; ifac1 < _nfacies; ifac1++)
    {
      for (Id ifac2 = 0; ifac2 < _nfacies; ifac2++)
      {
        auto nfifj = _getCount(ifac1 + 1, ifac2 + 1);
        if (nfifj <= 0) continue;

        /* Get the bounds */

        auto bounds = _rule->getThresh(ifac1 + 1);
        lower[0] = bounds[0];
        upper[0] = bounds[1];
        lower[2] = bounds[2];
        upper[2] = bounds[3];
        bounds = _rule->getThresh(ifac2 + 1);
        lower[1] = bounds[0];
        upper[1] = bounds[1];
        lower[2] = bounds[2];
        lower[3] = bounds[3];
        if (flag_reset)
        {
          mvndst4(
            lower.data(), upper.data(), correl.getValues().data(), _maxpts4,
            _abseps, _releps, &erval, &s, &inform);
          setStatProba(ifac1, ifac2, s);
        }
        else
          s = getStatProba(ifac1, ifac2);

        rj2 = -2. * log(s);
        S += static_cast<double>(nfifj) * rj2;

        /* Calculate the derivative */

        if (!flag_deriv) continue;
        grad[0] = -_st_ikl(0, 1, lower, upper, correl, maxpts2) / s;
        grad[1] = -_st_ikl(0, 3, lower, upper, correl, maxpts2) / s;
        grad[2] = -_st_ikl(1, 2, lower, upper, correl, maxpts2) / s;
        grad[3] = -_st_ikl(2, 3, lower, upper, correl, maxpts2) / s;
        gradgrad.fill(0.);
        for (Id i = 0; i < 4; i++)
          for (Id j = 0; j <= i; j++)
            gradgrad.setValue(i, j, grad[i] * grad[j]);

        hess.setValue(3, 0, -_st_d2_dkldij(lower, upper, correl));
        hess.setValue(2, 1, -_st_d2_dkldij(lower, upper, correl));
        hess.setValue(1, 0, _st_d2_dkldkj(0, 2, lower, upper, correl));
        hess.setValue(2, 0, _st_d2_dkldkj(1, 3, lower, upper, correl));
        hess.setValue(3, 1, _st_d2_dkldkj(3, 1, lower, upper, correl));
        hess.setValue(3, 2, _st_d2_dkldkj(2, 0, lower, upper, correl));
        hess.setValue(0, 0, _st_d2_dkldkl(0, 1, lower, upper, correl, maxpts4));
        hess.setValue(1, 1, _st_d2_dkldkl(0, 3, lower, upper, correl, maxpts4));
        hess.setValue(2, 2, _st_d2_dkldkl(1, 2, lower, upper, correl, maxpts4));
        hess.setValue(3, 3, _st_d2_dkldkl(2, 3, lower, upper, correl, maxpts4));

        for (Id i = 0; i < 4; i++)
        {
          Grad[i] += nfifj * grad[i];
          for (Id j = 0; j <= i; j++)
          {
            ggval = gradgrad.getValue(i, j);
            Hess.updValue(
              i, j, EOperator::SUBTRACT,
              nfifj * (hess.getValue(i, j) / s - ggval));
            JJ.updValue(i, j, EOperator::ADD, nfifj * ggval / rj2);
          }
        }
      }
    }
    return (S / 2.);
  }

  double CalcModelPGS::_calculGlobalScoreNoStat(
    Id flag_deriv,
    Id flag_reset,
    MatrixSymmetric& correl,
    VectorDouble& Grad,
    MatrixSymmetric& Hess,
    MatrixSymmetric& JJ,
    Id maxpts2,
    Id maxpts4)
  {
    double s, erval, dist;
    Id i1, i2, inform;
    VectorDouble grad(4);
    VectorDouble lower(4);
    VectorDouble upper(4);
    MatrixSymmetric hess(4);
    MatrixSymmetric gradgrad(4);

    double S = 0.;
    grad.fill(0.);
    hess.fill(0.);

    for (Id ipair = _ifirst; ipair < _ilast; ipair++)
    {
      _vorder.getIndices(ipair, &i1, &i2, &dist);
      Id ifac1 = static_cast<Id>(_db->getZVariable(i1, 0));
      Id ifac2 = static_cast<Id>(_db->getZVariable(i2, 0));
      double w1 = _db->getWeight(i1);
      double w2 = _db->getWeight(i2);

      /* Get the bounds */

      (void)_propdef.ruleThreshDefine(
        _db, _rule, ifac1, i1, 0, 0, 1, lower.data(), upper.data(), &lower[2],
        &upper[2]);
      (void)_propdef.ruleThreshDefine(
        _db, _rule, ifac2, i2, 0, 0, 1, &lower[1], &upper[1], &lower[3],
        &upper[3]);

      if (flag_reset)
      {
        mvndst4(
          lower.data(), upper.data(), correl.getValues().data(), _maxpts4,
          _abseps, _releps, &erval, &s, &inform);
        setMemInt(ipair, s);
      }
      else
        s = getMemInt(ipair);
      double rj2 = -2. * log(s);
      S += w1 * w2 * rj2;

      /* Calculate the derivative */

      if (!flag_deriv) continue;
      grad[0] = -_st_ikl(0, 1, lower, upper, correl, maxpts2) / s;
      grad[1] = -_st_ikl(0, 3, lower, upper, correl, maxpts2) / s;
      grad[2] = -_st_ikl(1, 2, lower, upper, correl, maxpts2) / s;
      grad[3] = -_st_ikl(2, 3, lower, upper, correl, maxpts2) / s;
      gradgrad.fill(0.);
      for (Id i = 0; i < 4; i++)
        for (Id j = 0; j <= i; j++) gradgrad.setValue(i, j, grad[i] * grad[j]);

      hess.setValue(3, 0, -_st_d2_dkldij(lower, upper, correl));
      hess.setValue(2, 1, -_st_d2_dkldij(lower, upper, correl));
      hess.setValue(1, 0, _st_d2_dkldkj(0, 2, lower, upper, correl));
      hess.setValue(2, 0, _st_d2_dkldkj(1, 3, lower, upper, correl));
      hess.setValue(3, 1, _st_d2_dkldkj(3, 1, lower, upper, correl));
      hess.setValue(3, 2, _st_d2_dkldkj(2, 0, lower, upper, correl));
      hess.setValue(0, 0, _st_d2_dkldkl(0, 1, lower, upper, correl, maxpts4));
      hess.setValue(1, 1, _st_d2_dkldkl(0, 3, lower, upper, correl, maxpts4));
      hess.setValue(2, 2, _st_d2_dkldkl(1, 2, lower, upper, correl, maxpts4));
      hess.setValue(3, 3, _st_d2_dkldkl(2, 3, lower, upper, correl, maxpts4));

      for (Id i = 0; i < 4; i++)
      {
        Grad[i] += w1 * w2 * grad[i];
        for (Id j = 0; j <= i; j++)
        {
          double ggval = gradgrad.getValue(i, j);
          Hess.updValue(
            i, j, EOperator::SUBTRACT,
            w1 * w2 * (hess.getValue(i, j) / s - ggval));
          JJ.updValue(i, j, EOperator::ADD, w1 * w2 * ggval / rj2);
        }
      }
    }
    return (S / 2.);
  }

  /****************************************************************************/
  /*!
   **  Global calculation
   **
   ** \return  The global score
   **
   ** \param[in]  flag_deriv   1 if the derivatives must be calculated
   ** \param[in]  flag_reset   1 to update the probability calculations
   ** \param[in]  params       Array of parameters
   **
   ** \param[out] Grad         Vector of cumulated gradients (Dimension= 4)
   ** \param[out] Hess         Matrix of cumulated Hessian (Dimension= 4*4)
   ** \param[out] JJ           Matrix of cumulated JJ (Dimension= 4*4)
   **
   *****************************************************************************/
  double CalcModelPGS::_calculGlobalScore(
    bool flag_deriv,
    bool flag_reset,
    const VectorDouble& params,
    VectorDouble& Grad,
    MatrixSymmetric& Hess,
    MatrixSymmetric& JJ)
  {
    MatrixSymmetric correl(4);
    double S = 0.;

    if (flag_deriv)
    {
      Grad.fill(0.);
      Hess.fill(0.);
      JJ.fill(0.);
    }

    _corpgs.buildCorrel(params, correl);

    if (_flagStat)
      S =
        _calculGlobalScoreStat(flag_deriv, flag_reset, correl, Grad, Hess, JJ);
    else
      S = _calculGlobalScoreNoStat(
        flag_deriv, flag_reset, correl, Grad, Hess, JJ);

    /* Modify the results due to the constraints on parameters */

    if (flag_deriv) _corpgs.updateConstraintsWithJJ(Grad, Hess, JJ);

    return (S / 2.);
  }

  double CalcModelPGS::_optimOneLagPgs(
    Id new_val,
    double maxiter,
    double delta0,
    double tolsort)
  {
    MatrixSymmetric invGn(4);
    MatrixSymmetric correl(4);
    MatrixSymmetric Hess(4);
    MatrixSymmetric Gn(4);
    MatrixSymmetric JJ(4);
    MatrixSymmetric d2(4);
    VectorDouble Grad(4);
    VectorDouble d1(4);
    VectorDouble gr(4);
    VectorDouble hgn(4);
    VectorDouble step(4);
    VectorDouble hsd(4);
    VectorDouble a(4);
    VectorDouble hgna(4);
    VectorDouble eigval(4);
    const MatrixSquare* eigvec = nullptr;

    double Spen = 0.;
    double Srpen = 0.;
    double penalize = 1000;
    bool barrier = false;

    /* Initializations */

    Id npar = _corpgs.getNpar();
    double delta = delta0;
    double mdiminution = 0.;

    if (new_val) _corpgs.initializeParams();

    VectorDouble param_temp = _corpgs.getParams();
    Grad.fill(0.);
    Hess.fill(0.);
    JJ.fill(0.);

    /* Calculate the score and the derivatives */

    double Sr =
      _calculGlobalScore(true, true, _corpgs.getParams(), Grad, Hess, JJ);
    double niter = 0.;
    double Snew = 0.;
    bool flag_sortie = false;
    bool flag_moved = true;

    while (!flag_sortie)
    {
      if (barrier)
      {
        _corpgs.buildCorrel(param_temp, correl);
        auto eigenvectors = EigenVectors(correl);
        eigval = eigenvectors.getEigenValues();
        eigvec = &eigenvectors.getEigenVectors();
        _corpgs.derivativeEigen(eigval[3], eigvec, d1, d2);
        Srpen = Sr - penalize * log(eigval[3]);
        VH::linearCombinationInPlace(1., Grad, -penalize / eigval[3], d1, Grad);
        AMatrix::linearCombinationInPlace(
          Hess, 0., npar, Hess, penalize / (eigval[3] * eigval[3]), d2);
        AMatrix::linearCombinationInPlace(
          JJ, 0., npar, JJ, penalize / (eigval[3] * eigval[3]), d2);
        penalize /= 2.;
      }
      niter++;
      double delta2 = delta * delta;
      if (flag_moved)
      {
        gr = Grad;
        Gn = Hess;
        if (!Gn.isDefinitePositive()) Gn = JJ;
        VH::linearCombinationInPlace(-1., gr, 0., VectorDouble(), hsd);
        invGn = Gn;
        if (invGn.invert()) messageAbort("st_optim_lag");
        hgn = AMatrix::product(invGn, hsd);
      }

      /* Determine the lag (hgn, alpha*hsd) or a convex combination of both */

      if (hgn.innerProduct(hgn, npar) <= delta2)
      {
        step = hgn;
      }
      else
      {
        double normgrad2 =
          gr.norm2(); // TODO: check that gr() shoumd be of dimension npar
        double alpha = normgrad2 / Gn.normVec(gr);
        double normgrad = sqrt(normgrad2);
        if (normgrad > (delta / alpha))
        {
          VH::linearCombinationInPlace(
            delta / normgrad, hsd, 0., VectorDouble(), step);
        }
        else
        {
          VH::linearCombinationInPlace(alpha, hsd, 0., VectorDouble(), a);
          VH::linearCombinationInPlace(1., hgn, -1., a, hgna);
          double c = a.innerProduct(hgn);
          double a2 = a.norm2();
          double hgna2 = hgna.norm2();
          double beta = 0.;
          if (c <= 0.)
            beta = (-c + sqrt(c * c + hgna2 * (delta2 - a2))) / hgna2;
          else
            beta = (delta2 - a2) / (c + sqrt(c * c + hgna2 * (delta2 - a2)));
          VH::linearCombinationInPlace(beta, hgn, 1. - beta, a, step);
        }
      }

      VH::linearCombinationInPlace(
        1., step, 1., _corpgs.getParams(), param_temp);
      _corpgs.buildCorrel(param_temp, correl);
      while (!correl.isDefinitePositive())
      {
        VH::linearCombinationInPlace(0.9, step, 0., step, step);
        VH::linearCombinationInPlace(
          1.0, step, 1., _corpgs.getParams(), param_temp);
        _corpgs.buildCorrel(param_temp, correl);
      }

      Snew = _calculGlobalScore(false, true, param_temp, Grad, Hess, JJ);

      if (barrier) Spen = Snew - penalize * log(eigval[3]);
      double rval = 0.;
      if (!FFFF(Snew))
      {

        mdiminution = Snew - Sr;
        if (barrier) mdiminution = Spen - Srpen;
        double stepgr = VH::innerProductCV(step, gr);
        double mdiminution_pred = stepgr + 0.5 * Gn.normVec(step);
        rval = mdiminution / mdiminution_pred;
        flag_moved = (mdiminution < 0);
      }
      else
      {
        flag_moved = 0;
        rval = 0.;
      }

      if (flag_moved)
      {
        Sr = Snew;
        Srpen = Spen;
        VectorDouble temp(4);
        VH::linearCombinationInPlace(1, param_temp, 0, VectorDouble(), temp);
        _corpgs.setParams(temp);
        Snew =
          _calculGlobalScore(true, false, _corpgs.getParams(), Grad, Hess, JJ);
        if (barrier)
        {
          _corpgs.derivativeEigen(eigval[3], eigvec, d1, d2);
          VH::linearCombinationInPlace(1, Grad, penalize / eigval[3], d1, Grad);
          AMatrix::linearCombinationInPlace(
            Hess, 0., npar, Hess, -penalize / (eigval[3] * eigval[3]), d2);
          AMatrix::linearCombinationInPlace(
            JJ, 0., npar, JJ, -penalize / (eigval[3] * eigval[3]), d2);
          penalize /= 2.;
        }
        if (rval > 0.75) delta = MAX(delta, 3. * sqrt(step.norm2()));
      }
      if (rval < 0.25) delta /= 2.;

      flag_sortie =
        (step.norm(0) < tolsort || niter == maxiter || Grad.norm(0) < 0.05
         || (fabs(mdiminution) < tolsort && flag_moved));
    }

    /* Returning arguments */

    if (OptDbg::query(EDbg::CONVERGE))
    {
      message("Lag %d - S = %lf Parameters =", _ipascur, Sr);
      for (Id i = 0; i < _corpgs.getNpar(); i++)
        message(" %lf", _corpgs.getParam(i));
      message("\n");
    }

    /* Store the trace */

    _tracepgs.addRow();
    _tracepgs.update(niter, Snew, 0, npar, Grad.data());

    return (Sr);
  }

  /****************************************************************************/
  /*!
   **  Evaluate the variogram of the underlying GRFs (assuming the two GRFs
   **  of the PGS model are correlated)
   **
   ** \param[in]  idir      Rank of the direction
   **
   *****************************************************************************/
  double CalcModelPGS::varcalcCorrelatedGRF(Id idir)
  {
    Id iad;

    Id opt_temp = _corpgs.getOptCorrel();
    double value = 0.;
    for (Id ilag = 0, nlag = _vario->getNLag(idir); ilag < nlag; ilag++)
    {
      mes_process("Inverting Variogram Lag", _vario->getNLag(idir), ilag);
      _ipascur = ilag;
      _tracepgs.addRow();
      if (_isLagUnused(idir, ilag)) continue;
      _vorder.getBounds(idir, ilag, &_ifirst, &_ilast);
      if (_ifirst >= _ilast) continue;

      if (opt_temp != 2) _corpgs.setOptCorrel(2);

      _optimOneLagPgs(1);
      _corpgs.setOptCorrel(opt_temp);
      value +=
        (_vario->getUtilizeByIndex(idir, _vario->getNLag(idir) + ilag)
         * _optimOneLagPgs(0));

      for (Id igrf = 0; igrf < _ngrf; igrf++)
        for (Id jgrf = 0; jgrf <= igrf; jgrf++)
        {
          iad = _vario->getAddressForGg(idir, igrf, jgrf, ilag, 1);
          _vario->setGgByIndex(idir, iad, _corpgs.paramExpand(igrf, jgrf, 1));
          iad = _vario->getAddressForGg(idir, igrf, jgrf, ilag, -1);
          _vario->setGgByIndex(idir, iad, _corpgs.paramExpand(igrf, jgrf, -1));
        }
    }
    return (value);
  }

  void CalcModelPGS::varcalcUncorrelatedGRF(Id idir)
  {
    Id iad;
    double testval, niter;

    /* Loop on the lags */

    for (Id ilag = 0, nlag = _vario->getNLag(idir); ilag < nlag; ilag++)
    {
      mes_process("Inverting Variogram Lag", _vario->getNLag(idir), ilag);
      _ipascur = ilag;
      _tracepgs.addRow();
      if (_isLagUnused(idir, ilag)) continue;
      _vorder.getBounds(idir, ilag, &_ifirst, &_ilast);
      if (_ifirst >= _ilast) continue;

      for (Id igrf = 0; igrf < _ngrf; igrf++)
      {
        _igrfcur = igrf;
        double result = golden_search(
          _goldenSearchFuncNoStat, static_cast<void*>(this), _epsGoldenSearch,
          -1., 1., &testval, &niter);
        _tracepgs.update(idir + 1, ilag + 1, 2 * igrf, 1, &testval);
        _tracepgs.update(idir + 1, ilag + 1, 2 * igrf + 1, 1, &niter);

        for (Id jgrf = 0; jgrf <= igrf; jgrf++)
        {
          double varloc = (igrf == jgrf) ? result : 0.;
          iad = _vario->getAddressForGg(idir, igrf, jgrf, ilag, 1);
          _vario->setGgByIndex(idir, iad, varloc);
          iad = _vario->getAddressForGg(idir, igrf, jgrf, ilag, -1);
          _vario->setGgByIndex(idir, iad, varloc);

          if (OptDbg::query(EDbg::CONVERGE))
            message(
              "Lag:%d - Grf:%d - Variogram(%d) = %lf\n", ilag, igrf, iad,
              varloc);
        }
      }
    }
  }

  /****************************************************************************/
  /*!
   **  Performing the variogram calculations
   **
   ** \return  Error return code
   **
   ** \param[out] flag_geometry 1 if Geometry must be established per direction
   **                           0 if Geometry is already calculated before
   **                             calling this function
   **
   *****************************************************************************/
  Id CalcModelPGS::_variopgsCalculNoRho(Id flag_geometry)
  {
    setRho(_rule->getRho());

    /* Loop on the directions */

    for (Id idir = 0, ndir = _vario->getNDir(); idir < ndir; idir++)
    {
      _idircur = idir;

      /* Establish the geometry */

      if (flag_geometry)
      {
        if (_variogramGeometryPgsCalcul(idir)) return (1);
        _variogramGeometryPgsCorrect(idir);
        if (_variogramGeometryPgsFinal()) return (1);
      }

      /* Set the value of C(0) */

      _variogramPatchC00(idir, _rule->getRho());
      if (_ngrf > 1 && (_optCorrel != 2 || _rule->getRho() != 0))
        varcalcCorrelatedGRF(idir);
      else
        varcalcUncorrelatedGRF(idir);

      /* Clear the geometry */

      if (flag_geometry) _vorder.clear();
    }
    return (0);
  }

  void CalcModelPGS::_makeSomeLagsInactive()
  {
    for (Id idir = 0, ndir = _vario->getNDir(); idir < ndir; idir++)
      for (Id ilag = 0, nlag = _vario->getNLag(idir); ilag < nlag; ilag++)
        _vario->setUtilizeByIndex(idir, _vario->getNLag(idir) + ilag, 1.);
  }

  void CalcModelPGS::_makeAllLagsActive()
  {
    for (Id idir = 0, ndir = _vario->getNDir(); idir < ndir; idir++)
      for (Id ilag = 0, nlag = _vario->getNLag(idir); ilag < nlag; ilag++)
        _vario->setUtilizeByIndex(idir, _vario->getNLag(idir) + ilag, 1.);
  }

  Id CalcModelPGS::_variopgsCalculRho()
  {
    double testval, niter;

    /* Calculate the geometry */

    for (Id idir = 0, ndir = _vario->getNDir(); idir < ndir; idir++)
    {
      if (_variogramGeometryPgsCalcul(idir)) return (1);
      _variogramGeometryPgsCorrect(idir);
    }
    if (_variogramGeometryPgsFinal()) return (1);

    _makeSomeLagsInactive();
    setRho(golden_search(
      _geoldenSearchFuncRho, static_cast<void*>(this), _epsGoldenSearch, -1.,
      1., &testval, &niter));
    _makeAllLagsActive();

    /* Perform the calculations with fixed rho */

    if (_variopgsCalculNoRho(0)) return (1);

    /* Clean the geometry */

    _vorder.clear();

    return (0);
  }

  Id CalcModelPGS::_varioPgsNoStat()
  {
    Id error = 1;
    bool flag_correl = _ngrf > 1 && (_optCorrel != 2 || _rule->getRho() != 0);

    if (_varioPgsVariable(1, 1, 0)) goto label_end;
    _corpgs.define(_optCorrel, _flagRho, _rule->getRho());
    _tracepgs.define(_ngrf, _corpgs.getNpar(), _flagRho, flag_correl);

    /* Infer the variogram of PGS */

    if (_varioPgsVariable(0, 1, 0)) goto label_end;
    if (!_flagRho)
    {
      setRho(_rule->getRho());
      if (_variopgsCalculNoRho(1)) goto label_end;
    }
    else
    {
      if (_variopgsCalculRho()) goto label_end;
    }

    /* Set the error return flag */

    error = 0;

  label_end:
    _tracepgs.extractTrace();
    (void)_varioPgsVariable(-1, 1, 0);
    return (error);
  }

  bool CalcModelPGS::_checkTestDiscret() const
  {
    // Avoiding Discretizing calculation is always valid
    if (!_useDiscrete) return true;

    // Case where the Discrete calculation option has been switch ON
    if (_rule->getModeRule() == ERule::STD)
    {
      // Only authorized when GRFs are not correlated
      if (_flagRho)
      {
        messerr("Calculations may not be performed using Discretized Version");
        messerr("when underlying GRFs are correlated");
        return false;
      }
    }
    else
    {
      messerr("Calculations may not be performed using Discretized version");
      messerr("when the Rule is not Standard (ERule::STD)");
      return false;
    }
    return true;
  }

  double CalcModelPGS::getProbaInd(
    double correl,
    double low[2],
    double up[2],
    Id iconf,
    Id maxpts)
  {
    double proba = TEST;
    if (!_useDiscrete)
    {
      if (correl == 0.)
      {
        double p1min = law_cdf_gaussian(low[0]);
        double p1max = law_cdf_gaussian(up[0]);
        double p2min = law_cdf_gaussian(low[1]);
        double p2max = law_cdf_gaussian(up[1]);
        proba = (p1max - p1min) * (p2max - p2min);
      }
      else
      {
        double err;
        Id ier;
        Id infin[2];
        infin[0] = mvndst_infin(low[0], up[0]);
        infin[1] = mvndst_infin(low[1], up[1]);
        mvndst(
          2, low, up, infin, &correl, maxpts, _abseps, _releps, &err, &proba,
          &ier);
      }
    }
    else
    {
      proba = _discretepgs->calculateByRank(iconf, low, up);
    }
    return proba;
  }

  /****************************************************************************/
  /*!
   **  Calculate the Gaussian variograms
   **
   ** \return  Error return code
   **
   ** \param[in]  db           Db structure
   ** \param[in]  varioparam   VarioParam structure for the GRFs
   ** \param[in]  ruleprop     RuleProp structure
   ** \param[in]  use_discrete 1 to use the discrete version of the calculation
   ** \param[in]  flag_rho     1 if the correlation coefficient must be regressed
   ** \param[in]  opt_correl   0 full model; 1 symmetrical; 2 residuals
   **
   ** \remarks This is simply a routine dispatching between the stationary function
   ** \remarks and the non-stationary one
   **
   *****************************************************************************/
  Vario* variogram_pgs(
    Db* db,
    const VarioParam* varioparam,
    const RuleProp* ruleprop,
    bool use_discrete,
    Id flag_rho,
    Id opt_correl)
  {
    CalcModelPGS modelpgs(db, varioparam, ruleprop);
    modelpgs.setRunType(1);
    modelpgs.setUseDb(true);
    modelpgs.setFlagRho(flag_rho);
    modelpgs.setOptCorrel(opt_correl);
    modelpgs.setUseDiscrete(use_discrete);

    Id error = (modelpgs.run()) ? 0 : 1;
    if (error) return nullptr;
    return modelpgs.getVario();
  }

  /****************************************************************************/
  /*!
   **  Evaluate the experimental variogram of indicators in PluriGaussian case
   **
   ** \return  Error return code
   **
   ** \param[in]  db         Db descriptor
   ** \param[in]  varioparam VarioParam structure
   ** \param[in]  ruleprop   RuleProp structure
   ** \param[in]  model1     First Model structure
   ** \param[in]  model2     Second Model structure (optional)
   ** \param[in]  use_discrete 1 to use the discrete version of the calculation
   **
   *****************************************************************************/
  Vario* model_pgs(
    Db* db,
    const VarioParam* varioparam,
    const RuleProp* ruleprop,
    const Model* model1,
    const Model* model2,
    bool use_discrete)
  {
    CalcModelPGS modelpgs(db, varioparam, ruleprop);
    modelpgs.setRunType(2);
    modelpgs.setUseDb(true);
    modelpgs.setModel1(model1);
    modelpgs.setModel2(model2);
    modelpgs.setUseDiscrete(use_discrete);

    Id error = (modelpgs.run()) ? 0 : 1;
    if (error) return nullptr;
    return modelpgs.getVario();
  }

} // namespace gstlrn
