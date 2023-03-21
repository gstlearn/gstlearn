/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "LithoRule/RuleProp.hpp"
#include "LithoRule/Rule.hpp"
#include "Db/Db.hpp"
#include "Basic/AException.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorHelper.hpp"
#include "geoslib_f_private.h"

RuleProp::RuleProp()
    : AStringable(),
      _flagStat(true),
      _propcst(),
      _dbprop(nullptr),
      _rules(),
      _ruleInternal(false)
{
}

RuleProp::RuleProp(const RuleProp& m)
  : AStringable(m),
    _flagStat(m._flagStat),
    _propcst(m._propcst),
    _dbprop(m._dbprop),
    _rules(m._rules),
    _ruleInternal(m._ruleInternal)
{
}

RuleProp& RuleProp::operator=(const RuleProp& m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _flagStat = m._flagStat;
    _propcst = m._propcst;
    _dbprop = m._dbprop;
    _rules = m._rules;
    _ruleInternal = m._ruleInternal;
  }
  return *this;
}

RuleProp::~RuleProp()
{
  _clear();
}

/**
 * This constructor is used in the exceptional case where the Rule is not yet defined
 * (typically when inferring the Rule)
 * @param dbprop  Db containing the Proportion information (ELoc::P fields)
 * @param propcst Vector of constant proportions
 */
int RuleProp::resetFromDb(const Db* dbprop, const VectorDouble& propcst)
{
  _clear();

  _flagStat = true;
  _dbprop = dbprop;
  _propcst = propcst;
  _ruleInternal = true;

  if (! _checkConsistency()) return 1;

  // A generic rule is created on the fly
  int nfacies = _getNFacies();
  _rules.push_back(Rule::createFromFaciesCount(nfacies));

  return 0;
}

int RuleProp::resetFromRule(const Rule* rule, const VectorDouble& propcst)
{
  _clear();

  _flagStat = true;
  _propcst = propcst;
  _ruleInternal = false;

  if (rule != nullptr) _rules.push_back(rule);
  if (! _checkConsistency()) return 1;
  return 0;
}

int RuleProp::resetFromRuleAndDb(const Rule* rule, const Db* dbprop)
{
  _clear();

  _flagStat = true;
  _dbprop = dbprop;
  _ruleInternal = false;

  _rules.push_back(rule);
  if (! _checkConsistency()) return 1;
  return 0;
}

int RuleProp::resetFromRules(const Rule* rule1,
                             const Rule* rule2,
                             const VectorDouble& propcst)
{
  _clear();

  _flagStat = true;
  _propcst = propcst;
  _ruleInternal = false;

  _rules.push_back(rule1);
  _rules.push_back(rule2);
  if (! _checkConsistency()) return 1;
  return 0;
}

int RuleProp::resetFromRulesAndDb(const Rule* rule1, const Rule* rule2, const Db* dbprop)
{
  _clear();

  _flagStat = true;
  _dbprop = dbprop;
  _ruleInternal = false;

  _rules.push_back(rule1);
  _rules.push_back(rule2);
  if (!_checkConsistency()) return 1;
  return 0;
}

void RuleProp::_clear()
{
  _dbprop = nullptr;
  if (_ruleInternal)
  {
    for (int ir = 0; ir < getRuleNumber(); ir++)
      delete _rules[ir];
  }
}

String RuleProp::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  if (getRuleNumber() <= 0) return sstr.str();

  // Stationary Flag
  if (_flagStat)
    mestitle(0,"RuleProp in Stationary Case");
  else
    mestitle(0,"RuleProp in Non-Stationary Case");

  // Constant proportions (Stationary case)
  if (_flagStat)
    sstr << "Constant Proportions" << VH::toString(_propcst) << std::endl;

  // Db file (Non-Stationary case)
  if (! _flagStat)
    sstr << _dbprop->toString(strfmt);

  // Rules
  for (int ir = 0; ir < getRuleNumber(); ir++)
  {
    const Rule* rule = _rules[ir];
    sstr << rule->toString(strfmt);
  }

  return sstr.str();
}

bool RuleProp::_checkConsistency()
{
  int nfacies = 0;

  // Check the number of facies against the Rule(s)
  if (getRuleNumber() > 0)
  {
    // In case of several rules, the number of facies is the product
    // of the number of facies per rule.
    int nfacrule = 1;
    for (int ir = 0; ir < getRuleNumber(); ir++)
      nfacrule *= _rules[ir]->getFaciesNumber();
    nfacies = nfacrule;
  }

  // Non-stationary case: proportions are provided using Dbprop
  if (_dbprop != nullptr)
  {
    _flagStat = false;
    _propcst.clear();

    // Check consistency of the number of facies
    int nfacdb = _dbprop->getFromLocatorNumber(ELoc::P);
    if (nfacies > 0 && nfacies != nfacdb)
    {
      messerr("Mismatch between:");
      messerr("- Number of Facies in Rule(s) (%d)",nfacies);
      messerr("- Number of Proportion fields in Db (%d)",nfacdb);
      return false;
    }
    return true;
  }

  // Stationary proportions provided by 'propcst'
  if (! _propcst.empty())
  {
    _flagStat = true;
    _dbprop = nullptr;

    // Check consistency of the number of facies
    int nfacprop = static_cast<int>(_propcst.size());
    if (nfacies > 0 && nfacies != nfacprop)
    {
      messerr("Mismatch between:");
      messerr("- Number of Facies in Rule(s) (%d)",nfacies);
      messerr("- Number of Proportion in Propcst (%d)",nfacprop);
      return false;
    }
    return true;
  }

  // Stationary case with proportions not provided
  if (nfacies <= 0)
  {
    messerr("No solution to determine the number of Facies");
    return false;
  }
  _flagStat = true;
  _dbprop = nullptr;
  _propcst.resize(nfacies, 1./(double) nfacies);
  return true;
}

bool RuleProp::_checkRuleRank(int rank) const
{
  int nrule = getRuleNumber();
  if (rank < 0 || rank >= nrule)
  {
    mesArg("Rule Rank",rank,nrule);
    return false;
  }
  return true;
}

int RuleProp::_getNFacies()
{
  // Check the number of facies against the Rule
  if (! _rules.empty())
  {
    int nfacies = 1;
    for (int ir = 0; ir < getRuleNumber(); ir++)
      nfacies *= _rules[ir]->getFaciesNumber();
    return nfacies;
  }

  // Non-stationary case: proportions are provided using Dbprop
  if (_dbprop != nullptr)
  {
    return _dbprop->getFromLocatorNumber(ELoc::P);
  }

  // Stationary proportions provided by 'propcst'
  if (! _propcst.empty())
  {
    return static_cast<int>(_propcst.size());
  }

  return 0;
}
const Rule* RuleProp::getRule(int rank) const
{
  if (! _checkRuleRank(rank)) return nullptr;
  return _rules[rank];
}

void RuleProp::addRule(const Rule* rule)
{
  _rules.push_back(rule);
}

void RuleProp::clearRule()
{
  _rules.clear();
}

int RuleProp::fit(Db* db, const VarioParam* varioparam, int ngrfmax, bool verbose)
{
  Rule* ruleFit = _rule_auto(db,varioparam,this,ngrfmax,verbose);
  if (ruleFit == nullptr) return 1;
  clearRule();
  addRule(ruleFit);
  return 0;
}

/**
 * Convert a set of Gaussian vectors into the corresponding Facies in a Db
 * @param db      Pointer to the Db structure (in/out)
 * @param namconv Naming convention
 * @return Error return code
 * @remarks The input variables must be locatorized Z or SIMU
 */
int RuleProp::gaussToCategory(Db* db, const NamingConvention& namconv) const
{
  if (_rules[0]->getModeRule() != ERule::STD)
  {
    messerr("This method is only available for ERule::STD type of Rule");
    return 1;
  }
  return _db_rule(db, this, nullptr, namconv);
}

/**
 * Derive the bounds variables for a Db (depending on the Category information of each sample)
 * @param db      Pointer to the Db structure (in/out)
 * @param namconv Naming convention
 * @return Error return code
 */
int RuleProp::categoryToThresh(Db *db, const NamingConvention& namconv) const
{
  if (_rules[0]->getModeRule() != ERule::STD)
  {
    messerr("This method is only available for ERule::STD type of Rule");
    return 1;
  }
  return _db_bounds(db, this, nullptr, namconv);
}

/**
 * Calculate all the thresholds at each sample of a Db
 * @param db      Pointer to the Db structure (in/out)
 * @param namconv Naming convention
 * @return Error return code
 */
int RuleProp::computeAllThreshes(Db *db, const NamingConvention& namconv) const
{
  if (_rules[0]->getModeRule() != ERule::STD)
  {
    messerr("This method is only available for ERule::STD type of Rule");
    return 1;
  }
  return _db_threshold(db, this, nullptr, namconv);
}

RuleProp* RuleProp::createFromDb(const Db* dbprop, const VectorDouble& propcst)
{
  RuleProp* ruleprop = new RuleProp;
  if (ruleprop->resetFromDb(dbprop, propcst))
  {
    messerr("Problem when creating from Db");
    delete ruleprop;
    return nullptr;
  }
  return ruleprop;
}
RuleProp* RuleProp::createFromRule(const Rule* rule,
                                                 const VectorDouble& propcst)
{
  RuleProp* ruleprop = new RuleProp;
  if (ruleprop->resetFromRule(rule, propcst))
  {
    messerr("Problem when creating from Rule & Proportions");
    delete ruleprop;
    return nullptr;
  }
  return ruleprop;
}
RuleProp* RuleProp::createFromRuleAndDb(const Rule* rule, const Db* dbprop)
{
  RuleProp* ruleprop = new RuleProp;
  if (ruleprop->resetFromRuleAndDb(rule, dbprop))
  {
    messerr("Problem when creating from Rule & Db");
    delete ruleprop;
    return nullptr;
  }
  return ruleprop;
}
RuleProp* RuleProp::createFromRules(const Rule* rule1,
                                    const Rule* rule2,
                                    const VectorDouble& propcst)
{
  RuleProp* ruleprop = new RuleProp;
  if (ruleprop->resetFromRules(rule1, rule2, propcst))
  {
    messerr("Problem when creating from Rules & Proportions");
    delete ruleprop;
    return nullptr;
  }
  return ruleprop;
}
RuleProp* RuleProp::createFromRulesAndDb(const Rule* rule1,
                                         const Rule* rule2,
                                         const Db* dbprop)
{
  RuleProp* ruleprop = new RuleProp;
  if (ruleprop->resetFromRulesAndDb(rule1, rule2, dbprop))
  {
    messerr("Problem when creating from Rules & Proportions");
    delete ruleprop;
    return nullptr;
  }
  return ruleprop;
}
