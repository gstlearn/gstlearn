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
#include "LithoRule/RuleProp.hpp"
#include "LithoRule/Rule.hpp"
#include "Db/Db.hpp"
#include "Basic/AException.hpp"
#include "Basic/AStringable.hpp"
#include "geoslib_f_private.h"

RuleProp::RuleProp()
    : _flagStat(true),
      _propcst(),
      _dbprop(nullptr),
      _rules(),
      _ruleInternal(false)
{
}

/**
 * This constructor is used in the exceptional case where the Rule is not yet defined
 * (typically when inferring the Rule)
 * @param dbprop  Db containing the Proportion information (LOC_P fields)
 * @param propcst Vector of constant proportions
 */
RuleProp::RuleProp(const Db* dbprop, const VectorDouble& propcst)
    : _flagStat(true),
      _propcst(propcst),
      _dbprop(dbprop),
      _rules(),
      _ruleInternal(true)
{
  if (! _checkConsistency())
    my_throw("Inconsistent arguments");

  // A generic rule is created on the fly
  int nfacies = _getNFacies();
  Rule rule = Rule(nfacies);
  _rules.push_back(new Rule(rule));
}

RuleProp::RuleProp(const Rule* rule, const VectorDouble& propcst)
    : _flagStat(true),
      _propcst(propcst),
      _dbprop(nullptr),
      _rules(),
      _ruleInternal(false)
{
  if (rule != (Rule *) NULL) _rules.push_back(rule);
  if (! _checkConsistency())
    my_throw("Inconsistent arguments");
}

RuleProp::RuleProp(const Rule* rule, const Db* dbprop)
    : _flagStat(true),
      _propcst(),
      _dbprop(dbprop),
      _rules(),
      _ruleInternal(false)
{
  _rules.push_back(rule);
  if (! _checkConsistency())
    my_throw("Inconsistent arguments");
}

RuleProp::RuleProp(const Rule* rule1,
                   const Rule* rule2,
                   const VectorDouble& propcst)
    : _flagStat(true),
      _propcst(propcst),
      _dbprop(nullptr),
      _rules(),
      _ruleInternal(false)
{
  _rules.push_back(rule1);
  _rules.push_back(rule2);
  if (! _checkConsistency())
    my_throw("Inconsistent arguments");
}

RuleProp::RuleProp(const Rule* rule1, const Rule* rule2, const Db* dbprop)
    : _flagStat(true),
      _propcst(),
      _dbprop(dbprop),
      _rules(),
      _ruleInternal(false)
{
  _rules.push_back(rule1);
  _rules.push_back(rule2);
  if (!_checkConsistency())
  my_throw("Inconsistent arguments");
}

RuleProp::RuleProp(const RuleProp& m)
  : _flagStat(m._flagStat),
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
  if (_ruleInternal)
  {
    for (int ir = 0; ir < getRuleNumber(); ir++)
      delete _rules[ir];
  }
}

String RuleProp::toString(int level) const
{
  std::stringstream sstr;

  // Stationary Flag
  if (_flagStat)
    mestitle(0,"RuleProp in Stationary Case");
  else
    mestitle(0,"RuleProp in Non-Stationary Case");

  // Constant proportions (Stationary case)
  if (_flagStat)
    sstr << "Constant Proportions" << ut_vector_string(_propcst) << std::endl;

  // Db file (Non-Stationary case)
  if (! _flagStat)
    sstr << _dbprop->toString(level);

  // Rules
  for (int ir = 0; ir < getRuleNumber(); ir++)
  {
    const Rule* rule = _rules[ir];
    sstr << rule->toString(level);
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
    int nfacdb = _dbprop->getFromLocatorNumber(LOC_P);
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
    int nfacprop = _propcst.size();
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
    return _dbprop->getFromLocatorNumber(LOC_P);
  }

  // Stationary proportions provided by 'propcst'
  if (! _propcst.empty())
  {
    return _propcst.size();
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
  Rule* ruleFit = rule_auto(db,varioparam,this,ngrfmax,verbose);
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
int RuleProp::gaussToCategory(Db* db, NamingConvention namconv) const
{
  if (_rules[0]->getModeRule() != ERule::STD)
  {
    messerr("This method is only available for ERule::STD type of Rule");
    return 1;
  }
  return db_rule(db, this, nullptr, namconv);
}

/**
 * Derive the bounds variables for a Db (depending on the Category information of each sample)
 * @param db      Pointer to the Db structure (in/out)
 * @param namconv Naming convention
 * @return Error return code
 */
int RuleProp::categoryToThresh(Db *db, NamingConvention namconv) const
{
  if (_rules[0]->getModeRule() != ERule::STD)
  {
    messerr("This method is only available for ERule::STD type of Rule");
    return 1;
  }
  return db_bounds(db, this, nullptr, namconv);
}

/**
 * Calculate all the thresholds at each sample of a Db
 * @param db      Pointer to the Db structure (in/out)
 * @param namconv Naming convention
 * @return Error return code
 */
int RuleProp::computeAllThreshes(Db *db, NamingConvention namconv) const
{
  if (_rules[0]->getModeRule() != ERule::STD)
  {
    messerr("This method is only available for ERule::STD type of Rule");
    return 1;
  }
  return db_threshold(db, this, nullptr, namconv);
}
