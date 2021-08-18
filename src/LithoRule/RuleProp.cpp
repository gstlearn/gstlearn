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
#include "geoslib_f_private.h"

RuleProp::RuleProp()
    : _flagStat(true),
      _propcst(),
      _dbprop(nullptr),
      _rule(nullptr),
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
      _rule(nullptr),
      _ruleInternal(true)
{
  if (! _checkConsistency())
    my_throw("Inconsistent arguments");

  // A generic rule is created on the fly
  int nfacies = _getNFacies();
  Rule rule = Rule();
  rule.init(nfacies);
  _rule = new Rule(rule);
}

RuleProp::RuleProp(const Rule* rule, const VectorDouble& propcst)
    : _flagStat(true),
      _propcst(propcst),
      _dbprop(nullptr),
      _rule(rule),
      _ruleInternal(false)
{
  if (! _checkConsistency())
    my_throw("Inconsistent arguments");
}

RuleProp::RuleProp(const Rule* rule, const Db* dbprop)
    : _flagStat(true),
      _propcst(),
      _dbprop(dbprop),
      _rule(rule),
      _ruleInternal(false)
{
  if (! _checkConsistency())
    my_throw("Inconsistent arguments");
}

RuleProp::RuleProp(const RuleProp& m)
  : _flagStat(m._flagStat),
    _propcst(m._propcst),
    _dbprop(m._dbprop),
    _rule(m._rule),
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
    _rule = m._rule;
    _ruleInternal = m._ruleInternal;
  }
  return *this;
}

RuleProp::~RuleProp()
{
  if (_ruleInternal) delete _rule;
}

std::string RuleProp::toString(int level) const
{
  std::stringstream sstr;
  sstr << "Class PGS" << std::endl;

  if (_rule != nullptr)
    sstr << _rule->toString(level);

  return sstr.str();
}

bool RuleProp::_checkConsistency()
{
  int nfacies = 0;

  // Check the number of facies against the Rule
  if (_rule != nullptr)
  {
    int nfacrule = _rule->getFaciesNumber();
    if (nfacies > 0 && nfacrule != nfacies)
    {
      messerr("Mismatch between:");
      message("- Number of facies passed as argument (%d)", nfacies);
      messerr("- Number of Facies in Rule (%d)",nfacrule);
      return false;
    }
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
      messerr("- Number of Facies in Rule (%d)",nfacies);
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
      messerr("- Number of Facies in Rule (%d)",nfacies);
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

int RuleProp::_getNFacies()
{
  // Check the number of facies against the Rule
  if (_rule != nullptr)
  {
    return _rule->getFaciesNumber();
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

int RuleProp::fit(Db* db, Vario* vario, int ngrfmax, bool verbose)
{
  Rule* ruleFit = rule_auto(db,vario,this,ngrfmax,verbose);
  if (ruleFit == nullptr) return 1;
  setRule(ruleFit);
  return 0;
}

/**
 * Convert a set of Gaussian vectors into the corresponding Facies in a Db
 * @param db      Pointer to the Db structure (in/out)
 * @param namconv Naming convention
 * @return Error return code
 * @remarks The input variables must be locatorized Z or SIMU
 */
int RuleProp::gaussToCategory(Db* db, NamingConvention namconv)
{
  if (_rule->getModeRule() != RULE_STD)
  {
    messerr("This method is only available for RULE_STD type of Rule");
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
int RuleProp::categoryToThresh(Db *db, NamingConvention namconv)
{
  if (_rule->getModeRule() != RULE_STD)
  {
    messerr("This method is only available for RULE_STD type of Rule");
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
int RuleProp::computeAllThreshs(Db *db, NamingConvention namconv)
{
  if (_rule->getModeRule() != RULE_STD)
  {
    messerr("This method is only available for RULE_STD type of Rule");
    return 1;
  }
  return db_threshold(db, this, nullptr, namconv);
}
