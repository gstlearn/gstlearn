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
#include "LithoRule/RuleProp.hpp"
#include "Basic/AStringable.hpp"
#include "Db/Db.hpp"
#include "LithoRule/Rule.hpp"
#include "PluriGaussian/CalcModelPGS.hpp"
#include "geoslib_f_private.h"
#include <algorithm>
#include <vector>

namespace gstlrn
{
  RuleProp::RuleProp()
    : AStringable()
    , _flagStat(true)
    , _propcst()
    , _dbprop(nullptr)
    , _rules()
  {
  }

  /**
   * This constructor is used in the exceptional case where the Rule is not yet defined
   * (typically when inferring the Rule)
   * @param dbprop  Db containing the Proportion information (ELoc::P fields)
   * @param propcst Vector of constant proportions
   */
  Id RuleProp::resetFromDb(const Db* dbprop, const VectorDouble& propcst)
  {
    _clear();

    _flagStat = true;
    _dbprop = dbprop;
    _propcst = propcst;

    if (!_checkConsistency()) return 1;

    // A generic rule is created on the fly
    auto nfacies = _getNFacies();
    _rules[0] = std::unique_ptr<Rule>(Rule::createFromFaciesCount(nfacies));

    return 0;
  }

  Id RuleProp::resetFromRule(const Rule* rule, const VectorDouble& propcst)
  {
    _clear();

    _flagStat = true;
    _propcst = propcst;

    if (rule != nullptr) _rules[0] = std::unique_ptr<Rule>(rule->clone());
    if (!_checkConsistency()) return 1;
    return 0;
  }

  Id RuleProp::resetFromRuleAndDb(const Rule* rule, const Db* dbprop)
  {
    _clear();

    _flagStat = true;
    _dbprop = dbprop;

    if (rule != nullptr) _rules[0] = std::unique_ptr<Rule>(rule->clone());
    if (!_checkConsistency()) return 1;
    return 0;
  }

  Id RuleProp::resetFromRules(
    const Rule* rule1,
    const Rule* rule2,
    const VectorDouble& propcst)
  {
    _clear();

    _flagStat = true;
    _propcst = propcst;

    if (rule1 != nullptr) _rules[0] = std::unique_ptr<Rule>(rule1->clone());
    if (rule2 != nullptr) _rules[1] = std::unique_ptr<Rule>(rule2->clone());
    if (!_checkConsistency()) return 1;
    return 0;
  }

  Id RuleProp::resetFromRulesAndDb(
    const Rule* rule1,
    const Rule* rule2,
    const Db* dbprop)
  {
    _clear();

    _flagStat = true;
    _dbprop = dbprop;

    if (rule1 != nullptr) _rules[0] = std::unique_ptr<Rule>(rule1->clone());
    if (rule2 != nullptr) _rules[1] = std::unique_ptr<Rule>(rule2->clone());
    if (!_checkConsistency()) return 1;
    return 0;
  }

  void RuleProp::_clear()
  {
    _dbprop = nullptr;
    clearRule();
  }

  String RuleProp::toString(const AStringFormat* strfmt) const
  {
    std::stringstream sstr;
    if (getNRule() <= 0) return sstr.str();

    // Stationary Flag
    if (_flagStat)
      mestitle(0, "RuleProp in Stationary Case");
    else
      mestitle(0, "RuleProp in Non-Stationary Case");

    // Proportions
    if (_flagStat)
      sstr << "- Constant Proportions" << toStrVector(String(), _propcst)
           << std::endl;
    else
      sstr << "- Non-Stationary Proportions are read from Db" << std::endl;
    // if (! _flagStat)
    //   sstr << _dbprop->toString(strfmt);

    // Rules
    for (Id ir = 0; ir < getNRule(); ir++)
    {
      sstr << _rules[ir]->toString(strfmt);
    }

    return sstr.str();
  }

  bool RuleProp::_checkConsistency()
  {
    Id nfacies = 0;

    // Check the number of facies against the Rule(s)
    if (getNRule() > 0)
    {
      // In case of several rules, the number of facies is the product
      // of the number of facies per rule.
      Id nfacrule = 1;
      for (Id ir = 0; ir < getNRule(); ir++)
        nfacrule *= _rules[ir]->getNFacies();
      nfacies = nfacrule;
    }

    // Non-stationary case: proportions are provided using Dbprop
    if (_dbprop != nullptr)
    {
      _flagStat = false;
      _propcst.clear();

      // Check consistency of the number of facies
      Id nfacdb = _dbprop->getNLoc(ELoc::P);
      if (nfacies > 0 && nfacies != nfacdb)
      {
        messerr("Mismatch between:");
        messerr("- Number of Facies in Rule(s) (%d)", nfacies);
        messerr("- Number of Proportion fields in Db (%d)", nfacdb);
        return false;
      }
      return true;
    }

    // Stationary proportions provided by 'propcst'
    if (!_propcst.empty())
    {
      _flagStat = true;
      _dbprop = nullptr;

      // Check consistency of the number of facies
      Id nfacprop = static_cast<Id>(_propcst.size());
      if (nfacies > 0 && nfacies != nfacprop)
      {
        messerr("Mismatch between:");
        messerr("- Number of Facies in Rule(s) (%d)", nfacies);
        messerr("- Number of Proportion in Propcst (%d)", nfacprop);
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
    _propcst.resize(nfacies, 1. / static_cast<double>(nfacies));
    return true;
  }

  bool RuleProp::_isRuleRankValid(Id rank) const
  {
    return checkArg("Rule Rank", rank, getNRule());
  }

  Id RuleProp::_getNFacies()
  {
    // Check the number of facies against the Rule
    if (!_rules.empty())
    {
      Id nfacies = 1;
      for (Id ir = 0; ir < getNRule(); ir++)
        nfacies *= _rules[ir]->getNFacies();
      return nfacies;
    }

    // Non-stationary case: proportions are provided using Dbprop
    if (_dbprop != nullptr)
    {
      return _dbprop->getNLoc(ELoc::P);
    }

    // Stationary proportions provided by 'propcst'
    if (!_propcst.empty())
    {
      return static_cast<Id>(_propcst.size());
    }

    return 0;
  }

  const Rule* RuleProp::getRule(Id rank) const
  {
    if (getNRule() <= 0) return nullptr;
    if (!_isRuleRankValid(rank)) return nullptr;
    return _rules[rank].get();
  }

  void RuleProp::addRule(const Rule& rule)
  {
    _rules[getNRule()] = std::unique_ptr<Rule>(rule.clone());
  }

  void RuleProp::clearRule()
  {
    _rules = {};
  }

  Id RuleProp::getNRule() const
  {
    return std::count_if(
      _rules.begin(), _rules.end(),
      [](const std::unique_ptr<Rule>& r) { return r != nullptr; });
  }

  /**
   * @brief Returns the list of potential Rule arrangements, sorted by incrasing scores
   *
   * @param db Descritpion of the input Data Base
   * @param varioparam Description of the manner to calculate the experimental indicator variograms
   * @param ngrfmax Maximim number of GRFs to be tested in the optimal Rules
   * @param use_discrete When True, use the numerical integration of Gaussian(s)
   * @param verbose Verbose flag
   *
   * @return List of potential Rule arrangements defined by pointers
   */
  std::vector<Rule> RuleProp::fit(
    Db* db,
    const VarioParam* varioparam,
    Id ngrfmax,
    bool use_discrete,
    bool verbose)
  {
    CalcModelPGS ruleauto(db, varioparam, this);
    ruleauto.setRunType(3);
    ruleauto.setUseDb(true);
    ruleauto.setNgrfMax(ngrfmax);
    ruleauto.setUseDiscrete(use_discrete);
    ruleauto.setVerbose(verbose);

    Id error = (ruleauto.run()) ? 0 : 1;

    std::vector<Rule> sortedRules;
    if (!error)
    {
      sortedRules = ruleauto.getSortedRules();
      auto& ruleOpt = sortedRules[0];
      clearRule();
      addRule(ruleOpt);
    }
    return sortedRules;
  }

  /**
   * Convert a set of Gaussian vectors into the corresponding Facies in a Db
   * @param db      Pointer to the Db structure (in/out)
   * @param namconv Naming convention
   * @return Error return code
   * @remarks The input variables must be locatorized Z or SIMU
   */
  Id RuleProp::gaussToCategory(Db* db, const NamingConvention& namconv) const
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
  Id RuleProp::categoryToThresh(Db* db, const NamingConvention& namconv) const
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
  Id RuleProp::computeAllThreshes(Db* db, const NamingConvention& namconv) const
  {
    if (_rules[0]->getModeRule() != ERule::STD)
    {
      messerr("This method is only available for ERule::STD type of Rule");
      return 1;
    }
    return _db_threshold(db, this, nullptr, namconv);
  }

  RuleProp*
    RuleProp::createFromDb(const Db* dbprop, const VectorDouble& propcst)
  {
    auto* ruleprop = new RuleProp;
    if (ruleprop->resetFromDb(dbprop, propcst))
    {
      messerr("Problem when creating from Db");
      delete ruleprop;
      return nullptr;
    }
    return ruleprop;
  }

  RuleProp*
    RuleProp::createFromRule(const Rule* rule, const VectorDouble& propcst)
  {
    auto* ruleprop = new RuleProp;
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
    auto* ruleprop = new RuleProp;
    if (ruleprop->resetFromRuleAndDb(rule, dbprop))
    {
      messerr("Problem when creating from Rule & Db");
      delete ruleprop;
      return nullptr;
    }
    return ruleprop;
  }

  RuleProp* RuleProp::createFromRules(
    const Rule* rule1,
    const Rule* rule2,
    const VectorDouble& propcst)
  {
    auto* ruleprop = new RuleProp;
    if (ruleprop->resetFromRules(rule1, rule2, propcst))
    {
      messerr("Problem when creating from Rules & Proportions");
      delete ruleprop;
      return nullptr;
    }
    return ruleprop;
  }

} // namespace gstlrn
