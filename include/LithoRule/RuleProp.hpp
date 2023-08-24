/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Basic/AStringable.hpp"
#include "LithoRule/Rule.hpp"
#include "Basic/NamingConvention.hpp"

#include <vector>

class Db;
class VarioParam;

class GSTLEARN_EXPORT RuleProp : public AStringable
{
public:
  RuleProp();
  RuleProp(const RuleProp& m);
  RuleProp& operator=(const RuleProp &m);
  virtual ~RuleProp();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int resetFromDb(const Db* dbprop,
                  const VectorDouble& propcst = VectorDouble());
  int resetFromRule(const Rule* rule,
                    const VectorDouble& propcst = VectorDouble());
  int resetFromRuleAndDb(const Rule* rule, const Db* dbprop);
  int resetFromRules(const Rule* rule1,
                     const Rule* rule2,
                     const VectorDouble& propcst = VectorDouble());
  int resetFromRulesAndDb(const Rule* rule1,
                          const Rule* rule2,
                          const Db* dbprop);

  static RuleProp* createFromDb(const Db* dbprop,
                                const VectorDouble& propcst = VectorDouble());
  static RuleProp* createFromRule(const Rule* rule,
                                  const VectorDouble& propcst = VectorDouble());
  static RuleProp* createFromRuleAndDb(const Rule* rule, const Db* dbprop);
  static RuleProp* createFromRules(const Rule* rule1,
                                   const Rule* rule2,
                                   const VectorDouble& propcst = VectorDouble());
  static RuleProp* createFromRulesAndDb(const Rule* rule1,
                                        const Rule* rule2,
                                        const Db* dbprop);

  const Db* getDbprop() const { return _dbprop; }
  void setDbprop(const Db* dbprop) { _dbprop = dbprop; }
  bool isFlagStat() const { return _flagStat; }
  void setFlagStat(bool flagStat) { _flagStat = flagStat; }
  const VectorDouble& getPropCst() const { return _propcst; }
  void setPropCst(const VectorDouble& propcst) { _propcst = propcst; }
  const Rule* getRule(int rank = 0) const;
  void addRule(const Rule* rule);
  void clearRule();
  int getRuleNumber() const { return static_cast<int>(_rules.size()); }

  int fit(Db* db,
          const VarioParam* varioparam,
          int ngrfmax = 1,
          bool verbose = false);
  int gaussToCategory(Db *db,
                      const NamingConvention &namconv = NamingConvention(
                          "Facies", true, true, true, ELoc::fromKey("FACIES"))) const;
  int categoryToThresh(Db *db, const NamingConvention& namconv = NamingConvention("Bounds")) const;
  int computeAllThreshes(Db *db, const NamingConvention& namconv = NamingConvention("Thresh")) const;

private:
  void _clear();
  bool _checkConsistency();
  bool _checkRuleRank(int rank) const;
  int _getNFacies();

private:
  bool _flagStat;
  VectorDouble _propcst;
  const Db* _dbprop;
  std::vector<const Rule*> _rules;
  bool _ruleInternal; // TRUE if a fictitious rule has been established internally
};
