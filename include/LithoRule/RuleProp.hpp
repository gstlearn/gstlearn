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

#include "Basic/AStringable.hpp"
#include "LithoRule/Rule.hpp"
#include "Basic/NamingConvention.hpp"

class Db;
class Vario;

class RuleProp : public AStringable
{
public:
  RuleProp();
  RuleProp(const Db* dbprop, const VectorDouble& propcst);
  RuleProp(const Rule* rule, const VectorDouble& propcst = VectorDouble());
  RuleProp(const Rule* rule, const Db* dbprop);
  RuleProp(const RuleProp& m);
  RuleProp& operator=(const RuleProp &m);
  virtual ~RuleProp();

  virtual std::string toString(int level = ITEST) const override;

  const Db* getDbprop() const { return _dbprop; }
  void setDbprop(const Db* dbprop) { _dbprop = dbprop; }
  bool isFlagStat() const { return _flagStat; }
  void setFlagStat(bool flagStat) { _flagStat = flagStat; }
  const VectorDouble& getPropCst() const { return _propcst; }
  void setPropCst(const VectorDouble& propcst) { _propcst = propcst; }
  const Rule* getRule() const { return _rule; }
  void setRule(const Rule* rule) { _rule = rule; }

  int fit(Db* db, Vario* vario, int ngrfmax = 1, bool verbose = false);
  int gaussToCategory(Db* db, NamingConvention namconv = NamingConvention("Facies",LOC_FACIES)) const;
  int categoryToThresh(Db *db, NamingConvention namconv = NamingConvention("Bounds")) const;
  int computeAllThreshes(Db *db, NamingConvention namconv = NamingConvention("Thresh")) const;

private:
  bool _checkConsistency();
  int _getNFacies();

private:
  bool _flagStat;
  VectorDouble _propcst;
  const Db* _dbprop;
  const Rule* _rule;
  bool _ruleInternal; // TRUE if a fictitious rule has been established internally
};
