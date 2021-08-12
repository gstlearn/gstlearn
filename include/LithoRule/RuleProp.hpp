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

class RuleProp : public AStringable
{
public:
  RuleProp();
  RuleProp(const Rule* rule, const VectorDouble& props = VectorDouble());
  RuleProp(const Rule* rule, const Db* dbprop);
  RuleProp(const RuleProp& m);
  RuleProp& operator=(const RuleProp &m);
  virtual ~RuleProp();

  virtual std::string toString(int level = ITEST) const override;

  const Db* getDbprop() const { return _dbprop; }
  void setDbprop(const Db* dbprop) { _dbprop = dbprop; }
  bool isFlagStat() const { return _flagStat; }
  void setFlagStat(bool flagStat) { _flagStat = flagStat; }
  const VectorDouble& getProps() const { return _props; }
  void setProps(const VectorDouble& props) { _props = props; }
  const Rule* getRule() const { return _rule; }
  void setRule(const Rule*& rule) { _rule = rule; }

private:
  bool _checkConsistency();

private:
  bool _flagStat;
  VectorDouble _props;
  const Db* _dbprop;
  const Rule* _rule;
};
