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
#include "Db/Db.hpp"
#include "Basic/AException.hpp"

RuleProp::RuleProp()
    : _flagStat(true),
      _props(),
      _dbprop(nullptr),
      _rule(nullptr)
{
}

RuleProp::RuleProp(const Rule* rule, const VectorDouble& props)
    : _flagStat(true),
      _props(props),
      _dbprop(nullptr),
      _rule(rule)
{
  if (_checkConsistency())
    my_throw("Inconsistent arguments");
}

RuleProp::RuleProp(const Rule* rule, const Db* dbprop)
    : _flagStat(true),
      _props(),
      _dbprop(dbprop),
      _rule(rule)
{
  if (_checkConsistency())
    my_throw("Inconsistent arguments");
}

RuleProp::RuleProp(const RuleProp& m)
  : _flagStat(m._flagStat),
    _props(m._props),
    _dbprop(m._dbprop),
    _rule(m._rule)
{
}

RuleProp& RuleProp::operator=(const RuleProp& m)
{
  if (this != &m)
  {
    _flagStat = m._flagStat;
    _props = m._props;
    _dbprop = m._dbprop;
    _rule = m._rule;
  }
  return *this;
}

RuleProp::~RuleProp()
{
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
  if (_rule == nullptr) return false;
  int nfacrule = _rule->getFaciesNumber();

  // Non-stationary case: proportions are provided using Dbprop
  if (_dbprop != nullptr)
  {
    _flagStat = false;
    _props.clear();

    // Check consistency of the number of facies
    int nfacdb = _dbprop->getFromLocatorNumber(LOC_P);
    if (nfacrule != nfacdb)
    {
      messerr("Mismatch between:");
      messerr("- Number of Facies in Rule (%d)",nfacrule);
      messerr("- Number of Proportion fields in Db (%d)",nfacdb);
      return false;
    }
    return true;
  }

  // Stationary proportions provided by 'props'
  if (! _props.empty())
  {
    _flagStat = true;
    _dbprop = nullptr;

    // Check consistency of the number of facies
    int nfacprop = _props.size();
    if (nfacrule != nfacprop)
    {
      messerr("Mismatch between:");
      messerr("- Number of Facies in Rule (%d)",nfacrule);
      messerr("- Number of Proportion in 'props' (%d)",nfacprop);
      return false;
    }
    return true;
  }

  // Stationary case with proportions not provided
  _flagStat = true;
  _dbprop = nullptr;
  _props.resize(nfacrule, 1./(double) nfacrule);
  return true;
}
