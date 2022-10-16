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
#include "Basic/OptDbg.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>

std::vector<EDbg> OptDbg::_dbg = std::vector<EDbg>();
int OptDbg::_currentIndex = -1;
int OptDbg::_reference = -1;

void OptDbg::reset()
{
  _dbg.clear();
}

bool OptDbg::query(const EDbg& option)
{
  if (force()) return true;
  for (auto e: _dbg)
  {
    if (e == option) return true;
  }
  return false;
}

bool OptDbg::queryByKey(const String& name)
{
  if (force()) return true;
  auto it = EDbg::getIterator();
  while (it.hasNext())
  {
    EDbg e = *it;
    if (e.getKey() == toUpper(name)) return query(e);
    it.toNext();
  }
  return false;
}

/**
 * Switching ON a option
 * @param option Description of the option (Keyword)
 */
void OptDbg::define(const EDbg& option)
{
  if (!query(option))
  {
    _dbg.push_back(option);
  }
}

/**
 * Switching OFF an option
 * @param option Description of the Option (Keyword)
 */
void OptDbg::undefine(const EDbg& option)
{
  if (query(option))
  {
    _dbg.erase(std::remove(_dbg.begin(), _dbg.end(), option), _dbg.end());
  }
}

/**
 * Switching ON a option
 * @param name Description of the option (Name)
 */
void OptDbg::defineByKey(const String& name)
{
  auto it = EDbg::getIterator();
  while (it.hasNext())
  {
    EDbg e = *it;
    if (e.getKey() == toUpper(name)) define(e);
    it.toNext();
  }
}

/**
 * Switching OFF a option
 * @param name Description of the option (Name)
 */
void OptDbg::undefineByKey(const String& name)
{
  auto it = EDbg::getIterator();
  while (it.hasNext())
  {
    EDbg e = *it;
    if (e.getKey() == toUpper(name)) undefine(e);
    it.toNext();
  }
}

void OptDbg::defineAll(void)
{
  auto it = EDbg::getIterator();
  while (it.hasNext())
  {
    EDbg e = *it;
    define(e);
    it.toNext();
  }
}

void OptDbg::undefineAll(void)
{
  auto it = EDbg::getIterator();
  while (it.hasNext())
  {
    EDbg e = *it;
    undefine(e);
    it.toNext();
  }
}

void OptDbg::display()
{
  std::stringstream sstr;

  sstr << toTitle(1,"List of Debug Options");
  auto it = EDbg::getIterator();
  while (it.hasNext())
  {
    EDbg e = *it;
    sstr << std::setw(30) << e.getDescr() <<
        "[ " << std::setw(9) << e.getKey() << "]" <<
        " : " << query(e) << std::endl;
    it.toNext();
  }

  if (_reference >= 0)
     sstr << "Index of the reference target under DEBUG = " << _reference << std::endl;

  sstr << "Use 'OptCst::define' or 'OptCst::setReference' to modify previous values" << std::endl;

  messageFlush(sstr.str());
}

bool OptDbg::force(void)
{
  if (_reference < 0) return false;
  if (_currentIndex != _reference) return false;
  return true;
}
