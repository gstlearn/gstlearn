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

/**
 * Check if a DEBUG option has been switched ON or NOT
 * @param option Type of the option to be searched for
 * @param discardForce When TRUE, does not consider the "forcing" option (see remarks)
 *
 * @remarks The use of gstlearn enables the user to FORCE the switch ON of one or several
 * debugging options. This is the case for example when you want to trace the processing
 * of a specific target (defined using setReference).
 * When this option is switched OFF, this specific case is not taken into account
 * @return TRUE if the option is switch ON, FALSE otherwise
 */
bool OptDbg::query(const EDbg& option, bool discardForce)
{
  DECLARE_UNUSED(discardForce);

  if (force()) return true;
#ifdef USE_BOOST_SPAN
  return std::any_of(_dbg.cbegin(), _dbg.cend(), [&option](const auto& e) { return e == option; });
#else
  return std::ranges::any_of(_dbg, [&option](const auto& e) { return e == option; });
#endif
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
        " : " << query(e, true) << std::endl;
    it.toNext();
  }

  if (_reference >= 0)
     sstr << "Index of the reference target under DEBUG = " << _reference << std::endl;

  sstr << "Use 'OptDbg::define' to modify the previous values" << std::endl;
  sstr << "Use 'OptDbg::setReference' to define the target index where all flags are turned ON" << std::endl;

  messageFlush(sstr.str());
}

bool OptDbg::force(void)
{
  if (_reference < 0) return false;
  if (_currentIndex != _reference) return false;
  return true;
}
