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
#include "Basic/OptCustom.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <map>

std::map<const String, double> OptCustom::_cst = std::map<const String, double>();

double OptCustom::query(const String& name, double valdef)
{
  for (auto e: _cst)
  {
    if (e.first == name) return e.second;
  }
  return valdef;
}

void OptCustom::define(const String& name, double value)
{
  // Check if the entry already exists
  for (auto &e: _cst)
  {
    if (e.first == name)
    {
      e.second = value;
      return;
    }
  }

  // Add the entry
  _cst.insert({name,value});
}

void OptCustom::undefine(const String& name)
{
  for (auto &e: _cst)
  {
    if (e.first == name)
    {
      _cst.erase(name);
      return;
    }
  }
}

void OptCustom::display(void)
{
  std::stringstream sstr;

  sstr << toTitle(1,"List of Custom Options");

  for (auto e: _cst)
  {
    sstr << std::setw(50) << e.first << " : " <<
        e.second << std::endl;
  }
  messageFlush(sstr.str());
}
