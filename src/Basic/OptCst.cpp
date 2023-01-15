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
#include "Enum/ECst.hpp"

#include "Basic/OptCst.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"
#include "Basic/Utilities.hpp"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <map>

/* This doesn't work with MinGw / G++ 8.3
std::map<const ECst, double> OptCst::_cst =
 {{ ECst::NTCAR,  10. },
  { ECst::NTDEC,   3. },
  { ECst::NTROW,   7. },
  { ECst::NTCOL,   7. },
  { ECst::NTBATCH, 7. },
  { ECst::NTNAME, 12. },
  { ECst::NTRANK,  3. },
  { ECst::NPROC,   0. };
*/
std::map<int, double> OptCst::_cst =
 {{ ECst::E_NTCAR,  10. },
  { ECst::E_NTDEC,   3. },
  { ECst::E_NTROW,   7. },
  { ECst::E_NTCOL,   7. },
  { ECst::E_NTBATCH, 7. },
  { ECst::E_NTNAME, 12. },
  { ECst::E_NTRANK,  3. },
  { ECst::E_NPROC,   0. }};

double OptCst::query(const ECst& option)
{
  for (auto e: _cst)
  {
    if (e.first == option.getValue()) return e.second;
  }
  return TEST;
}

double OptCst::queryByKey(const String& name)
{
  for (auto e: _cst)
  {
    if (ECst::fromValue(e.first).getKey() == toUpper(name))
      return e.second;
  }
  return TEST;
}


void OptCst::define(const ECst& option, double value)
{
  for (auto &e: _cst)
  {
    if (e.first == option.getValue())
    {
      e.second = value;
      return;
    }
  }
}

void OptCst::defineByKey(const String& name, double value)
{
  for (auto &e: _cst)
  {
    if (ECst::fromValue(e.first).getKey() == toUpper(name))
    {
      e.second = value;
      return;
    }
  }
}

void OptCst::display()
{
  std::stringstream sstr;

  sstr << toTitle(1,"List of Options of internal Constant values");

  for (auto e: _cst)
  {
    ECst ec = ECst::fromValue(e.first);
    sstr << std::setw(50) << ec.getDescr() <<
        " [" << std::setw(7) << ec.getKey() << "]" << " : " <<
        e.second << std::endl;
  }

  sstr << "Use 'OptCst::define' to modify previous values" << std::endl;

  messageFlush(sstr.str());
}
