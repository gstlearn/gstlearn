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
#include "Basic/EOptCst.hpp"
#include "Basic/OptCst.hpp"
#include "Basic/AStringable.hpp"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <map>

std::map<const ECst, double> OptCst::_cst = std::map<const ECst, double>({
                                       { ECst::NTCAR,  10. },
                                       { ECst::NTDEC,   3. },
                                       { ECst::NTROW,   7. },
                                       { ECst::NTCOL,   7. },
                                       { ECst::NTBATCH, 7. },
                                       { ECst::NTNAME, 12. },
                                       { ECst::NTRANK,  3. },
                                       { ECst::NPROC,   0. },
                                       { ECst::LOCMOD,  0. },
                                       { ECst::LOCNEW,  0. },
                                       { ECst::TOLINV,  1.e-20 },
                                       { ECst::TOLGEN,  1.e-20 },
                                       { ECst::EPSMAT,  2.3e-16 },
                                       { ECst::EPSSVD,  1.e-5 }
});

double OptCst::query(const ECst& option)
{
  for (auto e: _cst)
  {
    if (e.first == option) return e.second;
  }
  return TEST;
}

double OptCst::queryByKey(const String& name)
{
  for (auto e: _cst)
  {
    if (e.first.getKey() == name) return e.second;
  }
  return TEST;
}


void OptCst::define(const ECst& option, double value)
{
  for (auto &e: _cst)
  {
    if (e.first == option)
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
    if (e.first.getKey() == name)
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
    sstr << std::setw(50) << e.first.getDescr() <<
        " [" << std::setw(7) << e.first.getKey() << "]" << " : " <<
        e.second << std::endl;
  }

  sstr << "Use 'define' to modify previous values" << std::endl;

  messageFlush(sstr.str());
}
