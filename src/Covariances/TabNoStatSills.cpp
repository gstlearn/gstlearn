#include "Covariances/TabNoStatSills.hpp"
#include "Covariances/ParamId.hpp"
#include "Covariances/TabNoStat.hpp"
#include "Enum/EConsElem.hpp"
#include "geoslib_define.h"


TabNoStatSills::TabNoStatSills()
{
}

TabNoStatSills::TabNoStatSills(const TabNoStatSills& m)
  : TabNoStat(m)
{

}

TabNoStatSills& TabNoStatSills::operator=(const TabNoStatSills& m)
{
  if (this != &m)
  {
    TabNoStat::operator=(m);
  }
  return *this;
}

bool TabNoStatSills::_isValid(const EConsElem& econs) const
{
  return (econs == EConsElem::SILL);
}

String TabNoStatSills::toString(const AStringFormat* strfmt) const
{
  return toStringInside(strfmt, 0);
}
String TabNoStatSills::toStringInside(const AStringFormat* strfmt, int i) const
{
  std::stringstream sstr;
  if (empty()) return sstr.str();

  for (const auto& e: getTable())
  {
    sstr << std::to_string(i + 1) << " - ";
    sstr << e.first.toString(strfmt);
    sstr << e.second->toString(strfmt);
    i++;
  }
  return sstr.str();
}

bool TabNoStatSills::isDefinedForVariance() const
{
 return !empty();
}

int TabNoStatSills::getNSills() const
{
  return size();
}

TabNoStatSills::~TabNoStatSills()
{
}