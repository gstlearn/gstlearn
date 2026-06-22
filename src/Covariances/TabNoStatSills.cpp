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
#include "Covariances/TabNoStatSills.hpp"
#include "Covariances/TabNoStat.hpp"
#include "Enum/EConsElem.hpp"
#include "geoslib_define.h"

namespace gstlrn
{
  TabNoStatSills::TabNoStatSills() {}

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

  bool TabNoStatSills::isDefinedForVariance() const
  {
    return !empty();
  }

  Id TabNoStatSills::getNSills() const
  {
    // If not item is recorded, the returned rank is obviously "0"
    if (size() <= 0) return 0;

    // If some items are recorded, find the number of ones refering to "sills"
    Id count = 0;
    for (const auto& [paramId, noStatPtr]: getTable())
    {
      if (_isValid(paramId.getType())) count++;
    }
    return count;
  }

  TabNoStatSills::~TabNoStatSills() {}
} // namespace gstlrn
