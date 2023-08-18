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
#include <Geometry/BiTargetCheckDate.hpp>
#include "geoslib_f.h"

#include "Space/SpaceTarget.hpp"

BiTargetCheckDate::BiTargetCheckDate(double deltamin, double deltamax)
    : ABiTargetCheck(),
      _deltaMin(deltamin),
      _deltaMax(deltamax)
{
}

BiTargetCheckDate::BiTargetCheckDate(const BiTargetCheckDate &r)
    : ABiTargetCheck(r),
      _deltaMin(r._deltaMin),
      _deltaMax(r._deltaMax)
{
}

BiTargetCheckDate& BiTargetCheckDate::operator=(const BiTargetCheckDate &r)
{
  if (this != &r)
  {
    ABiTargetCheck::operator=(r);
    _deltaMin = r._deltaMin;
    _deltaMax = r._deltaMax;
  }
  return *this;
}

BiTargetCheckDate::~BiTargetCheckDate()
{
}

String BiTargetCheckDate::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << "- Date difference must lie between " << _deltaMin << " and " << _deltaMax << std::endl;

  return sstr.str();
}

bool BiTargetCheckDate::isOK(const SpaceTarget &T1, const SpaceTarget &T2) const
{
  double date1 = T1.getDate();
  double date2 = T2.getDate();

  if (FFFF(date1) || FFFF(date2)) return false;

  double delta = date2 - date1;
  if (delta <  _deltaMin) return false;
  if (delta >= _deltaMax) return false;
  return true;
}
