/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f.h"

#include "Geometry/BiPointCheckDistance.hpp"
#include "Space/SpacePoint.hpp"

BiPointCheckDistance::BiPointCheckDistance(double dmax)
    : ABiPointCheck(),
      _dmax(dmax)
{
}

BiPointCheckDistance::BiPointCheckDistance(const BiPointCheckDistance &r)
    : ABiPointCheck(r),
      _dmax(r._dmax)
{
}

BiPointCheckDistance& BiPointCheckDistance::operator=(const BiPointCheckDistance &r)
{
  if (this != &r)
  {
    ABiPointCheck::operator=(r);
    _dmax = r._dmax;
  }
  return *this;
}

BiPointCheckDistance::~BiPointCheckDistance()
{
}

String BiPointCheckDistance::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << ABiPointCheck::toString();

  sstr << "- Based on distance" << _dmax << std::endl;

  return sstr.str();
}

bool BiPointCheckDistance::isOK(const SpacePoint& P1, const SpacePoint& P2) const
{
  _dist = P1.getDistance(P2);
  return _dist <= _dmax;
}
