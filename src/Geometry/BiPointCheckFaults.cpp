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

#include "Geometry/BiPointCheckFaults.hpp"
#include "Space/SpacePoint.hpp"

BiPointCheckFaults::BiPointCheckFaults(const Faults* faults)
    : ABiPointCheck(),
      _faults(faults)
{
}

BiPointCheckFaults::BiPointCheckFaults(const BiPointCheckFaults &r)
    : ABiPointCheck(r),
      _faults(r._faults)
{
}

BiPointCheckFaults& BiPointCheckFaults::operator=(const BiPointCheckFaults &r)
{
  if (this != &r)
  {
    ABiPointCheck::operator=(r);
    _faults = r._faults;
  }
  return *this;
}

BiPointCheckFaults::~BiPointCheckFaults()
{
}

String BiPointCheckFaults::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << "- Separated by Faults" << std::endl;

  return sstr.str();
}

bool BiPointCheckFaults::isOK(const SpacePoint &P1,
                              const SpacePoint &P2,
                              int iech1,
                              int iech2) const
{
  return _faults->isSplitByFaultSP(P1, P2);
}
