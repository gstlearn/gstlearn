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
#include <Geometry/BiTargetCheckFaults.hpp>
#include "geoslib_f.h"

#include "Space/SpaceTarget.hpp"

BiTargetCheckFaults::BiTargetCheckFaults(const Faults* faults)
    : ABiTargetCheck(),
      _faults(faults)
{
}

BiTargetCheckFaults::BiTargetCheckFaults(const BiTargetCheckFaults &r)
    : ABiTargetCheck(r),
      _faults(r._faults)
{
}

BiTargetCheckFaults& BiTargetCheckFaults::operator=(const BiTargetCheckFaults &r)
{
  if (this != &r)
  {
    ABiTargetCheck::operator=(r);
    _faults = r._faults;
  }
  return *this;
}

BiTargetCheckFaults::~BiTargetCheckFaults()
{
}

BiTargetCheckFaults* BiTargetCheckFaults::create(const Faults* faults)
{
  return new BiTargetCheckFaults(faults);
}

String BiTargetCheckFaults::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  if (_faults != nullptr)
    sstr << "- Separated by Faults" << std::endl;

  return sstr.str();
}

bool BiTargetCheckFaults::isOK(const SpaceTarget &T1,
                               const SpaceTarget &T2) const
{
  if (_faults == nullptr) return true;
  return !_faults->isSplitByFaultSP(T1.getCoordAsSP(), T2.getCoordAsSP());
}
