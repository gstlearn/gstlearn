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

String BiTargetCheckFaults::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << "- Separated by Faults" << std::endl;

  return sstr.str();
}

bool BiTargetCheckFaults::isOK(const SpaceTarget &T1,
                              const SpaceTarget &T2) const
{
  return _faults->isSplitByFaultSP(T1.getCoordAsSP(), T2.getCoordAsSP());
}