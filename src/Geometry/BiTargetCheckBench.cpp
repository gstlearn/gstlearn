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
#include <Geometry/BiTargetCheckBench.hpp>
#include "geoslib_f.h"

#include "Space/SpacePoint.hpp"

BiTargetCheckBench::BiTargetCheckBench(int idim_bench, double width)
    : ABiTargetCheck(),
      _idimBench(idim_bench),
      _width(width)
{
}

BiTargetCheckBench::BiTargetCheckBench(const BiTargetCheckBench &r)
    : ABiTargetCheck(r),
      _idimBench(r._idimBench),
      _width(r._width)
{
}

BiTargetCheckBench& BiTargetCheckBench::operator=(const BiTargetCheckBench &r)
{
  if (this != &r)
  {
    ABiTargetCheck::operator=(r);
    _idimBench = r._idimBench;
    _width = r._width;
  }
  return *this;
}

BiTargetCheckBench::~BiTargetCheckBench()
{
}

BiTargetCheckBench* BiTargetCheckBench::create(int idim_bench, double width)
{
  return new BiTargetCheckBench(idim_bench, width);
}

String BiTargetCheckBench::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << "Bench width     = " << _width << std::endl;

  return sstr.str();
}

bool BiTargetCheckBench::isValid(const Db* dbin, const Db* dbout)
{
  _idimBench = dbin->getNDim() - 1;
  return true;
}

bool BiTargetCheckBench::isOK(const SpaceTarget &T1,
                             const SpaceTarget &T2) const
{
  /* Discard sample located outside the bench */

  return (ABS(T1.getCoord(_idimBench) - T2.getCoord(_idimBench)) <= _width);
}
