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

#include "Geometry/BiPointCheckBench.hpp"
#include "Space/SpacePoint.hpp"

BiPointCheckBench::BiPointCheckBench(int idim_bench, double width)
    : ABiPointCheck(),
      _idimBench(idim_bench),
      _width(width)
{
}

BiPointCheckBench::BiPointCheckBench(const BiPointCheckBench &r)
    : ABiPointCheck(r),
      _idimBench(r._idimBench),
      _width(r._width)
{
}

BiPointCheckBench& BiPointCheckBench::operator=(const BiPointCheckBench &r)
{
  if (this != &r)
  {
    ABiPointCheck::operator=(r);
    _idimBench = r._idimBench;
    _width = r._width;
  }
  return *this;
}

BiPointCheckBench::~BiPointCheckBench()
{
}

BiPointCheckBench* BiPointCheckBench::create(int idim_bench, double width)
{
  return new BiPointCheckBench(idim_bench, width);
}

/**
 * Printout
 * @param
 * @return
 *
 * @remark The printout is not performed here as the Checker is only set in 'attach'
 */
String BiPointCheckBench::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  return sstr.str();
}

bool BiPointCheckBench::isOK(const SpacePoint &P1,
                             const SpacePoint &P2,
                             int iech1,
                             int iech2) const
{
  /* Discard sample located outside the bench */

  return (ABS(P1.getCoord(_idimBench) - P2.getCoord(_idimBench)) <= _width);
}
