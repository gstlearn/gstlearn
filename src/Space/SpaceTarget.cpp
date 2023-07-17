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
#include "Space/SpaceTarget.hpp"
#include "Space/ASpace.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Utilities.hpp"

#include <iostream>
#include <math.h>

SpaceTarget::SpaceTarget(const ASpace* space)
: SpacePoint(space),
  _extend(),
  _code(TEST),
  _date(TEST)
{
  _initialize();
}

SpaceTarget::SpaceTarget(const SpaceTarget &r)
    : SpacePoint(r),
      _extend(r._extend),
      _code(r._code),
      _date(r._date)
{
}

SpaceTarget& SpaceTarget::operator=(const SpaceTarget& r)
{
  if (this != &r)
  {
    SpacePoint::operator=(r);
    _extend = r._extend;
    _code = r._code;
    _date = r._date;
  }
  return *this;
}

SpaceTarget::~SpaceTarget()
{
}

SpaceTarget* SpaceTarget::create(const VectorDouble &center,
                                 const VectorDouble &extend,
                                 double code,
                                 double date,
                                 const ASpace *space)
{
  SpaceTarget* st = new SpaceTarget();
  st->setCoord(center);
  st->setExtend(extend);
  st->setCode(code);
  st->setDate(date);
  return st;
}

void SpaceTarget::_initialize()
{
  // Fill the extension with zeroes
  if (_extend.empty())
    VH::fill(_extend, getNDim(), 0.);
}

String SpaceTarget::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << "- Center    = " << VH::toStringAsVD(getCoord());
  if (! _extend.empty())
    sstr << "- Extension = " << VH::toStringAsVD(_extend);
  else
    sstr << "- Extension = (undefined)" << std::endl;
  if (! FFFF(_code))
    sstr << "- Code      = " << _code << std::endl;
  else
    sstr << "- Code      = (undefined)" << std::endl;
  if (! FFFF(_date))
    sstr << "- Date      = " << _date << std::endl;
  else
    sstr << "- Date      = (undefined)" << std::endl;

  return sstr.str();
}
