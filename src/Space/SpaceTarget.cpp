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
#include "Space/SpaceTarget.hpp"
#include "Space/ASpace.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Utilities.hpp"

#include <iostream>
#include <math.h>


SpaceTarget::SpaceTarget(const std::shared_ptr<const ASpace>& space,
                         bool checkExtend,
                         bool checkCode,
                         bool checkDate)
  : SpacePoint(space)
  , _checkExtend(checkExtend)
  , _checkCode(checkCode)
  , _checkDate(checkDate)
  , _extend()
  , _code(TEST)
  , _date(TEST)
{
  _initialize();
}

SpaceTarget::SpaceTarget(const SpaceTarget& r)
  : SpacePoint(r)
  , _checkExtend(r._checkExtend)
  , _checkCode(r._checkCode)
  , _checkDate(r._checkDate)
  , _extend(r._extend)
  , _code(r._code)
  , _date(r._date)
{
}

SpaceTarget& SpaceTarget::operator=(const SpaceTarget& r)
{
  if (this != &r)
  {
    SpacePoint::operator=(r);
    _checkExtend = r._checkExtend;
    _checkCode   = r._checkCode;
    _checkDate = r._checkDate;
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
                                 const std::shared_ptr<const ASpace>& space)
{
  SpaceTarget* st = new SpaceTarget(space);
  st->setCoords(center);
  st->setExtend(extend);
  st->setCode(code);
  st->setDate(date);
  return st;
}

void SpaceTarget::_initialize()
{
  // Fill the extension with zeroes
  if (_extend.empty())
    VH::fill(_extend, 0., getNDim());
}

String SpaceTarget::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << "- Center    = " << VH::toStringAsSpan(getCoords());
  if (_checkExtend)
  {
    if (!_extend.empty())
      sstr << "- Extension = " << VH::toStringAsVD(_extend);
    else
      sstr << "- Extension = (undefined)" << std::endl;
  }
  if (_checkCode)
  {
    if (!FFFF(_code))
      sstr << "- Code      = " << _code << std::endl;
    else
      sstr << "- Code      = (undefined)" << std::endl;
  }
  if (_checkDate)
  {
    if (!FFFF(_date))
      sstr << "- Date      = " << _date << std::endl;
    else
      sstr << "- Date      = (undefined)" << std::endl;
  }
  return sstr.str();
}
