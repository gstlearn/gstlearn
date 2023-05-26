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

#include <iostream>
#include <math.h>

SpaceTarget::SpaceTarget(const ASpace* space)
: ASpaceObject(space),
  _center(),
  _extend()
{
  _initialize();
}

SpaceTarget::SpaceTarget(const SpaceTarget &r)
    : ASpaceObject(r),
      _center(r._center),
      _extend(r._extend)
{
}

SpaceTarget& SpaceTarget::operator=(const SpaceTarget& r)
{
  if (this != &r)
  {
    ASpaceObject::operator=(r);
    _center = r._center;
    _extend = r._extend;
  }
  return *this;
}

SpaceTarget::~SpaceTarget()
{
}

SpaceTarget* SpaceTarget::create(const VectorDouble &center,
                                 const VectorDouble &extend,
                                 const ASpace *space)
{
  SpaceTarget* st = new SpaceTarget();
  st->setCoord(center);
  st->setExtend(extend);
  return st;
}

void SpaceTarget::_initialize()
{
  // Initialize the point to the space origin
  if (_center.getCoord().empty()) _center = getOrigin();

  // Fill the extension with zeroes
  if (_extend.empty())
    VH::fill(_extend, getNDim(), 0.);
}

bool SpaceTarget::isConsistent(const ASpace* space) const
{
  return (space->getNDim() == _center.getNDim());
}

String SpaceTarget::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << "- Center    = " << VH::toString(_center.getCoord()) << std::endl;
  sstr << "- Extension = " << VH::toString(_extend) << std::endl;

  return sstr.str();
}
