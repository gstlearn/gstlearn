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
#include "Space/SpacePoint.hpp"
#include "Space/ASpace.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Utilities.hpp"

#include <iostream>
#include <math.h>

SpacePoint::SpacePoint(const ASpace* space)
: ASpaceObject(space),
  _coord()
{
  // Initialize the point to the space origin
  // TODO : Not true whatever the space
  _coord = getOrigin();
}

SpacePoint::SpacePoint(const SpacePoint& r)
: ASpaceObject(r)
, _coord(r._coord)
{
}

SpacePoint::SpacePoint(const VectorDouble& coord,
                       const ASpace* space)
: ASpaceObject(space),
  _coord()
{
  if (coord.size() == 0)
  {
    // Use a valid default SpacePoint (origin ?)
    // TODO : Not true whatever the space
    _coord = getOrigin();
  }
  else if (coord.size() != getNDim())
  {
    messerr("Error: Wrong number of coordinates. Space origin used.");
    // TODO : Not true whatever the space
    _coord = getOrigin();
  }
  else
  {
    _coord = coord;
  }
}

SpacePoint& SpacePoint::operator=(const SpacePoint& r)
{
  if (this != &r)
  {
    ASpaceObject::operator=(r);
    _coord = r._coord;
  }
  return *this;
}

SpacePoint::~SpacePoint()
{
}

void SpacePoint::setCoord(double coord)
{
  VH::fill(_coord, coord, static_cast<int> (_coord.size()));
}

void SpacePoint::setCoord(const VectorDouble& coord)
{
  if (getNDim() != coord.size())
    std::cout << "Error: Wrong number of coordinates. Point not modified." << std::endl;
  else
    _coord = coord;
}

bool SpacePoint::isConsistent(const ASpace* space) const
{
  return (space->getNDim() == _coord.size());
}

void SpacePoint::move(const VectorDouble& vec)
{
  getSpace()->move(*this, vec);
}

double SpacePoint::getDistance(const SpacePoint& pt) const
{
  return ASpaceObject::getDistance(*this, pt);
}

VectorDouble SpacePoint::getIncrement(const SpacePoint& pt) const
{
  return ASpaceObject::getIncrement(*this, pt);
}

String SpacePoint::toString(const AStringFormat* /*strfmt*/) const
{
  return VH::toStringAsVD(_coord);
}

void SpacePoint::setFFFF()
{
  setCoord(TEST);
}

bool SpacePoint::isFFFF() const
{
  for (int idim = 0, ndim = getNDim(); idim < ndim; idim++)
    if (! FFFF(_coord[idim])) return false;
  return true;
}

/**
 * Initialize point coordinates from angles
 *
 * TODO : initialize coordinates from angles for more than 2D & valid only for space RN ?
 * To be kept ?
 */
#include <iostream>
void SpacePoint::setCoordFromAngle(const VectorDouble& angles)
{
  if (getNDim() == 1 || angles.size() == 0)
  {
    my_throw("Inconsistent angles vector");
  }
  else if (getNDim() == 2)
  {
    if (angles.size() > 1)
    {
      std::cout << "Warning: Extra angle values ignored" << std::endl;
    }
    _coord[0] = cos( GV_PI * angles[0] / 180);
    _coord[1] = sin( GV_PI * angles[0] / 180);
  }
  else
  {
    my_throw("Not yet implemented");
  }
}
