#include "Space/SpacePoint.hpp"
#include "Space/ASpace.hpp"

#include "geoslib_f.h"

#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"

#include <iostream>

SpacePoint::SpacePoint(const ASpace* space)
: ASpaceObject(space),
  _coord()
{
  // Initialize the point to the space origin
  _coord = getOrigin().getCoord();
}

SpacePoint::SpacePoint(const VectorDouble& coord,
                       const ASpace* space)
: ASpaceObject(space),
  _coord()
{
  if (coord.size() != getNDim())
  {
    std::cout << "Error: Wrong number of coordinates. Space origin used." << std::endl;
    _coord = getOrigin().getCoord();
  }
  else
  {
    _coord = coord;
  }
}

SpacePoint::~SpacePoint()
{
}

void SpacePoint::setCoord(double coord)
{
  ut_vector_fill(_coord, coord, static_cast<int> (_coord.size()));
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

String SpacePoint::toString(int level) const
{
  return ut_vector_string(_coord);
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
