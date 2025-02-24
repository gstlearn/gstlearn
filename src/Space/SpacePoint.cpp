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
#include "Basic/AStringable.hpp"
#include "Space/ASpace.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Utilities.hpp"
#include "Space/ASpaceObject.hpp"
#include "geoslib_define.h"

#include <cmath>
#include <iostream>
#include <math.h>

SpacePoint::SpacePoint(const ASpaceSharedPtr& space)
  : ASpaceObject(space)
  , _coord()
  , _iech(-1)
  , _mode(1)
  , _isProjected(false)
{

  // Initialize the point to the space origin
  // TODO : Not true whatever the space
  _coord = getOrigin();
}

SpacePoint::SpacePoint(const SpacePoint& r)
  : ASpaceObject(r)
  , _coord(r._coord)
  , _iech(r._iech)
  , _mode(r._mode)
  , _isProjected(r._isProjected)
{
}

SpacePoint::SpacePoint(const VectorDouble& coord, int iech, const ASpaceSharedPtr& space)
  : ASpaceObject(space)
  , _coord(coord)
  , _iech(iech)
  , _mode(1)
  , _isProjected(false)
{
  if (coord.size() == 0 || coord.size() != getNDim())
  {
    // Use a valid default SpacePoint (origin ?)
    // TODO : Not true whatever the space
    messerr("Problem with the number of coordinates. \n");
    messerr("Point not created.\n");
    _coord = getOrigin();
  }
}

double SpacePoint::getCoord(int idim) const
{
  return _coord[idim];
}

constvect SpacePoint::getCoords() const
{
  return constvect(_coord.data(), getNDim());
}

SpacePoint& SpacePoint::operator=(const SpacePoint& r)
{
  if (this != &r)
  {
    ASpaceObject::operator=(r);
    _coord       = r._coord;
    _iech        = r._iech;
    _mode        = r._mode;
    _isProjected = r._isProjected;
  }
  return *this;
}

SpacePoint::~SpacePoint()
{
}

void SpacePoint::setCoord(double coord)
{
  _coord.fill(coord);
}

void SpacePoint::setCoords(const VectorDouble& coord)
{
  if ((int)getNDim() != (int)coord.size())
    std::cout << "Error: Wrong number of coordinates. Point not modified."
              << std::endl;
  else
    _coord = coord;
}

SpacePoint SpacePoint::spacePointOnSubspace(int ispace) const
{
  if (ispace < 0 || ispace >= (int)getNDim())
    return *this;

  /// TODO : Memory copies
  VectorDouble vec = getSpace()->projCoord(_coord, ispace);
  const auto sp = getSpace()->getComponent(ispace);
  SpacePoint p(vec, _iech, sp);
  p.setMode(_mode);
  return p;
}

void SpacePoint::setCoords(const double* coord, int size)
{
  if ((int)getNDim() != size)
    std::cout << "Error: Wrong number of coordinates. Point not modified." << std::endl;
  else
    for (int idim = 0; idim < size; idim++)
     _coord[idim] = coord[idim];
}

bool SpacePoint::isConsistent(const ASpace* space) const
{
  DECLARE_UNUSED(space)
  return (space->getNDim() == _coord.size());
}

void SpacePoint::move(const VectorDouble& vec)
{
  getSpace()->move(*this, vec);
}

double SpacePoint::getDistance(const SpacePoint& pt, int ispace) const
{
  return ASpaceObject::getDistance(*this, pt, ispace);
}

VectorDouble SpacePoint::getDistances(const SpacePoint& pt) const
{
  return ASpaceObject::getDistances(*this, pt);
}

VectorDouble SpacePoint::getIncrement(const SpacePoint& pt, int ispace) const
{
  return ASpaceObject::getIncrement(*this, pt, ispace);
}

String SpacePoint::toString(const AStringFormat* /*strfmt*/) const
{
  return VH::toStringAsSpan(constvect(_coord.data(),getNDim()));
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

double SpacePoint::getCosineToDirection(const SpacePoint &T2,
                                        const VectorDouble &codir) const
{
  double cosdir = 0.;
  double dn1 = 0.;
  double dn2 = 0.;
  VectorDouble delta = getIncrement(T2);
  for (int idim = 0; idim < (int)getNDim(); idim++)
  {
    cosdir += delta[idim] * codir[idim];
    dn1 += delta[idim] * delta[idim];
    dn2 += codir[idim] * codir[idim];
  }
  double prod = dn1 * dn2;
  if (prod <= 0.) return (1.);
  return (cosdir / sqrt(prod));
}

double SpacePoint::getOrthogonalDistance(const SpacePoint &P2,
                                         const VectorDouble &codir) const
{
  double dn1 = 0.;
  double dn2 = 0.;
  double v = 0.;
  double dproj = 0.;
  VectorDouble delta = getIncrement(P2);
  for (int idim = 0; idim < (int)getNDim(); idim++)
  {
    dproj += delta[idim] * codir[idim];
    dn1 += codir[idim] * codir[idim];
    dn2 += delta[idim] * delta[idim];
  }
  if (dn1 > 0.) v = sqrt(dn2 - dproj * dproj / dn1);
  return (v);
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

