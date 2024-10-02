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

SpacePoint::SpacePoint(const ASpace* space)
: ASpaceObject(space),
  _coord(new double[getNDim()]),
  _deleteCoord(true),
  _iech(-1)
{
  // Initialize the point to the space origin
  // TODO : Not true whatever the space
  for (int i = 0; i < (int)getNDim(); i++)
    _coord[i] = getOrigin().getVector()[i];
}

SpacePoint::SpacePoint(const SpacePoint& r)
: ASpaceObject(r)
,_coord(new double[getNDim()])
,_deleteCoord(true)
,_iech(-1)
{
  for (int i = 0; i < (int) getNDim(); i++)
    _coord[i] = r._coord[i];
}


double SpacePoint::getCoord(int idim) const
{
  return _coord[idim];
}

constvect SpacePoint::getCoords() const
{
  return constvect(_coord, getNDim());
}

SpacePoint::SpacePoint(vect coord, const ASpace* space,int iech)
: ASpaceObject(space)
{
  _coord = coord.data();
  _iech = iech;
}
SpacePoint::SpacePoint(const VectorDouble& coord,int iech,
                       const ASpace* space)
: ASpaceObject(space),
  _coord(new double[getNDim()]),
  _deleteCoord(true),
  _iech(iech)
{
  if (coord.size() == 0 || coord.size() != getNDim())
  {
    // Use a valid default SpacePoint (origin ?)
    // TODO : Not true whatever the spac
    messerr("Problem with the number of coordinates. \n");
    messerr("Point not created.\n");
    for (int i = 0; i < (int) getNDim(); i++)
      _coord[i] = getOrigin().getVector()[i];
  }
  else
  {
    for (int i = 0; i < (int)getNDim(); i++)
      _coord[i] = coord[i];
  }
}

SpacePoint& SpacePoint::operator=(const SpacePoint& r)
{
  _coord = new double[getNDim()];
  _deleteCoord = true;
  if (this != &r)
  {
    ASpaceObject::operator=(r);
    for (int i = 0; i < (int) getNDim(); i++)
      _coord[i] = r._coord[i];
  }
  return *this;
}

SpacePoint::~SpacePoint()
{
  if (_deleteCoord)
    delete [] _coord;
}

void SpacePoint::setCoord(double coord)
{
  std::fill(_coord,_coord + getNDim(),coord);
}

void SpacePoint::setCoords(const VectorDouble& coord)
{

  if ((int)getNDim() != (int)coord.size())
    std::cout << "Error: Wrong number of coordinates. Point not modified." << std::endl;
  else
  {
    for (int idim = 0; idim < (int)getNDim(); idim++)
      _coord[idim] = coord[idim];
  }
}

SpacePoint SpacePoint::projection(int ispace) const
{
  if (ispace < 0)
    return *this;
  else
  {
    int ndim = getNDim(ispace);
    SpacePoint p(vect(_coord,_coord+ndim),getSpace()->getComponent(ispace),_iech);
    return p;
  }
}
void SpacePoint::setCoords(const double* coord, int size)
{
 
  if ((int)getNDim() != size)
    std::cout << "Error: Wrong number of coordinates. Point not modified."
              << std::endl;
  else
    for (int idim = 0; idim < size; idim++)
     _coord[idim] = coord[idim];
}

bool SpacePoint::isConsistent(const ASpace* space) const
{
  DECLARE_UNUSED(space)
  //return (space->getNDim() == _coord.size());
  return true;
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

double SpacePoint::getDistance1D(const SpacePoint &pt, int idim) const
{
  return ASpaceObject::getDistance1D(*this, pt, idim);
}

VectorDouble SpacePoint::getIncrement(const SpacePoint& pt, int ispace) const
{
  return ASpaceObject::getIncrement(*this, pt, ispace);
}

String SpacePoint::toString(const AStringFormat* /*strfmt*/) const
{
  return VH::toStringAsSpan(constvect(_coord,getNDim()));
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

