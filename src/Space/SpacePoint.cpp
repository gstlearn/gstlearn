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
#include "geoslib_define.h"

#include <cmath>
#include <iostream>
#include <math.h>

SpacePoint::SpacePoint(const ASpace* space)
: ASpaceObject(space),
  _coord(new double[getNDim()]),
  _deleteCoord(true)
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
{
  for (int i = 0; i < (int) getNDim(); i++)
    _coord[i] = r._coord[i];
}


double SpacePoint::getCoord(int idim) const
{
  int offset = getSpace()->getCurrentOffset();
  return _coord[idim + offset];
}

constvect SpacePoint::getCoords() const
{
  int offset = getSpace()->getCurrentOffset();
  int ndim = getSpace()->getCurrentNDim();
  return constvect(_coord+ offset, ndim);
}

SpacePoint::SpacePoint(const VectorDouble& coord,
                       const ASpace* space)
: ASpaceObject(space),
  _coord(new double[getNDim()]),
  _deleteCoord(true)
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
  int offset = getSpace()->getCurrentOffset();
  int ndim = getSpace()->getCurrentNDim();
  std::fill(_coord+offset,_coord+ ndim,coord);
}

void SpacePoint::setCoords(const VectorDouble& coord)
{
  int offset = getSpace()->getCurrentOffset();
  int ndim = getSpace()->getCurrentNDim();
  if (ndim != (int)coord.size())
    std::cout << "Error: Wrong number of coordinates. Point not modified." << std::endl;
  else
  {
    for (int idim = 0; idim < ndim; idim++)
      _coord[idim+offset] = coord[idim];
  }
}

void SpacePoint::setCoords(const double* coord, int size)
{
  int offset = getSpace()->getCurrentOffset();
  int ndim = getSpace()->getCurrentNDim();
  if (ndim != size)
    std::cout << "Error: Wrong number of coordinates. Point not modified."
              << std::endl;
  else
    for (int idim = 0; idim < size; idim++)
     _coord[idim+offset] = coord[idim];
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
  int ndim = getSpace()->getCurrentNDim();
  int offset = getSpace()->getCurrentOffset();
  if (getSpace()->getNComponents() > 1)
    message("Rank of the current space %d\n", getSpace()->getSpaceRankView());
  return VH::toStringAsSpan(constvect(_coord + offset,ndim));
}

void SpacePoint::setFFFF()
{
  setCoord(TEST);
}

bool SpacePoint::isFFFF() const
{

  for (int idim = 0, ndim = getNDim(0); idim < ndim; idim++)
    if (! FFFF(_coord[idim])) return false;
  return true;
}

double SpacePoint::getCosineToDirection(const SpacePoint &T2,
                                        const VectorDouble &codir) const
{
  int viewRank = getSpace()->getSpaceRankView();
  int ndim = getSpace()->getCurrentNDim();
  double cosdir = 0.;
  double dn1 = 0.;
  double dn2 = 0.;
  VectorDouble delta = getIncrement(T2,viewRank);
  for (int idim = 0; idim < ndim; idim++)
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
  int viewRank = getSpace()->getSpaceRankView();
  int ndim = getSpace()->getCurrentNDim();
  double dn1 = 0.;
  double dn2 = 0.;
  double v = 0.;
  double dproj = 0.;
  VectorDouble delta = getIncrement(P2,viewRank);
  for (int idim = 0; idim < ndim; idim++)
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
  int offset = getSpace()->getCurrentOffset();
  int ndim = getSpace()->getCurrentNDim();
  if (ndim == 1 || angles.size() == 0)
  {
    my_throw("Inconsistent angles vector");
  }
  else if (ndim == 2)
  {
    if (angles.size() > 1)
    {
      std::cout << "Warning: Extra angle values ignored" << std::endl;
    }
    _coord[offset] = cos( GV_PI * angles[0] / 180);
    _coord[offset + 1] = sin( GV_PI * angles[0] / 180);
  }
  else
  {
    my_throw("Not yet implemented");
  }
}

