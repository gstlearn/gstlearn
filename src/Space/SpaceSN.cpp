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
#include <Geometry/GeometryHelper.hpp>
#include "Space/SpaceSN.hpp"
#include "Space/SpacePoint.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorHelper.hpp"

SpaceSN::SpaceSN(unsigned int ndim, double radius)
: ASpace(ndim),
  _radius(radius)
{
  if (ndim != 2)
    my_throw("SN is only implemented for ndim=2 (sphere)");
}

SpaceSN::SpaceSN(const SpaceSN& r)
: ASpace(r),
  _radius(r._radius)
{
}

SpaceSN& SpaceSN::operator=(const SpaceSN& r)
{
  if (this != &r)
  {
    ASpace::operator=(r);
    _radius = r._radius;
  }
  return *this;
}

SpaceSN::~SpaceSN()
{
}

bool SpaceSN::isEqual(const ASpace* space) const
{
  if (!ASpace::isEqual(space)) return false;
  const SpaceSN* s = dynamic_cast<const SpaceSN*>(space);
  if (s == nullptr || _radius != s->_radius) return false;
  return true;
}

void SpaceSN::move(SpacePoint& p1,
                   const VectorDouble& vec) const
{
  /// TODO : SpaceSN::move
  p1.setCoord(VH::add(p1.getCoord(), vec));
}

double SpaceSN::getDistance(const SpacePoint& p1,
                            const SpacePoint& p2) const
{
  double long1 = p1.getCoord(0);
  double lat1  = p1.getCoord(1);
  double long2 = p2.getCoord(0);
  double lat2  = p2.getCoord(1);
  double dist = GH::geodeticAngularDistance(long1, lat1, long2, lat2, _radius);
  return dist;
}

double SpaceSN::getDistance(const SpacePoint& /*p1*/,
                            const SpacePoint& /*p2*/,
                            const Tensor& /*tensor*/) const
{
  /// TODO : SpaceSN::getDistance
  return 0.;
}

double SpaceSN::getFrequentialDistance(const SpacePoint& /*p1*/,
                              const SpacePoint& /*p2*/,
                              const Tensor& /*tensor*/) const
{
  return 0.;
}

VectorDouble SpaceSN::getIncrement(const SpacePoint& p1,
                                   const SpacePoint& p2) const
{
  /// TODO : SpaceSN::getIncrement
  return VH::subtract(p1.getCoord(), p2.getCoord());
}

String SpaceSN::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  sstr << "Space Type      = " << getType().getKey() << std::endl;
  sstr << "Space Dimension = " << getNDim() << std::endl;
  sstr << "Sphere Radius   = " << _radius << std::endl;
  return sstr.str();
}

