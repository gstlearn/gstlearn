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
#include "Space/SpaceSN.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpacePoint.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorHelper.hpp"
#include "Geometry/GeometryHelper.hpp"

SpaceSN::SpaceSN(unsigned int ndim, double radius)
    : ASpace(ndim),
      _radius(radius)
{
  if (ndim != 2)
  my_throw("SN is only implemented for ndim=2 (sphere)");
}

SpaceSN::SpaceSN(const SpaceSN &r)
    : ASpace(r),
      _radius(r._radius)
{
}

SpaceSN& SpaceSN::operator=(const SpaceSN &r)
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

String SpaceSN::_toString(const AStringFormat* strfmt, int idx) const
{
  std::stringstream sstr;
  sstr << ASpace::_toString(strfmt, idx);
  if (strfmt->getLevel() == 1)
  {
    if (idx < 0)
    {
      sstr << "Sphere Radius   = " << _radius << std::endl;
    }
    else
    {
      sstr << "Sphere Radius   [" << idx << "] = " << _radius << std::endl;
    }
  }
  return sstr.str();
}

bool SpaceSN::_isEqual(const ASpace *space) const
{
  if (!ASpace::_isEqual(space)) return false;
  const SpaceSN *s = dynamic_cast<const SpaceSN*>(space);
  if (s == nullptr || _radius != s->_radius) return false;
  return true;
}

void SpaceSN::_move(SpacePoint &p1, const VectorDouble &vec) const
{
  /// TODO : SpaceSN::_move
  p1.setCoord(VH::add(p1.getCoord(), vec));
}

double SpaceSN::_getDistance(const SpacePoint &p1, const SpacePoint &p2) const
{
  return GH::geodeticAngularDistance(p1.getCoord(0),
                                     p1.getCoord(1),
                                     p2.getCoord(0),
                                     p2.getCoord(1),
                                     _radius);
}

double SpaceSN::_getDistance(const SpacePoint& p1,
                             const SpacePoint& p2,
                             const Tensor& tensor) const
{
  /// TODO : SpaceSN::_getDistance with tensor
  DECLARE_UNUSED(tensor);
  return GH::geodeticAngularDistance(p1.getCoord(0),
                                     p1.getCoord(1),
                                     p2.getCoord(0),
                                     p2.getCoord(1), 
                                     _radius);
}

double SpaceSN::_getDistance1D(const SpacePoint &p1,
                               const SpacePoint &p2,
                               int idim) const
{
  /// TODO : SpaceSN::_getDistance1D
  DECLARE_UNUSED(p1);
  DECLARE_UNUSED(p2);
  DECLARE_UNUSED(idim);
  return 0;
}

double SpaceSN::_getFrequentialDistance(const SpacePoint& p1,
                                        const SpacePoint& p2,
                                        const Tensor& tensor) const
{
  /// TODO : SpaceSN::_getFrequentialDistance
  DECLARE_UNUSED(p1);
  DECLARE_UNUSED(p2);
  DECLARE_UNUSED(tensor);
  return 0.;
}

VectorDouble SpaceSN::_getIncrement(const SpacePoint &p1,
                                    const SpacePoint &p2) const
{
  _getIncrementInPlace(p1, p2, _work1);
  return _work1;
}

/// Return the increment vector between two space points in a given vector
void SpaceSN::_getIncrementInPlace(const SpacePoint &p1,
                                   const SpacePoint &p2,
                                   VectorDouble &ptemp) const
{
  /// TODO : SpaceSN::_getIncrementInPlace
  ptemp = VH::subtract(p1.getCoord(), p2.getCoord());
}