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
#include "Basic/AException.hpp"
#include "Geometry/GeometryHelper.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpacePoint.hpp"
#include <memory>

namespace gstlrn
{
SpaceSN::SpaceSN(size_t ndim, double radius)
  : ASpace(ndim)
  , _radius(radius)
{
  if (ndim != 2)
    my_throw("SN is only implemented for ndim=2 (sphere)");
}

SpaceSN::SpaceSN(const SpaceSN& r)
  : ASpace(r)
  , _radius(r._radius)
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

ASpaceSharedPtr SpaceSN::create(Id ndim, double radius)
{
  return std::shared_ptr<SpaceSN>(new SpaceSN(ndim, radius));
}

String SpaceSN::toString(const AStringFormat* strfmt, Id idx) const
{
  std::stringstream sstr;
  sstr << ASpace::toString(strfmt, idx);
  if (strfmt == nullptr || strfmt->getLevel() == 1)
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

bool SpaceSN::isEqual(const ASpace* space) const
{
  if (!ASpace::isEqual(space)) return false;
  const auto* s = dynamic_cast<const SpaceSN*>(space);
  return s != nullptr && _radius == s->_radius;
}

void SpaceSN::_move(SpacePoint& p1, const VectorDouble& vec) const
{
  /// TODO : SpaceSN::_move
  auto offset = static_cast<Id>(getOffset());
  auto ndim   = static_cast<Id>(getNDim());
  for (Id i = offset; i < ndim + offset; i++)
  {
    p1.setCoord(i, p1.getCoord(i) + vec[i]);
  }
}

double SpaceSN::_getDistance(const SpacePoint& p1,
                             const SpacePoint& p2,
                             Id ispace) const
{
  DECLARE_UNUSED(ispace)
  auto offset = static_cast<Id>(getOffset());
  return GH::geodeticAngularDistance(p1.getCoord(offset),
                                     p1.getCoord(offset + 1),
                                     p2.getCoord(offset),
                                     p2.getCoord(offset + 1),
                                     _radius);
}

double SpaceSN::_getDistance(const SpacePoint& p1,
                             const SpacePoint& p2,
                             const Tensor& tensor,
                             Id ispace) const
{
  DECLARE_UNUSED(ispace)
  /// TODO : SpaceSN::_getDistance with tensor
  DECLARE_UNUSED(tensor);
  auto offset = static_cast<Id>(getOffset());
  return GH::geodeticAngularDistance(p1.getCoord(offset),
                                     p1.getCoord(offset + 1),
                                     p2.getCoord(offset),
                                     p2.getCoord(offset + 1),
                                     _radius);
}

double SpaceSN::_getFrequentialDistance(const SpacePoint& p1,
                                        const SpacePoint& p2,
                                        const Tensor& tensor,
                                        Id ispace) const
{
  DECLARE_UNUSED(ispace)
  /// TODO : SpaceSN::_getFrequentialDistance
  DECLARE_UNUSED(p1);
  DECLARE_UNUSED(p2);
  DECLARE_UNUSED(tensor);
  return 0.;
}

VectorDouble SpaceSN::_getIncrement(const SpacePoint& p1,
                                    const SpacePoint& p2,
                                    Id ispace) const
{
  DECLARE_UNUSED(ispace)
  _getIncrementInPlace(p1, p2, _work1);
  return _work1;
}

void SpaceSN::_getIncrementInPlace(const SpacePoint& p1,
                                   const SpacePoint& p2,
                                   VectorDouble& ptemp,
                                   Id ispace) const
{
  DECLARE_UNUSED(ispace)
  /// TODO : SpaceSN::_getIncrementInPlace
  Id j        = 0;
  auto offset = static_cast<Id>(getOffset());
  auto ndim   = static_cast<Id>(getNDim());
  for (Id i = offset; i < ndim + offset; i++)
    ptemp[j++] = p2.getCoord(i) - p1.getCoord(i);
}
} // namespace gstlrn
