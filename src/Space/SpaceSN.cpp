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
#include "Geometry/GeometryHelper.hpp"
#include <memory>

SpaceSN::SpaceSN(unsigned int ndim, double radius)
  : ASpace(ndim)
  , _radius(radius)
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

ASpaceSharedPtr SpaceSN::create(int ndim, double radius)
{
  return std::shared_ptr<SpaceSN>(new SpaceSN(ndim, radius));
}

String SpaceSN::toString(const AStringFormat* strfmt, int idx) const
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

bool SpaceSN::isEqual(const ASpace *space) const
{
  if (!ASpace::isEqual(space)) return false;
  const SpaceSN* s = dynamic_cast<const SpaceSN*>(space);
  return s != nullptr && _radius == s->_radius;
}

void SpaceSN::_move(SpacePoint &p1, const VectorDouble &vec) const
{
  /// TODO : SpaceSN::_move
  unsigned int offset = getOffset();
  unsigned int ndim   = getNDim();
  for (unsigned int i = offset; i < ndim + offset; i++)
  {
    p1.setCoord(i, p1.getCoord(i) + vec[i]);
  }
}

double SpaceSN::_getDistance(const SpacePoint& p1,
                             const SpacePoint& p2,
                             int ispace) const
{
  DECLARE_UNUSED(ispace)
  unsigned int offset = getOffset();
  return GH::geodeticAngularDistance(p1.getCoord(offset),
                                     p1.getCoord(offset + 1),
                                     p2.getCoord(offset),
                                     p2.getCoord(offset + 1),
                                     _radius);
}

double SpaceSN::_getDistance(const SpacePoint& p1,
                             const SpacePoint& p2,
                             const Tensor& tensor,
                             int ispace) const
{
  DECLARE_UNUSED(ispace)
  /// TODO : SpaceSN::_getDistance with tensor
  DECLARE_UNUSED(tensor);
  unsigned int offset = getOffset();
  return GH::geodeticAngularDistance(p1.getCoord(offset),
                                     p1.getCoord(offset + 1),
                                     p2.getCoord(offset),
                                     p2.getCoord(offset + 1), 
                                     _radius);
}

double SpaceSN::_getFrequentialDistance(const SpacePoint& p1,
                                        const SpacePoint& p2,
                                        const Tensor& tensor,
                                        int ispace) const
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
                                    int ispace) const
{
  DECLARE_UNUSED(ispace)
  _getIncrementInPlace(p1, p2, _work1);
  return _work1;
}

void SpaceSN::_getIncrementInPlace(const SpacePoint& p1,
                                   const SpacePoint& p2,
                                   VectorDouble& ptemp,
                                   int ispace) const
{
  DECLARE_UNUSED(ispace)
  /// TODO : SpaceSN::_getIncrementInPlace
  int j = 0;
  unsigned int offset = getOffset();
  unsigned int ndim   = getNDim();
  for (unsigned int i = offset; i < ndim + offset; i++)
    ptemp[j++] = p2.getCoord(i) - p1.getCoord(i);
}