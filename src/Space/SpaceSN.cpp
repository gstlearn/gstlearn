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

SpaceSN::SpaceSN(unsigned int ndim, double radius, bool addtime)
    : ASpace(ndim, addtime),
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

bool SpaceSN::_isEqual(const ASpace *space) const
{
  if (!ASpace::_isEqual(space)) return false;
  const SpaceSN* s = dynamic_cast<const SpaceSN*>(space);
  return s != nullptr && _radius == s->_radius;
}

void SpaceSN::_move(SpacePoint &p1, const VectorDouble &vec) const
{
  /// TODO : SpaceSN::_move
  for (unsigned int i = _iDimOffset; i < _nDim + _iDimOffset; i++)
  {
    p1.setCoord(i, p1.getCoord(i) + vec[i]);
  }
}

double SpaceSN::_getDistance(const SpacePoint &p1, 
                             const SpacePoint &p2) const
{
  return GH::geodeticAngularDistance(p1.getCoord(_iDimOffset),
                                     p1.getCoord(_iDimOffset + 1),
                                     p2.getCoord(_iDimOffset),
                                     p2.getCoord(_iDimOffset + 1),
                                     _radius);
}

double SpaceSN::_getDistance(const SpacePoint& p1,
                             const SpacePoint& p2,
                             const Tensor& tensor) const
{
  /// TODO : SpaceSN::_getDistance with tensor
  DECLARE_UNUSED(tensor);
  return GH::geodeticAngularDistance(p1.getCoord(_iDimOffset),
                                     p1.getCoord(_iDimOffset + 1),
                                     p2.getCoord(_iDimOffset),
                                     p2.getCoord(_iDimOffset + 1), 
                                     _radius);
}

double SpaceSN::_getDistance1D(double c1, double c2) const
{
  /// TODO : SpaceSN::_getDistance1D
  DECLARE_UNUSED(c1);
  DECLARE_UNUSED(c2);
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
  int j = 0;
  for (unsigned int i = _iDimOffset; i < _nDim + _iDimOffset; i++)
    ptemp[j++] = p2.getCoord(i) - p1.getCoord(i);
}