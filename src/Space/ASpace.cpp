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
#include "Space/ASpace.hpp"
#include "Basic/AStringable.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"

#include <iostream>
#include <sstream>

namespace gstlrn
{
ASpace::ASpace(size_t ndim)
  : AStringable()
  , _nDim(ndim)
  , _origin(VectorDouble(ndim, 0.))
  , _offset(0)
  , _work1(ndim)
  , _work2(ndim)
{
}

ASpace::ASpace(const ASpace& r)
  : AStringable(r)
  , _nDim(r._nDim)
  , _origin(r._origin)
  , _offset(r._offset)
  , _work1(r._nDim) // No need to copy the contents, just allocate
  , _work2(r._nDim)
{
}

ASpace& ASpace::operator=(const ASpace& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _nDim = r._nDim;
    _origin = r._origin;
    _offset = r._offset;
    _work1 = r._work1;
    _work2 = r._work2;
  }
  return *this;
}

ASpace::~ASpace() 
{
  //messerr("coucou");
}

/// Interface for AStringable
String ASpace::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  sstr << toString(strfmt, -1);
  if (strfmt != nullptr && strfmt->getLevel() == 0) sstr << std::endl;
  return sstr.str();
}

/// Update the origin of the space
void ASpace::setOrigin(const VectorDouble& origin) 
{
  if (origin.size() != getNDim())
  {
    std::cout << "Error: Inconsistent space origin. Origin not changed." << std::endl;
    return;
  }
  _origin = origin;
}

/// Get the number of dimensions
size_t ASpace::getNDim(Id ispace) const
{
  DECLARE_UNUSED(ispace)
  return _nDim;
}

/// Get the offset index for coordinates
size_t ASpace::getOffset(Id ispace) const
{
  DECLARE_UNUSED(ispace)
  return _offset;
}

/// Return the space origin coordinates
const VectorDouble& ASpace::getOrigin(Id ispace) const
{
  DECLARE_UNUSED(ispace)
  return _origin;
}

ASpaceSharedPtr ASpace::getDefaultSpaceIfNull(const ASpaceSharedPtr& space)
{
  if (space == nullptr) return getDefaultSpaceSh();
  return space;
}

/// Get the number of space components
size_t ASpace::getNComponents() const
{
  return 1;
}

/// Return the space component at index ispace
ASpaceSharedPtr ASpace::getComponent(Id ispace) const
{
  DECLARE_UNUSED(ispace)
  return ASpaceSharedPtr(this);
}

/// Dump a space in a string (given the space index)
String ASpace::toString(const AStringFormat* strfmt, Id ispace) const
{
  std::stringstream sstr;
  if (strfmt != nullptr && strfmt->getLevel() == 0)
  {
    sstr << getType().getKey() << "(" << getNDim() << ")";
  }
  else
  {
    if (ispace < 0)
    {
      sstr << "Space Type      = " << getType().getKey() << std::endl;
      sstr << "Space Dimension = " << getNDim() << std::endl;
    }
    else
    {
      sstr << "Space Type      [" << ispace << "] = " << getType().getKey() << std::endl;
      sstr << "Space Dimension [" << ispace << "] = " << getNDim() << std::endl;
    }
  }
  return sstr.str();
}

/// Return true if the given space is equal to me (same dimension and space definition)
bool ASpace::isEqual(const ASpace* space) const
{
  return (getNComponents() == space->getNComponents() &&
          getType() == space->getType() &&
          getNDim() == space->getNDim() &&
          getOrigin() == space->getOrigin() &&
          getOffset() == space->getOffset());
}

/// Return all the distances (one by space component) between two space points
VectorDouble ASpace::getDistances(const SpacePoint& p1,
                                  const SpacePoint& p2) const
{
  VectorDouble dis;
  if (p1.getNDim() != p2.getNDim())
  {
    std::cout << "Error: Inconsistent point dimension. Return empty distances"
              << std::endl;
    return dis;
  }
  dis.push_back(_getDistance(p1, p2));
  return dis;
}

/// Move the given space point by the given vector
void ASpace::move(SpacePoint& p1, const VectorDouble& vec) const
{
  if (vec.size() <= 0 ||
      vec.size() < getOffset() + getNDim() ||
      vec.size() != p1.getNDim())
  {
    std::cout << "Error: Inconsistent vector dimension. Point not moved."
              << std::endl;
    return;
  }
  _move(p1, vec);
}

VectorDouble ASpace::getUnitaryVector() const
{
  VectorDouble uni;
  uni.resize(getNDim(), 0.);
  uni[0] = 1;
  return uni;
}
/// Return the distance between two space points
double ASpace::getDistance(const SpacePoint &p1,
                           const SpacePoint &p2,
                           Id ispace) const
{
  if (p1.getNDim() != p2.getNDim())
  {
    std::cout << "Error: Inconsistent space dimension. Return TEST."
              << std::endl;
    return TEST;
  }
  return _getDistance(p1, p2, ispace);
}

/// Return the distance between two space points with the given tensor
double ASpace::getDistance(const SpacePoint& p1,
                           const SpacePoint& p2,
                           const Tensor& tensor,
                           Id ispace) const
{
  if (p1.getNDim() != p2.getNDim())
  /// TODO : test Tensor dimension
  {
    std::cout << "Error: Inconsistent space dimension. Return TEST."
              << std::endl;
    return TEST;
  }
  return _getDistance(p1, p2, tensor, ispace);
}

/// Return the distance in frequential domain between two space points with the given tensor
double ASpace::getFrequentialDistance(const SpacePoint& p1,
                                      const SpacePoint& p2,
                                      const Tensor& tensor,
                                      Id ispace) const
{
  if (p1.getNDim() != p2.getNDim())
  /// TODO : test tensor size
  {
    std::cout << "Error: Inconsistent point dimensions. Return TEST."
              << std::endl;
    return TEST;
  }
  return _getFrequentialDistance(p1, p2, tensor, ispace);
}

/// Return the increment vector between two space points
VectorDouble ASpace::getIncrement(const SpacePoint& p1,
                                  const SpacePoint& p2,
                                  Id ispace) const
{
  if (p1.getNDim() != p2.getNDim())
  /// TODO : test tensor size
  {
    std::cout << "Error: Inconsistent point dimensions. Return empty vector."
              << std::endl;
    return VectorDouble();
  }
  return _getIncrement(p1, p2, ispace);
}

void ASpace::getIncrementInPlace(const SpacePoint& p1,
                                 const SpacePoint& p2,
                                 VectorDouble& ptemp,
                                 Id ispace) const
{
  if (p1.getNDim() != p2.getNDim())
  /// TODO : test tensor size
  {
    std::cout << "Error: Inconsistent point dimensions. Return empty vector."
              << std::endl;
    ptemp.clear();
  }
  _getIncrementInPlace(p1, p2, ptemp, ispace);
}

VectorDouble ASpace::projCoord(const VectorDouble& coord, Id ispace) const
{
  if (ispace < 0 || ispace >= static_cast<Id>(getNComponents())) return coord;
  auto sp = getComponent(ispace);
  auto first       = coord.cbegin() + sp->getOffset();
  auto last        = first          + sp->getNDim();
  /// TODO : Memory copies !
  return VectorDouble(first, last);
}

/////////////////////////////////////////////////////////
}