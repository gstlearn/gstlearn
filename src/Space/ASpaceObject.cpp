#include "Space/ASpaceObject.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpaceRN.hpp"
#include "Space/SpaceSN.hpp"
#include "Basic/Tensor.hpp"

#include "Basic/AException.hpp"

#include <iostream>

ASpace* ASpaceObject::_defaultSpace = nullptr;

ASpaceObject::ASpaceObject(const ASpace* space)
: _space(nullptr)
{
  if (nullptr == space)
    // If the object is created without space, clone the global default space
    _space = cloneDefaultSpace();
  else
    // else duplicate the provided pointer
    _space = dynamic_cast<const ASpace*>(space->clone());
}

ASpaceObject::ASpaceObject(const ASpace& space)
: _space(nullptr)
{
  _space = dynamic_cast<const ASpace*>(space.clone());
}

ASpaceObject::ASpaceObject(const ASpaceObject& r)
: _space(nullptr)
{
  // Always duplicate!
  _space = dynamic_cast<const ASpace*>(r._space->clone());
}

ASpaceObject& ASpaceObject::operator=(const ASpaceObject& r)
{
  if (this != &r)
  {
    // Delete the previous space
    if (nullptr != _space)
      delete _space;
    // Clone the space of the object to be copied
    _space = dynamic_cast<const ASpace*>(r._space->clone());
  }
  return *this;
}

ASpaceObject::~ASpaceObject()
{
  delete _space;
}
/**
 * Factory for defining the unique default global space
 * (optional parameter can be used for sphere radius for example)
 *
 * @param type Space type (RN, SN, ...)
 * @param ndim Number of dimension
 * @param param Optional space parameter
 */
void ASpaceObject::defineDefaultSpace(SpaceType type,
                                      unsigned int ndim,
                                      double param)
{
  if (nullptr != _defaultSpace)
    delete _defaultSpace;

  switch (type)
  {
    case SPACE_SN:
    {
      _defaultSpace = new SpaceSN(ndim, param);
      break;
    }
    case SPACE_RN:
    {
      _defaultSpace = new SpaceRN(ndim);
      break;
    }
    default:
    {
      my_throw("Unknown space type!");
    }
  }
}

const ASpace* ASpaceObject::cloneDefaultSpace()
{
  if (nullptr == _defaultSpace)
    defineDefaultSpace(SPACE_RN, 2);

  return (dynamic_cast<const ASpace*>(_defaultSpace->clone()));
}

VectorDouble ASpaceObject::getUnitaryVector() const
{
  VectorDouble uni;
  uni.resize(getNDim(), 0.);
  uni[0] = 1;
  return uni;
}

unsigned int ASpaceObject::getNDim() const
{
  return (_space->getNDim());
}

const VectorDouble& ASpaceObject::getOrigin() const
{
  return (_space->getOrigin());
}

double ASpaceObject::getDistance(const SpacePoint& p1, const SpacePoint& p2) const
{
  return (_space->getDistance(p1, p2));
}

VectorDouble ASpaceObject::getIncrement(const SpacePoint& p1, const SpacePoint& p2) const
{
  return (_space->getIncrement(p1, p2));
}

