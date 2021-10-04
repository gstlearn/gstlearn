#include "Space/ASpaceObject.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpaceRN.hpp"
#include "Space/SpaceSN.hpp"
#include "Basic/Tensor.hpp"

#include "Basic/AException.hpp"

#include <iostream>

ASpace* ASpaceObject::_globalSpace = nullptr;
bool ASpaceObject::_fakeSpace = true;

ASpaceObject::ASpaceObject(const ASpace* space)
: _space(space) /// TODO : ASpace shared pointer
{
  if (nullptr == _space)
    _space = getGlobalSpace();
}

ASpaceObject::ASpaceObject(const ASpaceObject& r)
: _space(r._space) /// TODO : ASpace shared pointer
{

}
ASpaceObject& ASpaceObject::operator=(const ASpaceObject& r)
{
  if (this != &r)
  {
    /// TODO : ASpace shared pointer
    _space = r._space;
  }
  return *this;
}

ASpaceObject::~ASpaceObject()
{
}

const ASpace* ASpaceObject::createGlobalSpace(SpaceType type,
                                              unsigned int ndim,
                                              double param)
{
  if (nullptr != _globalSpace)
  {
    std::cout << "Cannot recreate global space context! Should never occur!" << std::endl;
    std::cout << "Creation aborted: previous Global Space is returned" << std::endl;
    return _globalSpace;
  }
  switch (type)
  {
    case SPACE_SN:
    {
      _globalSpace = new SpaceSN(ndim, param);
      break;
    }
    case SPACE_RN:
    {
      _globalSpace = new SpaceRN(ndim);
      break;
    }
    default:
    {
      my_throw("Unknown space type!");
    }
  }
  _fakeSpace = false;
  return _globalSpace;
}

const ASpace* ASpaceObject::getGlobalSpace()
{
  if (nullptr == _globalSpace)
  {
//    std::cout << "Creating default global space: SpaceRN 2D..." << std::endl;
//    std::cout << "Call ASpaceObject::createGlobalSpace to avoid this message!" << std::endl;
    createGlobalSpace(SPACE_RN, 2);
    _fakeSpace = true;
  }
  return _globalSpace;
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
  return (getSpace()->getNDim());
}

const SpacePoint& ASpaceObject::getOrigin() const
{
  return (getSpace()->getOrigin());
}

double ASpaceObject::getDistance(const SpacePoint& p1, const SpacePoint& p2) const
{
  return (getSpace()->getDistance(p1, p2));
}

VectorDouble ASpaceObject::getIncrement(const SpacePoint& p1, const SpacePoint& p2) const
{
  return (getSpace()->getIncrement(p1, p2));
}

bool ASpaceObject::isSpaceDimensionValid(int ndim)
{
  // If the (local) Space dimension is not defined, skip the test
  if (ndim <= 0) return false;

  // If (local) Space dimension is defined
  if (hasGlobalSpace())
  {
    // Compare the Local and the Global Space
    if ((unsigned int)(ndim) != ASpaceObject::getGlobalSpace()->getNDim())
    {
      std::cout << "Inconsistency in the Space Dimension" << std::endl;
      std::cout << "Local Space Dimension  =" << ndim << std::endl;
      std::cout << "Global Space Dimension =" << ASpaceObject::getGlobalSpace()->getNDim() << std::endl;
      return false;
    }
  }
  else
  {
    // Global Space is not defined yet, it is time to declare it
    createGlobalSpace(SPACE_RN, ndim);
    return true;
  }
  return true;
}
