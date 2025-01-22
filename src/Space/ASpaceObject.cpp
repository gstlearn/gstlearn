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
#include "geoslib_define.h"

#include "Space/ASpaceObject.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpaceRN.hpp"
#include "Space/SpaceSN.hpp"
#include "Basic/AException.hpp"
#include <memory>

/// Unique default global space
static std::shared_ptr<const ASpace> defaultSpace = nullptr;

ASpaceObject::ASpaceObject(const ASpaceSharedPtr& space)
  : AStringable(),
    _space(ASpace::getDefaultSpaceIfNull(space))
{

}


ASpaceObject::ASpaceObject(const ASpaceObject& r)
  : AStringable(r),
    _space(r._space)
{
}

bool ASpaceObject::isConsistent(std::shared_ptr<const ASpace> space) const
{
  return (isConsistent(space.get()));
}
ASpaceObject& ASpaceObject::operator=(const ASpaceObject& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _space = r._space;
  }
  return *this;
}

ASpaceObject::~ASpaceObject()
{

}

/// AStringable interface
String ASpaceObject::toString(const AStringFormat* /*strfmt*/) const
{
  messerr("ASpaceObject: 'toString' not yet implemented");
  return "";
}

VectorDouble ASpaceObject::getUnitaryVector() const
{
  VectorDouble uni;
  uni.resize(getNDim(), 0.);
  uni[0] = 1;
  return uni;
}

unsigned int ASpaceObject::getNDim(int ispace) const
{
  return (_space->getNDim(ispace));
}

const VectorDouble& ASpaceObject::getOrigin(int ispace) const
{
  if (_space == nullptr)
    return _dummy;
  return (_space->getOrigin(ispace));
}

double ASpaceObject::getDistance(const SpacePoint& p1,
                                 const SpacePoint& p2,
                                 int ispace) const
{
  return (_space->getDistance(p1, p2, ispace));
}

VectorDouble ASpaceObject::getDistances(const SpacePoint& p1,
                                        const SpacePoint& p2) const
{
  return (_space->getDistances(p1, p2));
}

VectorDouble ASpaceObject::getIncrement(const SpacePoint& p1,
                                        const SpacePoint& p2,
                                        int ispace) const
{
  return (_space->getIncrement(p1, p2, ispace));
}

/**
 * Modify the Space dimension of an already created item
 * (To be used only during creation ... in particular when reading NF)
 * @param ndim
 */
void ASpaceObject::setNDim(int ndim)
{
  if (_space->getType() != ESpaceType::RN)
    my_throw("Object is not in Space RN");

  _space = std::shared_ptr<const ASpace>(new SpaceRN(ndim));
}

/**
 * Factory for defining the unique default global space
 * (optional parameter can be used for sphere radius for example)
 *
 * @param type Space type (RN, SN, ...)
 * @param ndim Number of dimensions
 * @param param Optional space parameter (ex: radius of the sphere)
 */
void defineDefaultSpace(const ESpaceType& type, unsigned int ndim, double param)
{

  switch (type.getValue())
  {
    case ESpaceType::E_SN:
    {
      ndim = 2;
      if (param <= 0.) param = EARTH_RADIUS;
      defaultSpace = std::shared_ptr<const ASpace>(new SpaceSN(ndim, param));
      break;
    }
    case ESpaceType::E_RN:
    {
      defaultSpace = std::shared_ptr<const ASpace>(new SpaceRN(ndim));
      break;
    }
    default:
    {
      my_throw("Unknown space type!");
    }
  }
}

/**
 * @brief Defining the default space from another one
 * 
 * @param space 
 */
void setDefaultSpace(const std::shared_ptr<const ASpace> &space)
{
  defaultSpace = space;
}


ESpaceType getDefaultSpaceType()
{
  if (nullptr == defaultSpace)
    defineDefaultSpace(ESpaceType::RN, 2);
  return defaultSpace->getType();
}

int getDefaultSpaceDimension()
{
  if (nullptr == defaultSpace)
    defineDefaultSpace(ESpaceType::RN, 2);
  return defaultSpace->getNDim();
}

const ASpace* getDefaultSpace()
{
  return getDefaultSpaceSh().get();
}


std::shared_ptr<const ASpace> getDefaultSpaceSh()
{
  if (nullptr == defaultSpace)
    defineDefaultSpace(ESpaceType::RN, 2);
  return defaultSpace;
}

bool isDefaultSpaceSphere()
{
  return (getDefaultSpaceType() == ESpaceType::SN);
}

