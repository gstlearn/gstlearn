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

/// Unique default global space
static ASpace* defaultSpace = nullptr;

ASpaceObject::ASpaceObject(const ASpace* space)
  : AStringable(),
    _space(nullptr)
{
  if (nullptr == space)
    // If the object is created without space, clone the global default space
    _space = cloneDefaultSpace();
  else
    // else duplicate the provided pointer
    _space = dynamic_cast<const ASpace*>(space->clone());
}

ASpaceObject::ASpaceObject(const ASpace& space)
  : AStringable(),
    _space(nullptr)
{
  _space = dynamic_cast<const ASpace*>(space.clone());
}

ASpaceObject::ASpaceObject(const ASpaceObject& r)
  : AStringable(r),
    _space(nullptr)
{
  // Always duplicate
  _space = dynamic_cast<const ASpace*>(r._space->clone());
}

ASpaceObject& ASpaceObject::operator=(const ASpaceObject& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
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
  // Always delete the space (always cloned before)
  delete _space;
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

double ASpaceObject::getDistance1D(const SpacePoint& p1,
                                   const SpacePoint& p2,
                                   int idim) const
{
  return (_space->getDistance1D(p1, p2, idim));
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
// TODO: this function should be removed as dangerous
void ASpaceObject::setNDim(int ndim)
{
  if (_space->getType() != ESpaceType::RN)
    my_throw("Object is not in Space RN");

  delete _space;
  _space = new SpaceRN(ndim);
}

/**
 * Factory for defining the unique default global space
 * (optional parameter can be used for sphere radius for example)
 *
 * @param type Space type (RN, SN, ...)
 * @param ndim Number of dimensions
 * @param param Optional space parameter (ex: radius of the sphere)
 * @param addtime Optional add time dimension (composit space)
 */
void defineDefaultSpace(ESpaceType type, unsigned int ndim, double param, bool addtime)
{
  if (nullptr != defaultSpace)
    delete defaultSpace;

  switch (type.getValue())
  {
    case ESpaceType::E_SN:
    {
      ndim = 2;
      if (param <= 0.) param = EARTH_RADIUS;
      defaultSpace = new SpaceSN(ndim, param, addtime);
      break;
    }
    case ESpaceType::E_RN:
    {
      defaultSpace = new SpaceRN(ndim, addtime);
      break;
    }
    default:
    {
      my_throw("Unknown space type!");
    }
  }
}

const ASpace* cloneDefaultSpace()
{
  if (nullptr == defaultSpace)
    defineDefaultSpace(ESpaceType::RN, 2);

  return (dynamic_cast<const ASpace*>(defaultSpace->clone()));
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
  if (nullptr == defaultSpace)
    defineDefaultSpace(ESpaceType::RN, 2);
  return defaultSpace;
}

bool isDefaultSpaceSphere()
{
  return (getDefaultSpaceType() == ESpaceType::SN);
}

