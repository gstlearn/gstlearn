#include "Space/ASpace.hpp"
#include "Space/SpacePoint.hpp"

#include "geoslib_f.h"

#include <iostream>

ASpace::ASpace(unsigned int ndim)
: _nDim(ndim),
  _origin(this)
{
  _origin.setCoord(VectorDouble(ndim, 0.));
}

ASpace::ASpace(const ASpace& r)
: _nDim(r._nDim),
  _origin(r._origin)
{

}
ASpace& ASpace::operator=(const ASpace& r)
{
  if (this != &r)
  {
    _nDim = r._nDim;
    _origin = r._origin;
  }
  return *this;
}

ASpace::~ASpace()
{
}

std::string ASpace::toString(int level) const
{
  std::stringstream sstr;
  sstr << "Space Type      = " << getType() << std::endl;
  sstr << "Space Dimension = " << getNDim() << std::endl;
  return sstr.str();
}

void ASpace::setOrigin(const SpacePoint& origin)
{
  if (!origin.isConsistent(this))
    std::cout << "Error: Inconsistent space origin. Nothing changed." << std::endl;
  else
    _origin = origin;
}

bool ASpace::isEqual(const ASpace* space) const
{
  if (space == this) return true;
  if (getOrigin() == space->getOrigin()) return true;
  return false;
}
