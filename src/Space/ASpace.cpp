#include "Space/ASpace.hpp"
#include "Space/SpacePoint.hpp"

#include <iostream>

ASpace::ASpace(unsigned int ndim)
    : AStringable(),
      _nDim(ndim),
      _origin(VectorDouble(ndim, 0.))
{
  if (ndim <= 0)
  {
    messerr(">>> Creating a Space with dimension 0 should be forbidden");
  }
}

ASpace::ASpace(const ASpace& r)
    : AStringable(r),
      _nDim(r._nDim),
      _origin(r._origin)
{
}

ASpace& ASpace::operator=(const ASpace& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _nDim = r._nDim;
    _origin = r._origin;
  }
  return *this;
}

ASpace::~ASpace()
{
}

String ASpace::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  sstr << "Space Type      = " << getType().getKey() << std::endl;
  sstr << "Space Dimension = " << getNDim() << std::endl;
  return sstr.str();
}

void ASpace::setOrigin(const VectorDouble& origin)
{
  /// TODO : not true whatever the space
  if (_origin.size() != getNDim())
    std::cout << "Error: Inconsistent space origin. Nothing changed." << std::endl;
  else
    _origin = origin;
}

bool ASpace::isEqual(const ASpace* space) const
{
  if (getType() == space->getType() && getOrigin() == space->getOrigin())
    return true;
  return false;
}
