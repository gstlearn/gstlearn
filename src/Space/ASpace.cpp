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
#include "Space/SpacePoint.hpp"

#include <iostream>
#include <vector>

ASpace::ASpace(unsigned int ndim)
    : AStringable(),
      _nDim(ndim),
      _origin(VectorDouble(ndim, 0.)),
      _work1(ndim),
      _work2(ndim)
{
  if (ndim <= 0)
  {
    _nDim = 2;
  }
}

ASpace::ASpace(const ASpace& r)
    : AStringable(r),
      _nDim(r._nDim),
      _origin(r._origin),
      _work1(r._nDim), // No need to copy the contents, just allocated
      _work2(r._nDim)
{
}

ASpace& ASpace::operator=(const ASpace& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _nDim = r._nDim;
    _origin = r._origin;
    _work1 = r._work1;
    _work2 = r._work2;
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


void ASpace::_getIncrementInPlaceVect(const SpacePoint &p1,
                                      const std::vector<SpacePoint> &pv,
                                      VectorVectorDouble &res) const
{
	int np = (int)res.size();
	for(int i = 0; i<np;i++)
	{
		_getIncrementInPlace(p1,pv[i],res[i]);
	}
}

