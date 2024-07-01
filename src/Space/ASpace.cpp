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
#include "Space/SpaceRN.hpp"

#include <iostream>
#include <sstream>
#include <vector>

ASpace::ASpace(unsigned int ndim, bool addTime)
    : AStringable(),
      _comps(),
      _nDim(ndim),
      _origin(VectorDouble(ndim, 0.)),
      _globalOrigin(VectorDouble(ndim, 0.)),
      _work1(ndim),
      _work2(ndim)
{
  if (ndim <= 0)
  {
    _nDim = 2;
  }

  if(addTime)
  {
    // Time dimension is Euclidean
    // neglecting Eistein's relativity theory :-)
    _comps.push_back(new SpaceRN(1));
    _globalOrigin.push_back(0.);
  }
}

ASpace::ASpace(const ASpace& r)
    : AStringable(r),
      _comps(),
      _nDim(r._nDim),
      _origin(r._origin),
      _globalOrigin(r._globalOrigin),
      _work1(r._nDim), // No need to copy the contents, just allocate
      _work2(r._nDim)
{
  for(auto c : r._comps)
  {
    _comps.push_back(dynamic_cast<ASpace*>(c->clone()));
  }
}

ASpace& ASpace::operator=(const ASpace& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _nDim = r._nDim;
    _origin = r._origin;
    _globalOrigin = r._globalOrigin;
    _work1 = r._work1;
    _work2 = r._work2;
    for(auto c : r._comps)
    {
      _comps.push_back(dynamic_cast<ASpace*>(c->clone()));
    }
  }
  return *this;
}

ASpace::~ASpace()
{
  for(auto c : _comps)
  {
    delete c;
  }
}

void ASpace::addSpaceComponent(const ASpace* comp)
{
  ASpace* sp = dynamic_cast<ASpace*>(comp->clone());
  _comps.push_back(sp);
  const VectorDouble& o = sp->getOrigin();
  _comps.push_back(sp);
  _globalOrigin.insert(_globalOrigin.end(), o.begin(), o.end());
}

String ASpace::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  int idx = getNComponents() > 1 ? 1 : -1;
  sstr << _toString(strfmt, idx);
  for(auto c : _comps)
  {
    if (strfmt->getLevel() == 0)
      sstr << " + ";
    sstr << c->_toString(strfmt, idx++);
  }
  return sstr.str();
}

void ASpace::setOrigin(const VectorDouble& origin)
{
  /// TODO : not true whatever the space
  if (origin.size() != getNDim())
  {
    std::cout << "Error: Inconsistent space origin. Nothing changed." << std::endl;
    return;
  }
  _globalOrigin = origin;
  auto first = origin.cbegin();
  auto last = origin.cbegin() + _origin.size();
  _origin.assign(first, last);
  for(auto c : _comps)
  {
    first = last;
    last = last + c->getNDim();
    c->setOrigin(VectorDouble(first, last));
  }
}

unsigned int ASpace::getNDim() const
{
  unsigned int ndim(_nDim);
  for(auto c : _comps)
  {
    ndim += c->getNDim();
  }
  return ndim;
}

unsigned int ASpace::getNComponents() const
{
  return _comps.size() + 1;
}

const VectorDouble& ASpace::getOrigin() const
{
  return _globalOrigin;
}

bool ASpace::isEqual(const ASpace* space) const
{
  bool equal = getNComponents() == space->getNComponents() && _isEqual(space);
  for(unsigned int i = 0; equal && i < _comps.size(); i++)
  {
    equal &= _comps[i]->isEqual(space->_comps[i]);
  }
  return equal;
}

void ASpace::move(SpacePoint& p1, const VectorDouble& vec) const
{
  _move(p1, vec);
}

/// Return the distance between two space points
double ASpace::getDistance(const SpacePoint &p1,
                           const SpacePoint &p2) const
{
  return _getDistance(p1, p2);
}

/// Return the distance between two space points with the given tensor
double ASpace::getDistance(const SpacePoint& p1,
                           const SpacePoint& p2,
                           const Tensor& tensor) const
{
  return _getDistance(p1, p2, tensor);
}

/// Return the distance along one direction between two space points
double ASpace::getDistance1D(const SpacePoint &p1,
                             const SpacePoint &p2,
                             int idim) const
{
  return _getDistance1D(p1, p2, idim);
}

/// Return the distance in frequential domain between two space points with the given tensor
double ASpace::getFrequentialDistance(const SpacePoint& p1,
                                      const SpacePoint& p2,
                                      const Tensor& tensor) const
{
  return _getFrequentialDistance(p1, p2, tensor);
}

/// Return the increment vector between two space points
VectorDouble ASpace::getIncrement(const SpacePoint& p1,
                                  const SpacePoint& p2) const
{
  return _getIncrement(p1, p2);
}

/////////////////////////////////////////////////////////

String ASpace::_toString(const AStringFormat* strfmt, int idx) const
{
  std::stringstream sstr;
  if (strfmt->getLevel() == 0)
  {
    sstr << getType().getKey() << "(" << getNDim() << ")";
  }
  else
  {
    if (idx < 0)
    {
      sstr << "Space Type      = " << getType().getKey() << std::endl;
      sstr << "Space Dimension = " << getNDim() << std::endl;
    }
    else
    {
      sstr << "Space Type      [" << idx << "] = " << getType().getKey() << std::endl;
      sstr << "Space Dimension [" << idx << "] = " << getNDim() << std::endl;
    }
  }
  return sstr.str();
}

bool ASpace::_isEqual(const ASpace* space) const
{
  return (getNDim()   == space->getNDim() &&
          getType()   == space->getType() &&
          getOrigin() == space->getOrigin());
}
