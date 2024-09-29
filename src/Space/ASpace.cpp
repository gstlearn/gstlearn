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
      _nDim(ndim),
      _origin(VectorDouble(ndim, 0.)),
      _iDimOffset(0),
      _dimStart(),
      _comps(),
      _globalNDim(ndim),
      _globalOrigin(VectorDouble(ndim, 0.)),
      _work1(ndim),
      _work2(ndim)
{
  if (ndim <= 0)
  {
    _nDim = 2;
  }

  _dimStart.push_back(0);

  if(addTime)
  {
    // Time dimension is Euclidean
    // neglecting Eistein's relativity theory :-)
    SpaceRN* ts = new SpaceRN(1);
    addSpaceComponent(ts); // ts is cloned, so delete it
    _dimStart.push_back(2);
    delete ts;
  }
}

ASpace::ASpace(const ASpace& r)
    : AStringable(r),
      _nDim(r._nDim),
      _origin(r._origin),
      _iDimOffset(r._iDimOffset),
      _dimStart(r._dimStart),
      _comps(),
      _globalNDim(r._globalNDim),
      _globalOrigin(r._globalOrigin),
      _work1(r._globalNDim), // No need to copy the contents, just allocate
      _work2(r._globalNDim)
{
  for(auto* c : r._comps)
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
    _iDimOffset = r._iDimOffset;
    _dimStart = r._dimStart;
    _globalNDim = r._globalNDim;
    _globalOrigin = r._globalOrigin;
    _work1 = r._work1;
    _work2 = r._work2;
    for(auto* c : r._comps)
    {
      _comps.push_back(dynamic_cast<ASpace*>(c->clone()));
    }
  }
  return *this;
}

ASpace::~ASpace()
{
  for(auto* c : _comps)
  {
    delete c;
  }
}

void ASpace::addSpaceComponent(const ASpace* comp)
{
  ASpace* sp = dynamic_cast<ASpace*>(comp->clone());
  sp->_setDimOffset(_globalNDim);
  _dimStart.push_back(_globalNDim);
  _comps.push_back(sp);
  _globalNDim += sp->getNDim(0);
  const VectorDouble& o = sp->getOrigin(0);
  _globalOrigin.insert(_globalOrigin.end(), o.begin(), o.end());
}

unsigned int ASpace::getNComponents() const
{
  return (int)_comps.size() + 1;
}

String ASpace::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  int idx = getNComponents() > 1 ? 1 : -1;
  sstr << _toString(strfmt, idx);
  for(auto* c : _comps)
  {
    if (strfmt != nullptr && strfmt->getLevel() == 0)
      sstr << " + ";
    idx++;
    sstr << c->_toString(strfmt, idx);
  }
  if (strfmt != nullptr && strfmt->getLevel() == 0)
    sstr << std::endl;
  return sstr.str();
}

void ASpace::setOrigin(const VectorDouble& origin)
{
  /// TODO : not true whatever the space
  if (origin.size() != _globalNDim)
  {
    std::cout << "Error: Inconsistent space origin. Origin not changed." << std::endl;
    return;
  }
  _globalOrigin = origin;
  auto first = origin.cbegin();
  auto last = origin.cbegin() + _origin.size();
  _origin.assign(first, last);
  for(auto* c : _comps)
  {
    first = last;
    last = last + c->getNDim();
    c->setOrigin(VectorDouble(first, last));
  }
}

unsigned int ASpace::getNDim(int ispace) const
{
  if (ispace < 0 || ispace >= (int)getNComponents())
    return _globalNDim;
  if (ispace == 0)
    return _nDim;
  return _comps[ispace - 1]->getNDim(0);
}

const VectorDouble& ASpace::getOrigin(int ispace) const
{
  if (ispace < 0 || ispace >= (int)getNComponents())
    return _globalOrigin;
  if (ispace == 0)
    return _origin;
  return _comps[ispace-1]->getOrigin(0);
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
  if (vec.size() != _globalNDim)
  {
    std::cout << "Error: Inconsistent vector dimension. Point not moved." << std::endl;
    return;
  }
  _move(p1, vec);
  for (auto* sp : _comps)
  {
    sp->_move(p1, vec);
  }
}

/// Return the distance between two space points
double ASpace::getDistance(const SpacePoint &p1,
                           const SpacePoint &p2,
                           int ispace) const
{
  if (ispace <= 0 || ispace >= (int)getNComponents())
    return _getDistance(p1, p2);
  return _comps[ispace-1]->_getDistance(p1, p2);
}

/// Return the distance between two space points with the given tensor
double ASpace::getDistance(const SpacePoint& p1,
                           const SpacePoint& p2,
                           const Tensor& tensor,
                           int ispace) const
{
  if (ispace <= 0 || ispace >= (int)getNComponents())
    return _getDistance(p1, p2, tensor);
  return _comps[ispace-1]->_getDistance(p1, p2, tensor);
}

/// Return all the distances (one by space component) between two space points
VectorDouble ASpace::getDistances(const SpacePoint &p1,
                                  const SpacePoint &p2) const
{
  VectorDouble dis;
  dis.push_back(_getDistance(p1, p2));
  for (auto* sp : _comps)
  {
    dis.push_back(sp->_getDistance(p1, p2));
  }
  return dis;
}

/// Return the distance along one direction between two space points
double ASpace::getDistance1D(const SpacePoint &p1,
                             const SpacePoint &p2,
                             unsigned int idim) const
{
  if (idim >= _globalNDim)
  {
    std::cout << "Error: Inconsistent space dimension. Return TEST."
              << std::endl;
    return TEST;
  }

  if (idim < _nDim)
    return _getDistance1D(p1.getCoord(idim), p2.getCoord(idim));

  unsigned int jdim = idim - _nDim;
  for (auto* sp : _comps)
  {
    unsigned int ndim = sp->getNDim(0);
    if (jdim < ndim)
      return sp->_getDistance1D(p1.getCoord(idim), p2.getCoord(idim));
    jdim -= ndim;
  }

  return TEST;
}

/// Return the distance in frequential domain between two space points with the given tensor
double ASpace::getFrequentialDistance(const SpacePoint& p1,
                                      const SpacePoint& p2,
                                      const Tensor& tensor,
                                      int ispace) const
{
  if (ispace <= 0 || ispace >= (int)getNComponents())
    return _getFrequentialDistance(p1, p2, tensor);
  return _comps[ispace - 1]->_getFrequentialDistance(p1, p2, tensor);
}

/// Return the increment vector between two space points
VectorDouble ASpace::getIncrement(const SpacePoint& p1,
                                  const SpacePoint& p2,
                                  int ispace) const
{
  if (ispace <= 0 || ispace >= (int)getNComponents())
    return _getIncrement(p1, p2);
  return _comps[ispace - 1]->_getIncrement(p1, p2);
}

/////////////////////////////////////////////////////////

String ASpace::_toString(const AStringFormat* strfmt, int idx) const
{
  std::stringstream sstr;
  if (strfmt != nullptr && strfmt->getLevel() == 0)
  {
    sstr << getType().getKey() << "(" << getNDim(0) << ")";
  }
  else
  {
    if (idx < 0)
    {
      sstr << "Space Type      = " << getType().getKey() << std::endl;
      sstr << "Space Dimension = " << getNDim(0) << std::endl;
    }
    else
    {
      sstr << "Space Type      [" << idx << "] = " << getType().getKey() << std::endl;
      sstr << "Space Dimension [" << idx << "] = " << getNDim(0) << std::endl;
    }
  }
  return sstr.str();
}

bool ASpace::_isEqual(const ASpace* space) const
{
  return (getNDim()      == space->getNDim()   &&
          getType()      == space->getType()   &&
          getOrigin()    == space->getOrigin() &&
          getDimOffset() == space->getDimOffset());
}
