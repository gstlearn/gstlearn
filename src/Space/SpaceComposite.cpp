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
#include "Space/SpaceComposite.hpp"

#include "Space/SpacePoint.hpp"

#include <iostream>
#include <sstream>
#include <vector>

SpaceComposite::SpaceComposite()
  : ASpace(0)
  , _comps()
{
}

SpaceComposite::SpaceComposite(const std::vector<const ASpace*>& vectspace)
  : ASpace(0)
  , _comps()
{
  for (const auto* sp : vectspace)
  {
    addSpaceComponent(sp);
  }
}

SpaceComposite::SpaceComposite(const SpaceComposite& r)
  : ASpace(r)
  , _comps()
{
  for(auto* c : r._comps)
  {
    _comps.push_back(dynamic_cast<ASpace*>(c->clone()));
  }
}

SpaceComposite& SpaceComposite::operator=(const SpaceComposite& r)
{
  if (this != &r)
  {
    _destroyComponents();
    ASpace::operator=(r);
    for(auto* c : r._comps)
    {
      _comps.push_back(dynamic_cast<ASpace*>(c->clone()));
    }
  }
  return *this;
}

SpaceComposite::~SpaceComposite()
{
  _destroyComponents();
}

void SpaceComposite::setOrigin(const VectorDouble& origin)
{
  if (origin.size() != ASpace::getNDim())
  {
    std::cout << "Error: Inconsistent space origin. Origin not changed." << std::endl;
    return;
  }
  _origin = origin;
  auto first = origin.cbegin();
  auto last = origin.cbegin();
  for(auto* c : _comps)
  {
    first = last;
    last = last + c->getNDim();
    c->setOrigin(VectorDouble(first, last));
  }
}

unsigned int SpaceComposite::getNDim(int ispace) const
{
  if (ispace < 0 || ispace >= (int)getNComponents())
    return ASpace::getNDim();
  return _comps[ispace]->getNDim();
}

unsigned int SpaceComposite::getOffset(int ispace) const
{
  if (ispace < 0 || ispace >= (int)getNComponents())
    return ASpace::getOffset();
  return _comps[ispace]->getOffset();
}

const VectorDouble& SpaceComposite::getOrigin(int ispace) const
{
  if (ispace < 0 || ispace >= (int)getNComponents())
    return ASpace::getOrigin();
  return _comps[ispace]->getOrigin();
}

unsigned int SpaceComposite::getNComponents() const
{
  return (int)_comps.size();
}

const ASpace* SpaceComposite::getComponent(int ispace) const
{
  if (ispace < 0 || ispace >= (int)getNComponents())
    return ASpace::getComponent(); // Return this if wrong ispace
  return _comps[ispace];
}

String SpaceComposite::toString(const AStringFormat* strfmt, int ispace) const
{
  DECLARE_UNUSED(ispace)
  std::stringstream sstr;
  sstr << ASpace::toString(strfmt, -1);
  if (strfmt != nullptr && strfmt->getLevel() == 0) sstr << ": ";
  unsigned int nc = getNComponents();
  for (unsigned int idx = 0; idx < nc; idx++)
  {
    const auto* c = getComponent(idx);
    sstr << c->toString(strfmt, idx);
    if (idx < nc - 1 && strfmt != nullptr && strfmt->getLevel() == 0) sstr << " + ";
  }
  if (strfmt != nullptr && strfmt->getLevel() == 0) sstr << std::endl;
  return sstr.str();
}

bool SpaceComposite::isEqual(const ASpace* space) const
{
  if (!ASpace::isEqual(space)) return false;
  unsigned int nc = getNComponents();
  for (unsigned int idx = 0; idx < nc; idx++)
  {
    const auto* c1 = getComponent(idx);
    const auto* c2 = space->getComponent(idx);
    if (!c1->isEqual(c2)) return false;
  }
  return true;
}

VectorDouble SpaceComposite::getDistances(const SpacePoint& p1,
                                          const SpacePoint& p2) const
{
  VectorDouble dis;
  if (p1.getNDim() != p2.getNDim())
  {
    std::cout << "Error: Inconsistent point dimension. Return empty distances"
              << std::endl;
    return dis;
  }
  for (auto* sp: _comps)
  {
    dis.push_back(sp->getDistance(p1, p2));
  }
  return dis;
}

/////////////////////////////////////////////////////////

void SpaceComposite::addSpaceComponent(const ASpace* comp)
{
  ASpace* sp = dynamic_cast<ASpace*>(comp->clone());
  sp->ASpace::setOffset(getNDim()); // TODO : I want this to be private and me a friend of ASpace
  _comps.push_back(sp);
  _nDim += sp->getNDim();
  const VectorDouble& o = sp->getOrigin(0);
  _origin.insert(_origin.end(), o.begin(), o.end());
  _work1.resize(_nDim);
  _work2.resize(_nDim);
}

/////////////////////////////////////////////////////////

void SpaceComposite::_move(SpacePoint& p1, const VectorDouble& vec) const
{
  for (auto* sp: _comps)
  {
    sp->move(p1, vec);
  }
}

/// Return the distance between two space points
double SpaceComposite::_getDistance(const SpacePoint& p1,
                                    const SpacePoint& p2,
                                    int ispace) const
{
  if (ispace < 0 || ispace >= (int)getNComponents())
    return getDistances(p1, p2).norm(); // Return the norm of sub-distances vector
  return _comps[ispace]->getDistance(p1, p2);
}

/// Return the distance between two space points with the given tensor
double SpaceComposite::_getDistance(const SpacePoint& p1,
                                    const SpacePoint& p2,
                                    const Tensor& tensor,
                                    int ispace) const
{
  if (ispace < 0 || ispace >= (int)getNComponents())
  {
    std::cout << "Error: Inconsistent space dimension. Return TEST."
              << std::endl;
    return TEST;
  }
  return _comps[ispace]->getDistance(p1, p2, tensor);
}

/// Return the distance in frequential domain between two space points with the
/// given tensor
double SpaceComposite::_getFrequentialDistance(const SpacePoint& p1,
                                              const SpacePoint& p2,
                                              const Tensor& tensor,
                                              int ispace) const
{
  if (ispace < 0 || ispace >= (int)getNComponents())
  {
    std::cout << "Error: Inconsistent space dimension. Return TEST."
              << std::endl;
    return TEST;
  }
  return _comps[ispace]->getFrequentialDistance(p1, p2, tensor);
}

/// Return the increment vector between two space points
VectorDouble SpaceComposite::_getIncrement(const SpacePoint& p1,
                                           const SpacePoint& p2,
                                           int ispace) const
{
  _getIncrementInPlace(p1, p2, _work1, ispace);
  return _work1;
}

/// Return the increment vector between two space points in a given vector
void SpaceComposite::_getIncrementInPlace(const SpacePoint& p1,
                                          const SpacePoint& p2,
                                          VectorDouble& ptemp,
                                          int ispace) const
{
  ptemp.clear();
  if (ispace < 0 || ispace >= (int)getNComponents())
  {
    for (auto* sp: _comps)
    {
      VectorDouble inc = sp->getIncrement(p1, p2);
      ptemp.insert(ptemp.begin(), inc.begin(), inc.end());
    }
  }
  else
  {
    ptemp = _comps[ispace]->getIncrement(p1, p2);
  }
}

/// Destroy components and empty vector
void SpaceComposite::_destroyComponents()
{
  for (auto* c: _comps)
  {
    delete c;
  }
}