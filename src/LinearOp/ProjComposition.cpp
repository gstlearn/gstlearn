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
#include "LinearOp/ProjComposition.hpp"

ProjComposition::ProjComposition(std::vector<const IProj*> projs) : IProj()
{
  // Nothing to do except checking compatibility. Incompatible sizes are a
  // developer bug and return codes of other functions are never really
  // checked, but at least we can say we're doing it!

  // Two special cases: empty list is treated later as an error and a single
  // element is weird but why not (but it's just handled normally below)?
  if (projs.size() == 0) return;

  // Check compatibility and allocate work arrays. Use raw pointers to make
  // iterating easier then if everything is OK convert to std::unique_ptr<>.
  _works.resize(projs.size()-1);
  int idx = 0;
  const IProj* p1 = projs[idx++];
  while (idx < projs.size())
  {
    // Check that this operator is compatible with the previous one.
    const IProj* p2 = projs[idx];
    if (p1->getNPoint() != p2->getNApex())
    {
      // Abort. Delete everything first.
      for (auto p : projs) delete p;
      return;
    }
    _works[idx-1].resize(p1->getNPoint());
    p1 = p2;
    idx++;
  }
  for (auto p : projs) _projs.push_back(std::unique_ptr<const IProj>(p));
}

int ProjComposition::_addPoint2mesh(const constvect in, vect out) const
{
  if (_projs.size() == 0) return -1;
  if (_projs.size() == 1) return _projs[0]->addPoint2mesh(in, out);

  // Call point2mesh() on 'in' to initialise the temporary result to 0, but
  // use addPoint2Mesh() for the rest to preserve what's already in 'out'.
  size_t idx = _projs.size()-1;

  // Unroll a bit the loop because first/last use different in/out arrays.
  // This also means size == 1 is a special case (treated above).
  int ret = _projs[idx]->point2mesh(in, _works[idx-1]);
  if (ret != 0) return ret;
  idx--;

  while (idx > 0)
  {
    ret = _projs[idx]->addPoint2mesh(_works[idx], _works[idx-1]);
    if (ret != 0) return ret;
    idx--;
  }

  return _projs[idx]->addPoint2mesh(_works[0], out);
}

int ProjComposition::_addMesh2point(const constvect in, vect out) const {
  if (_projs.size() == 0) return -1;
  if (_projs.size() == 1) return _projs[0]->addMesh2point(in, out);

  size_t idx = 0;

  int ret = _projs[0]->mesh2point(in, _works[0]);
  if (ret != 0) return ret;
  idx++;

  while (idx < _projs.size() - 1)
  {
    ret = _projs[idx]->addMesh2point(_works[idx-1], _works[idx]);
    if (ret != 0) return ret;
    idx++;
  }

  return _projs[idx]->addMesh2point(_works[idx-1], out);
}
