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

namespace gstlrn
{
  ProjComposition::ProjComposition(std::vector<const IProj*> projs)
    : IProj()
  {
    // Nothing to do except checking compatibility. Incompatible sizes are a
    // developer bug and return codes of other functions are never really
    // checked, but at least we can say we're doing it!

    // Two special cases: empty list is treated later as an error and a single
    // element is weird but why not (but it's just handled normally below)?
    if (projs.size() == 0) return;

    // Check compatibility. Use raw pointers to make iterating easier.
    size_t idx = 0;
    const IProj* p1 = projs[idx++];
    while (idx < projs.size())
    {
      // Check that this operator is compatible with the previous one.
      const IProj* p2 = projs[idx];
      if (p1->getNPoint() != p2->getNApex())
      {
        // Abort.
        return;
      }
      p1 = p2;
      idx++;
    }
    for (const auto* p: projs) _projs.emplace_back(*p);
  }

  Id ProjComposition::setWorkArrays(vect work1, vect work2)
  {
    // Set work arrays (typically a space reused by other operators to save
    // memory). This function just checks that they are large enough.

    if (_projs.size() < 2) return 0; // no work array needed

    // Max of all intermediate sizes. One of the two must be that size, the
    // other (if needed) could be smaller but the proper size depends on how all
    // operators are chained (it might not be the second largest size!). Usually
    // both sizes will not be too different anyway so use same for both.
    size_t max_size = 0;
    for (size_t i = 0; i < _projs.size() - 1; ++i)
    {
      max_size = std::max(max_size, (size_t)_projs[i].get().getNPoint());
    }

    if (work1.size() < max_size
        || (_projs.size() > 2 && work2.size() < max_size))
      return -1;

    _work1 = work1;
    _work2 = work2;
    return 0;
  }

  Id ProjComposition::initWorkArrays(vect& work1, vect& work2) const
  {
    if (_projs.size() < 2) return 0; // no work array needed

    size_t max_size = 0;
    for (size_t i = 0; i < _projs.size() - 1; ++i)
    {
      max_size = std::max(max_size, (size_t)_projs[i].get().getNPoint());
    }

    if (_work1.size() > 0)
    {
      if (_work1.size() < max_size) return -1;
      work1 = _work1;
    }
    else
    {
      if (_w1.size() < max_size) _w1.resize(max_size);
      work1 = _w1;
    }

    if (_projs.size() == 2) return 0;

    if (_work2.size() > 0)
    {
      if (_work2.size() < max_size) return -1;
      work2 = _work2;
    }
    else
    {
      if (_w2.size() < max_size) _w2.resize(max_size);
      work2 = _w2;
    }

    return 0;
  }

  Id ProjComposition::_addPoint2mesh(const constvect in, vect out) const
  {
    if (_projs.size() == 0) return -1;
    if (_projs.size() == 1) return _projs[0].get().addPoint2mesh(in, out);

    vect work1, work2;
    Id ret = initWorkArrays(work1, work2);
    if (ret != 0) return ret;

    // Call point2mesh() to initialise temporary results to 0, but use
    // addPoint2Mesh() at the end to preserve what's already in 'out'.
    size_t idx = _projs.size() - 1;

    // Unroll a bit the loop because first/last use different in/out arrays.
    // This also means size == 1 is a special case (treated above).
    ret = _projs[idx].get().point2mesh(in, work1);
    if (ret != 0) return ret;
    idx--;

    vect win = work1, wout = work2;
    while (idx > 0)
    {
      ret = _projs[idx].get().point2mesh(win, wout);
      if (ret != 0) return ret;
      idx--;
      std::swap(win, wout);
    }

    return _projs[idx].get().addPoint2mesh(win, out);
  }

  Id ProjComposition::_addMesh2point(const constvect in, vect out) const
  {
    if (_projs.size() == 0) return -1;
    if (_projs.size() == 1) return _projs[0].get().addMesh2point(in, out);

    vect work1, work2;
    Id ret = initWorkArrays(work1, work2);
    if (ret != 0) return ret;

    size_t idx = 0;
    ret = _projs[idx].get().mesh2point(in, work1);
    if (ret != 0) return ret;
    idx++;

    vect win = work1, wout = work2;
    while (idx < _projs.size() - 1)
    {
      ret = _projs[idx].get().mesh2point(win, wout);
      if (ret != 0) return ret;
      idx++;
      std::swap(win, wout);
    }

    return _projs[idx].get().addMesh2point(win, out);
  }
} // namespace gstlrn
