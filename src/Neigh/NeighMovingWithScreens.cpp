/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "Neigh/NeighMovingWithScreens.hpp"

NeighMovingWithScreens::NeighMovingWithScreens(
  int ndim, 
  const VectorDouble &screens,
  bool flag_xvalid
) : NeighMoving(ndim, flag_xvalid)
{
  if (! screens.empty()) setScreens(screens);
}

int NeighMovingWithScreens::reset(
  int ndim,
  bool flag_xvalid,
  const VectorDouble &screens,
  int nmaxi,
  double radius,
  int nmini,
  int nsect,
  int nsmax,
  VectorDouble coeffs,
  VectorDouble angles,
  double distcont
)
{
  auto res = NeighMoving::reset(
      ndim, flag_xvalid, nmaxi, radius, nmini, nsect, nsmax, coeffs, angles, distcont
  );
  if (res != 0) return res;
  setScreens(screens);
  return 0;
}

NeighMovingWithScreens* NeighMovingWithScreens::create(
  int ndim,
  bool flag_xvalid,
  const VectorDouble &screens,
  int nmaxi,
  double radius,
  int nmini,
  int nsect,
  int nsmax,
  VectorDouble coeffs,
  VectorDouble angles,
  double distcont
)
{
  NeighMovingWithScreens* neighM = new NeighMovingWithScreens();
  if (neighM->reset(ndim, flag_xvalid, screens, nmaxi, radius, nmini,
                    nsect, nsmax, coeffs, angles, distcont))
  {
    messerr("Problem when creating Moving NeighMovingWithScreens");
    delete neighM;
    neighM = nullptr;
  }
  return neighM;
}

/**
 * Sets the list of screens. Each screen is a 2D segment defined by two points.
 * @param screens All segments (pairs of 2D points) concatenated into a single
 * array. Vector size must be a multiple of 4 (2 * 2) to be valid
 */
void NeighMovingWithScreens::setScreens(const VectorDouble &screens)
{
  constexpr auto ndim = 2; // FIXME Consider N-D screens?
  constexpr auto size_seg = 2 * ndim;
  if (screens.size() % size_seg != 0) {
    messerr(
      "Invalid number of input coordinates (%i). Should be a multiple of (2 * %i)",
      screens.size(), ndim
    );
    return;
  }
  _screens = VectorScreens();
  _screens.reserve(screens.size() / size_seg);
  for (auto itr = screens.cbegin(); itr != screens.cend(); itr += size_seg) {
    auto p0 = VectorDouble(itr, itr + ndim);
    auto p1 = VectorDouble(itr + ndim, itr + size_seg);
    if (p0 != p1) _screens.emplace_back(p0, p1);
  }
}

// TODO Still to implement: reuse NeighMoving methods and add _screen
// String NeighMovingWithScreens::toString(const AStringFormat* strfmt) const
// int NeighMovingWithScreens::dumpToNF(const String& neutralFilename, bool verbose) const
// NeighMovingWithScreens* NeighMovingWithScreens::createFromNF(const String& neutralFilename, bool verbose)
// int NeighMovingWithScreens::_deserialize(std::istream& is, bool verbose)
// int NeighMovingWithScreens::_serialize(std::ostream& os, bool verbose) const
