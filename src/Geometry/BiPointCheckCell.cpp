/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f.h"

#include "Geometry/BiPointCheckCell.hpp"
#include "Space/SpacePoint.hpp"

BiPointCheckCell::BiPointCheckCell(const DbGrid* dbgrid)
    : ABiPointCheck(),
      _dbgrid(dbgrid)
{
}

BiPointCheckCell::BiPointCheckCell(const BiPointCheckCell &r)
    : ABiPointCheck(r),
      _dbgrid(r._dbgrid)
{
}

BiPointCheckCell& BiPointCheckCell::operator=(const BiPointCheckCell &r)
{
  if (this != &r)
  {
    ABiPointCheck::operator=(r);
    _dbgrid = r._dbgrid;
  }
  return *this;
}

BiPointCheckCell::~BiPointCheckCell()
{
}

/**
 * Print the context of the Cell rejection
 * @param strfmt Printing format
 * @return String describing the Checker option
 *
 * @remark The printout is not performed here as this checker can only
 * @remark be instantiated after 'attach'. It is directly performed in NeighMoving.
 */
String BiPointCheckCell::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  return sstr.str();
}

bool BiPointCheckCell::isOK(const SpacePoint &P1,
                            const SpacePoint &P2,
                            int iech1,
                            int iech2) const
{
  // Get the coordinates of the sample
  VectorDouble coor = P2.getCoord();

  // Identify the dimensions of the cell
  VectorDouble dxsPerCell = _dbgrid->getBlockExtensions(iech1);

  // Check if the sample belongs to the cell
  return _dbgrid->getGrid().sampleBelongsToCell(coor, iech1, dxsPerCell);
}
