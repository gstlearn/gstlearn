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
#include "Geometry/BiTargetCheckCell.hpp"
#include "Db/DbGrid.hpp"
#include "Space/SpaceTarget.hpp"

BiTargetCheckCell::BiTargetCheckCell(const DbGrid* dbgrid)
    : ABiTargetCheck(),
      _dbgrid(dbgrid)
{
}

BiTargetCheckCell::BiTargetCheckCell(const BiTargetCheckCell &r)
    : ABiTargetCheck(r),
      _dbgrid(r._dbgrid)
{
}

BiTargetCheckCell& BiTargetCheckCell::operator=(const BiTargetCheckCell &r)
{
  if (this != &r)
  {
    ABiTargetCheck::operator=(r);
    _dbgrid = r._dbgrid;
  }
  return *this;
}

BiTargetCheckCell::~BiTargetCheckCell()
{
}

BiTargetCheckCell* BiTargetCheckCell::create(const DbGrid* dbgrid)
{
  return new BiTargetCheckCell(dbgrid);
}

String BiTargetCheckCell::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << "Reject samples which do not belong to target Block" << std::endl;

  return sstr.str();
}

bool BiTargetCheckCell::isOK(const SpaceTarget &T1,
                             const SpaceTarget &T2) const
{
  // Check if the sample belongs to the cell
  bool valOK = _dbgrid->getGrid().sampleBelongsToCell(T2.getCoord(), T1.getCoord(), T1.getExtend());
  return valOK;
}

bool BiTargetCheckCell::isValid(const Db* dbin, const Db* dbout)
{
  DECLARE_UNUSED(dbin);
  if (! dbout->isGrid()) return false;
  _dbgrid = dynamic_cast<const DbGrid*>(dbout);
  return true;
}
