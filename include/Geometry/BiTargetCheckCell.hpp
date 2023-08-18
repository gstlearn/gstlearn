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
#pragma once

#include <Geometry/ABiTargetCheck.hpp>
#include "gstlearn_export.hpp"

#include "Faults/Faults.hpp"

class GSTLEARN_EXPORT BiTargetCheckCell: public ABiTargetCheck
{
public:
  BiTargetCheckCell(const DbGrid* dbgrid = nullptr);
  BiTargetCheckCell(const BiTargetCheckCell& r);
  BiTargetCheckCell& operator=(const BiTargetCheckCell& r);
  virtual ~BiTargetCheckCell();

  virtual bool isOK(const SpaceTarget &T1,
                    const SpaceTarget &T2) const override;
  virtual bool isValid(const Db* dbin, const Db* dbout) override;

  /// ICloneable Interface
  IMPLEMENT_CLONING(BiTargetCheckCell)

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static BiTargetCheckCell* create(const DbGrid* dbgrid = nullptr);

private:
  const DbGrid* _dbgrid;
};
