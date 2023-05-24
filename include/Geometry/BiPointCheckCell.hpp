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
#pragma once

#include "gstlearn_export.hpp"

#include "Geometry/ABiPointCheck.hpp"
#include "Faults/Faults.hpp"

class GSTLEARN_EXPORT BiPointCheckCell: public ABiPointCheck
{
public:
  BiPointCheckCell(const DbGrid* dbgrid);
  BiPointCheckCell(const BiPointCheckCell& r);
  BiPointCheckCell& operator=(const BiPointCheckCell& r);
  virtual ~BiPointCheckCell();

  virtual bool isOK(const SpacePoint &P1,
                    const SpacePoint &P2,
                    int iech1 = -1,
                    int iech2 = -1) const override;

  /// ICloneable Interface
  IMPLEMENT_CLONING(BiPointCheckCell)

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

private:
  const DbGrid* _dbgrid;
};
