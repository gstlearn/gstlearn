/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include <Geometry/ABiTargetCheck.hpp>
#include "gstlearn_export.hpp"

#include "Faults/Faults.hpp"

class GSTLEARN_EXPORT BiTargetCheckFaults: public ABiTargetCheck
{
public:
  BiTargetCheckFaults(const Faults* faults);
  BiTargetCheckFaults(const BiTargetCheckFaults& r);
  BiTargetCheckFaults& operator=(const BiTargetCheckFaults& r);
  virtual ~BiTargetCheckFaults();

  /// ICloneable Interface
  IMPLEMENT_CLONING(BiTargetCheckFaults)

  virtual bool isOK(const SpaceTarget &T1,
                    const SpaceTarget &T2) const override;

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  const Faults* getFaults() const { return _faults; }

private:
  const Faults* _faults;
};
