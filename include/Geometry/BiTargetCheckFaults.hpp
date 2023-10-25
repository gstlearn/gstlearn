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

#include "gstlearn_export.hpp"

#include "Geometry/ABiTargetCheck.hpp"
#include "Faults/Faults.hpp"

class GSTLEARN_EXPORT BiTargetCheckFaults: public ABiTargetCheck
{
public:
  BiTargetCheckFaults(const Faults* faults);
  BiTargetCheckFaults(const BiTargetCheckFaults& r);
  BiTargetCheckFaults& operator=(const BiTargetCheckFaults& r);
  virtual ~BiTargetCheckFaults();

  /// ICloneable Interface
  //IMPLEMENT_CLONING(BiTargetCheckFaults)

  virtual bool isOK(const SpaceTarget &T1,
                    const SpaceTarget &T2) const override;

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static BiTargetCheckFaults* create(const Faults* faults);

private:
  const Faults* _faults;
};
