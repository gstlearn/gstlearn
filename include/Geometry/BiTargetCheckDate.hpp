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

class GSTLEARN_EXPORT BiTargetCheckDate: public ABiTargetCheck
{
public:
  BiTargetCheckDate(double deltamin, double deltamax);
  BiTargetCheckDate(const BiTargetCheckDate& r);
  BiTargetCheckDate& operator=(const BiTargetCheckDate& r);
  virtual ~BiTargetCheckDate();

  /// ICloneable Interface
  IMPLEMENT_CLONING(BiTargetCheckDate)

  virtual bool isOK(const SpaceTarget &T1, const SpaceTarget &T2) const override;

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static BiTargetCheckDate* create(double deltamin, double deltamax);

private:
  double _deltaMin;
  double _deltaMax;
};
