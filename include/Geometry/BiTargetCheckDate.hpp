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

  double getDeltaMax() const { return _deltaMax; }
  void setDeltaMax(double deltaMax) { _deltaMax = deltaMax; }
  double getDeltaMin() const { return _deltaMin; }
  void setDeltaMin(double deltaMin) { _deltaMin = deltaMin; }

private:
  double _deltaMin;
  double _deltaMax;
};
