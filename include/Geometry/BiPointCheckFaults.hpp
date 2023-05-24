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

class GSTLEARN_EXPORT BiPointCheckFaults: public ABiPointCheck
{
public:
  BiPointCheckFaults(const Faults* faults);
  BiPointCheckFaults(const BiPointCheckFaults& r);
  BiPointCheckFaults& operator=(const BiPointCheckFaults& r);
  virtual ~BiPointCheckFaults();

  virtual bool isOK(const SpacePoint& P1, const SpacePoint& P2) const override;
  virtual bool isDistanceCalculated() const override { return false; }
  virtual double getDistance() const override { return TEST; }

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

private:
  const Faults* _faults;
};
