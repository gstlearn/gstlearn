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

class GSTLEARN_EXPORT BiPointCheckDistance: public ABiPointCheck
{
public:
  BiPointCheckDistance(double dmax);
  BiPointCheckDistance(const BiPointCheckDistance& r);
  BiPointCheckDistance& operator=(const BiPointCheckDistance& r);
  virtual ~BiPointCheckDistance();

  virtual bool isOK(const SpacePoint& P1, const SpacePoint& P2) const override;
  virtual bool isDistanceCalculated() const override { return true; }
  virtual double getDistance() const override { return _dist; }

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  double getDmax() const { return _dmax; }
  double getDist() const { return _dist; }

private:
  double _dmax;

  mutable double _dist;
};
