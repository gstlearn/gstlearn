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

#include "Basic/AStringable.hpp"

class GSTLEARN_EXPORT ABiPointCheck: public AStringable
{
public:
  ABiPointCheck();
  ABiPointCheck(const ABiPointCheck& r);
  ABiPointCheck& operator=(const ABiPointCheck& r);
  virtual ~ABiPointCheck();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  virtual bool isOK(const SpacePoint& P1, const SpacePoint& P2) const = 0;
  virtual bool isDistanceCalculated() const { return false; }
  virtual double getDistance() const { return TEST; }
};
