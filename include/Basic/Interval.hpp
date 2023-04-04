/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"

class GSTLEARN_EXPORT Interval : public AStringable
{
public:
  Interval(double vmin = TEST,
           double vmax = TEST,
           bool   mininc = true,
           bool   maxinc = false);
  Interval(const Interval &m);
  Interval& operator=(const Interval &m);
  virtual ~Interval();

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  void init(double vmin = TEST,
            double vmax = TEST,
            bool mininc = true,
            bool maxinc = false);
  bool isMinDefined() const { return (! FFFF(_vmin)); }
  bool isMaxDefined() const { return (! FFFF(_vmax)); }
  bool isInside(double value) const;
  bool isOutsideBelow(double value) const;
  bool isOutsideAbove(double value) const;

  double getVmax() const { return _vmax; }
  void   setVmax(double vmax) { _vmax = vmax; }
  double getVmin() const { return _vmin; }
  void   setVmin(double vmin) { _vmin = vmin; }
  bool   getMaxIncluded() const { return _maxIncluded; }
  void   setMaxIncluded(bool maxIncluded) { _maxIncluded = maxIncluded; }
  bool   getMinIncluded() const { return _minIncluded; }
  void   setMinIncluded(bool minIncluded) { _minIncluded = minIncluded; }

  VectorDouble getBounds() const;

  bool isValid() const;
  bool isDisjoint(const Interval& m) const;

private:
  void _modifyUnbounded();
  bool _isValidInterval(void);

private:
  double _vmin;
  double _vmax;
  bool _minIncluded;
  bool _maxIncluded;
};
