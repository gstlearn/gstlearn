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

#include "Basic/AStringable.hpp"

namespace gstlrn
{
class GSTLEARN_EXPORT FracDesc: public AStringable
{
public:
  FracDesc();
  FracDesc(const FracDesc& r);
  FracDesc& operator=(const FracDesc& r);
  virtual ~FracDesc();

  /// Interface for AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  Id getNPoint() const { return static_cast<Id>(_x.size()); }

  Id getFamily() const { return _family; }
  void setFamily(Id family) { _family = family; }
  double getOrient() const { return _orient; }
  void setOrient(double orient) { _orient = orient; }
  double getXXF(Id i) const { return _x[i]; }
  double getYYF(Id i) const { return _y[i]; }
  void setXXF(Id i, double value) { _x[i] = value; }
  void setYYF(Id i, double value) { _y[i] = value; }

  void addPoint(double x, double y);
  double fractureExtension(double cote, double dcote) const;

private:
  Id _family;
  double _orient;
  VectorDouble _x;
  VectorDouble _y;
};
} // namespace gstlrn