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
#include "geoslib_define.h"

#include "Anamorphosis/AnamContinuous.hpp"
#include "Basic/ASerializable.hpp"

class GSTLEARN_EXPORT AnamUser: public AnamContinuous
{
private:
  double (*_y2z_function)(double);
  double (*_z2y_function)(double);

public:
  AnamUser();
  AnamUser(const AnamUser &m);
  AnamUser& operator= (const AnamUser &m);
  virtual ~AnamUser();

  /// ICloneable Interface
  IMPLEMENT_CLONING(AnamUser)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// AAnam Interface
  const EAnam& getType() const override { return EAnam::fromKey("EXTERNAL"); }
  bool isChangeSupportDefined() const override { return false; }

  /// AnamContinuous Interface
  void   calculateMeanAndVariance() override;
  double transformToRawValue(double h) const override;
  double rawToTransformValue(double h) const override;

  void setY2zFunction(double (*y2z_function)(double)) { _y2z_function = y2z_function; }
  void setZ2yFunction(double (*z2y_function)(double)) { _z2y_function = z2y_function; }

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "AnamUser"; }
};
