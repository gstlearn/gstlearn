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

namespace gstlrn
{
class GSTLEARN_EXPORT BiTargetCheckCode: public ABiTargetCheck
{
public:
  BiTargetCheckCode(Id optcode=1, double tolcode=EPSILON6);
  BiTargetCheckCode(const BiTargetCheckCode& r);
  BiTargetCheckCode& operator=(const BiTargetCheckCode& r);
  virtual ~BiTargetCheckCode();

  /// ICloneable Interface
  IMPLEMENT_CLONING(BiTargetCheckCode)

  bool isOK(const SpaceTarget &T1, const SpaceTarget &T2) const override;

  /// Interface to AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  static BiTargetCheckCode* create(Id optcode=1, double tolcode=EPSILON6);

private:
  Id _optCode;
  double _tolCode;
};
}