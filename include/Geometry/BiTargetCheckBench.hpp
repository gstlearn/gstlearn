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
class GSTLEARN_EXPORT BiTargetCheckBench: public ABiTargetCheck
{
public:
  BiTargetCheckBench(Id idim_bench, double width);
  BiTargetCheckBench(const BiTargetCheckBench& r);
  BiTargetCheckBench& operator=(const BiTargetCheckBench& r);
  virtual ~BiTargetCheckBench();

  bool isOK(const SpaceTarget &T1, const SpaceTarget &T2) const override;
  bool isValid(const Db* dbin, const Db* dbout) override;

  /// ICloneable Interface
  IMPLEMENT_CLONING(BiTargetCheckBench)

  /// Interface to AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  static BiTargetCheckBench* create(Id idim_bench, double width);

  double getWidth() const { return _width; }

private:
  Id       _idimBench;
  double    _width;
};
}