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

class GSTLEARN_EXPORT BiTargetCheckBench: public ABiTargetCheck
{
public:
  BiTargetCheckBench(int idim_bench, double width);
  BiTargetCheckBench(const BiTargetCheckBench& r);
  BiTargetCheckBench& operator=(const BiTargetCheckBench& r);
  virtual ~BiTargetCheckBench();

  virtual bool isOK(const SpaceTarget &T1,
                    const SpaceTarget &T2) const override;
  virtual bool isValid(const Db* dbin, const Db* dbout) override;

  /// ICloneable Interface
  //IMPLEMENT_CLONING(BiTargetCheckBench)

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static BiTargetCheckBench* create(int idim_bench, double width);

  double getWidth() const { return _width; }

private:
  int       _idimBench;
  double    _width;
};
