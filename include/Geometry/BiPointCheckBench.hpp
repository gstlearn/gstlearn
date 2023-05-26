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
#include "Faults/Faults.hpp"

class GSTLEARN_EXPORT BiPointCheckBench: public ABiPointCheck
{
public:
  BiPointCheckBench(int idim_bench, double width);
  BiPointCheckBench(const BiPointCheckBench& r);
  BiPointCheckBench& operator=(const BiPointCheckBench& r);
  virtual ~BiPointCheckBench();

  virtual bool isOK(const SpacePoint &P1,
                    const SpacePoint &P2,
                    int iech1 = -1,
                    int iech2 = -1) const override;

  /// ICloneable Interface
  IMPLEMENT_CLONING(BiPointCheckBench)

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static BiPointCheckBench* create(int idim_bench, double width);

  double getWidth() const { return _width; }

private:
  int       _idimBench;
  double    _width;
};
