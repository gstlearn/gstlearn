/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Drifts/ADriftElem.hpp"

class GSTLEARN_EXPORT DriftXY : public ADriftElem
{
public:
  DriftXY(const CovContext& ctxt = CovContext());
  DriftXY(const DriftXY &r);
  DriftXY& operator= (const DriftXY &r);
  virtual ~DriftXY();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftXY)

  String getDriftSymbol() const override { return "xy"; }
  String getDriftName() const override { return "Drift XY"; }
  int getOrderIRF() const override { return 2; }
  int getNDim() const { return 2; }
  double eval(const Db* db, int iech) const override;
};

