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

class GSTLEARN_EXPORT DriftZ3 : public ADriftElem
{
public:
  DriftZ3(const CovContext& ctxt = CovContext());
  DriftZ3(const DriftZ3 &r);
  DriftZ3& operator= (const DriftZ3 &r);
  virtual ~DriftZ3();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftZ3)

  String getDriftSymbol() const override { return "z3"; }
  String getDriftName() const override { return "Drift Z^3"; }
  int getOrderIRF() const override { return 3; }
  int getNDim() const { return 3; }
  double eval(const Db* db, int iech) const override;
};

