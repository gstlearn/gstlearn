/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Drifts/ADriftElem.hpp"

class GSTLEARN_EXPORT DriftZ2 : public ADriftElem
{
public:
  DriftZ2(const CovContext& ctxt = CovContext());
  DriftZ2(const DriftZ2 &r);
  DriftZ2& operator= (const DriftZ2 &r);
  virtual ~DriftZ2();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftZ2)

  String getDriftSymbol() const override { return "z2"; }
  String getDriftName() const override { return "Drift Z^2"; }
  int getOrderIRF() const override { return 2; }
  int getNDim() const { return 3; }
  double eval(const Db* db, int iech) const override;
};

