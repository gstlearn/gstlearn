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

class GSTLEARN_EXPORT DriftX2Y : public ADriftElem
{
public:
  DriftX2Y(const CovContext& ctxt = CovContext());
  DriftX2Y(const DriftX2Y &r);
  DriftX2Y& operator= (const DriftX2Y &r);
  virtual ~DriftX2Y();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftX2Y)

  String getDriftSymbol() const override { return "x2y"; }
  String getDriftName() const override { return "Drift X^2Y"; }
  int getOrderIRF() const override { return 3; }
  int getNDim() const { return 2; }
  double eval(const Db* db, int iech) const override;
};

