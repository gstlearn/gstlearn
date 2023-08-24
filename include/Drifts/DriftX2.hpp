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

class GSTLEARN_EXPORT DriftX2 : public ADriftElem
{
public:
  DriftX2(const CovContext& ctxt = CovContext());
  DriftX2(const DriftX2 &r);
  DriftX2& operator= (const DriftX2 &r);
  virtual ~DriftX2();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftX2)

  String getDriftSymbol() const override { return "x2"; }
  String getDriftName() const override { return "Drift X^2"; }
  int getOrderIRF() const override { return 2; }
  int getNDim() const { return 1; }
  double eval(const Db* db, int iech) const override;
};

