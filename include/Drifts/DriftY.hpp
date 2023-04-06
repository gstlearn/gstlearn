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
#include "Drifts/ADriftElem.hpp"

class GSTLEARN_EXPORT DriftY : public ADriftElem
{
public:
  DriftY(const CovContext& ctxt = CovContext());
  DriftY(const DriftY &r);
  DriftY& operator= (const DriftY &r);
  virtual ~DriftY();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftY)

  String getDriftSymbol() const override { return "y"; }
  String getDriftName() const override { return "Drift Y"; }
  int getOrderIRF() const override { return 1; }
  int getNDim() const { return 2; }
  double eval(const Db* db, int iech) const override;
};

