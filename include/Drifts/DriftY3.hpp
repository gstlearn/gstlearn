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

class GSTLEARN_EXPORT DriftY3 : public ADriftElem
{
public:
  DriftY3(const CovContext& ctxt = CovContext());
  DriftY3(const DriftY3 &r);
  DriftY3& operator= (const DriftY3 &r);
  virtual ~DriftY3();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftY3)

  String getDriftSymbol() const override { return "y3"; }
  String getDriftName() const override { return "Drift Y^3"; }
  int getOrderIRF() const override { return 3; }
  int getNDim() const { return 2; }
  double eval(const Db* db, int iech) const override;
};

