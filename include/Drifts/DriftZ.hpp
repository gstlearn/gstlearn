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

class GSTLEARN_EXPORT DriftZ : public ADriftElem
{
public:
  DriftZ(const CovContext& ctxt = CovContext());
  DriftZ(const DriftZ &r);
  DriftZ& operator= (const DriftZ &r);
  virtual ~DriftZ();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftZ)

  String getDriftSymbol() const override { return "z"; }
  String getDriftName() const override { return "Drift Z"; }
  int getOrderIRF() const override { return 1; }
  int getNDim() const { return 3; }
  double eval(const Db* db, int iech) const override;
};

