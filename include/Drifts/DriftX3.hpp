/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Drifts/ADriftElem.hpp"

class GSTLEARN_EXPORT DriftX3 : public ADriftElem
{
public:
  DriftX3(const CovContext& ctxt = CovContext());
  DriftX3(const DriftX3 &r);
  DriftX3& operator= (const DriftX3 &r);
  virtual ~DriftX3();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftX3)

  String getDriftSymbol() const override { return "x3"; }
  String getDriftName() const override { return "Drift X^3"; }
  int getOrderIRF() const override { return 3; }
  int getNDim() const { return 1; }
  double eval(const Db* db, int iech) const override;
};

