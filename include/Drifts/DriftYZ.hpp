/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Drifts/ADriftElem.hpp"

class GSTLEARN_EXPORT DriftYZ : public ADriftElem
{
public:
  DriftYZ(const CovContext& ctxt = CovContext());
  DriftYZ(const DriftYZ &r);
  DriftYZ& operator= (const DriftYZ &r);
  virtual ~DriftYZ();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftYZ)

  String getDriftSymbol() const override { return "yz"; }
  String getDriftName() const override { return "Drift YZ"; }
  int getOrderIRF() const override { return 2; }
  int getNDim() const { return 3; }
  double eval(const Db* db, int iech) const override;
};

