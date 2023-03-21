/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Drifts/ADriftElem.hpp"

class GSTLEARN_EXPORT DriftXY2 : public ADriftElem
{
public:
  DriftXY2(const CovContext& ctxt = CovContext());
  DriftXY2(const DriftXY2 &r);
  DriftXY2& operator= (const DriftXY2 &r);
  virtual ~DriftXY2();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftXY2)

  String getDriftSymbol() const override { return "xy2"; }
  String getDriftName() const override { return "Drift XY^2"; }
  int getOrderIRF() const override { return 3; }
  int getNDim() const { return 2; }
  double eval(const Db* db, int iech) const override;
};

