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

class GSTLEARN_EXPORT DriftY2 : public ADriftElem
{
public:
  DriftY2(const CovContext& ctxt = CovContext());
  DriftY2(const DriftY2 &r);
  DriftY2& operator= (const DriftY2 &r);
  virtual ~DriftY2();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftY2)

  String getDriftSymbol() const override { return "y2"; }
  String getDriftName() const override { return "Drift Y^2"; }
  int getOrderIRF() const override { return 2; }
  int getNDim() const { return 2; }
  double eval(const Db* db, int iech) const override;
};

