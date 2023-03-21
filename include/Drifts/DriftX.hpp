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

class GSTLEARN_EXPORT DriftX : public ADriftElem
{
public:
  DriftX(const CovContext& ctxt = CovContext());
  DriftX(const DriftX &r);
  DriftX& operator= (const DriftX &r);
  virtual ~DriftX();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftX)

  String getDriftSymbol() const override { return "x"; }
  String getDriftName() const override { return "Drift X"; }
  int getOrderIRF() const override { return 1; }
  int getNDim() const { return 1; }
  double eval(const Db* db, int iech) const override;
};
