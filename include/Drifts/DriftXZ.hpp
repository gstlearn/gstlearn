/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Drifts/ADriftElem.hpp"

class GSTLEARN_EXPORT DriftXZ : public ADriftElem
{
public:
  DriftXZ(const CovContext& ctxt = CovContext());
  DriftXZ(const DriftXZ &r);
  DriftXZ& operator= (const DriftXZ &r);
  virtual ~DriftXZ();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftXZ)

  String getDriftSymbol() const override { return "xz"; }
  String getDriftName() const override { return "Drift XZ"; }
  int getOrderIRF() const override { return 2; }
  int getNDim() const { return 3; }
  double eval(const Db* db, int iech) const override;
};

