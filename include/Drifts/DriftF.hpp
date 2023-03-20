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

class GSTLEARN_EXPORT DriftF : public ADriftElem
{
public:
  DriftF(const CovContext& ctxt = CovContext());
  DriftF(const DriftF &r);
  DriftF& operator= (const DriftF &r);
  virtual ~DriftF();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftF)

  String getDriftSymbol() const override { return "f"; }
  String getDriftName() const override { return "External Drift"; }
  int getOrderIRF() const override { return 0; }
  bool getDriftExternal() const override { return true; }
  double eval(const Db* db, int iech) const override;
};

