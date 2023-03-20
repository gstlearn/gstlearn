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

class GSTLEARN_EXPORT Drift1 : public ADriftElem
{
public:
  Drift1(const CovContext& ctxt = CovContext());
  Drift1(const Drift1 &r);
  Drift1& operator= (const Drift1 &r);
  virtual ~Drift1();

  /// ICloneable interface
  IMPLEMENT_CLONING(Drift1)

  String getDriftSymbol() const override { return "1"; }
  String getDriftName() const override { return "Universality Condition"; }
  int getOrderIRF() const override { return 0; }
  double eval(const Db* db, int iech) const override;
};
