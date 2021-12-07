/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Drifts/ADriftElem.hpp"

class GSTLEARN_EXPORT Drift1 : public ADriftElem
{
public:
  Drift1(const CovContext& ctxt);
  Drift1(const Drift1 &r);
  Drift1& operator= (const Drift1 &r);
  virtual ~Drift1();

  IClonable* clone() const override;

  String getDriftSymbol() const override { return "1"; }
  String getDriftName() const override { return "Universality Condition"; }
  int getOrderIRF() const override { return 0; }
  double eval(const Db* db, int iech) const override;
};
