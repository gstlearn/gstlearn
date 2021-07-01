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

#include "Drifts/ADriftElem.hpp"

class DriftXY2 : public ADriftElem
{
public:
  DriftXY2(const CovContext& ctxt);
  DriftXY2(const DriftXY2 &r);
  DriftXY2& operator= (const DriftXY2 &r);
  virtual ~DriftXY2();

  IClonable* clone() const override;

  String getDriftSymbol() const override { return "xy2"; }
  String getDriftName() const override { return "Drift XY^2"; }
  int getOrderIRF() const override { return 3; }
  double eval(const Db* db, int iech) const override;
};

