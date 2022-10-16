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

