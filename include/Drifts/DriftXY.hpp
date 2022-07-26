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

class GSTLEARN_EXPORT DriftXY : public ADriftElem
{
public:
  DriftXY(const CovContext& ctxt);
  DriftXY(const DriftXY &r);
  DriftXY& operator= (const DriftXY &r);
  virtual ~DriftXY();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftXY)

  String getDriftSymbol() const override { return "xy"; }
  String getDriftName() const override { return "Drift XY"; }
  int getOrderIRF() const override { return 2; }
  int getNDim() const { return 2; }
  double eval(const Db* db, int iech) const override;
};

