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

class GSTLEARN_EXPORT DriftX3 : public ADriftElem
{
public:
  DriftX3(const CovContext& ctxt);
  DriftX3(const DriftX3 &r);
  DriftX3& operator= (const DriftX3 &r);
  virtual ~DriftX3();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftX3)

  String getDriftSymbol() const override { return "x3"; }
  String getDriftName() const override { return "Drift X^3"; }
  int getOrderIRF() const override { return 3; }
  int getNDim() const { return 1; }
  double eval(const Db* db, int iech) const override;
};

