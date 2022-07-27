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

class GSTLEARN_EXPORT DriftZ : public ADriftElem
{
public:
  DriftZ(const CovContext& ctxt);
  DriftZ(const DriftZ &r);
  DriftZ& operator= (const DriftZ &r);
  virtual ~DriftZ();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftZ)

  String getDriftSymbol() const override { return "z"; }
  String getDriftName() const override { return "Drift Z"; }
  int getOrderIRF() const override { return 1; }
  int getNDim() const { return 3; }
  double eval(const Db* db, int iech) const override;
};

