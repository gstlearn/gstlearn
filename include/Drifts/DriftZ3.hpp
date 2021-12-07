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

class GSTLEARN_EXPORT DriftZ3 : public ADriftElem
{
public:
  DriftZ3(const CovContext& ctxt);
  DriftZ3(const DriftZ3 &r);
  DriftZ3& operator= (const DriftZ3 &r);
  virtual ~DriftZ3();

  IClonable* clone() const override;

  String getDriftSymbol() const override { return "z3"; }
  String getDriftName() const override { return "Drift Z^3"; }
  int getOrderIRF() const override { return 3; }
  double eval(const Db* db, int iech) const override;
};

