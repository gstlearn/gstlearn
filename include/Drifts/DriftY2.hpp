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

class GSTLEARN_EXPORT DriftY2 : public ADriftElem
{
public:
  DriftY2(const CovContext& ctxt);
  DriftY2(const DriftY2 &r);
  DriftY2& operator= (const DriftY2 &r);
  virtual ~DriftY2();

  IClonable* clone() const override;

  String getDriftSymbol() const override { return "y2"; }
  String getDriftName() const override { return "Drift Y^2"; }
  int getOrderIRF() const override { return 2; }
  double eval(const Db* db, int iech) const override;
};

