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

class GSTLEARN_EXPORT DriftY3 : public ADriftElem
{
public:
  DriftY3(const CovContext& ctxt);
  DriftY3(const DriftY3 &r);
  DriftY3& operator= (const DriftY3 &r);
  virtual ~DriftY3();

  IClonable* clone() const override;

  String getDriftSymbol() const override { return "y3"; }
  String getDriftName() const override { return "Drift Y^3"; }
  int getOrderIRF() const override { return 3; }
  int getNDim() const { return 2; }
  double eval(const Db* db, int iech) const override;
};

