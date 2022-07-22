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

class GSTLEARN_EXPORT DriftY : public ADriftElem
{
public:
  DriftY(const CovContext& ctxt);
  DriftY(const DriftY &r);
  DriftY& operator= (const DriftY &r);
  virtual ~DriftY();

  ICloneable* clone() const override;

  String getDriftSymbol() const override { return "y"; }
  String getDriftName() const override { return "Drift Y"; }
  int getOrderIRF() const override { return 1; }
  int getNDim() const { return 2; }
  double eval(const Db* db, int iech) const override;
};

