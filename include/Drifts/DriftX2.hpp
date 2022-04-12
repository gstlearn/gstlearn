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

class GSTLEARN_EXPORT DriftX2 : public ADriftElem
{
public:
  DriftX2(const CovContext& ctxt);
  DriftX2(const DriftX2 &r);
  DriftX2& operator= (const DriftX2 &r);
  virtual ~DriftX2();

  IClonable* clone() const override;

  String getDriftSymbol() const override { return "x2"; }
  String getDriftName() const override { return "Drift X^2"; }
  int getOrderIRF() const override { return 2; }
  int getNDim() const { return 1; }
  double eval(const Db* db, int iech) const override;
};

