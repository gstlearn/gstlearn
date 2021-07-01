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

class DriftX : public ADriftElem
{
public:
  DriftX(const CovContext& ctxt);
  DriftX(const DriftX &r);
  DriftX& operator= (const DriftX &r);
  virtual ~DriftX();

  IClonable* clone() const override;

  String getDriftSymbol() const override { return "x"; }
  String getDriftName() const override { return "Drift X"; }
  int getOrderIRF() const override { return 1; }
  double eval(const Db* db, int iech) const override;
};
