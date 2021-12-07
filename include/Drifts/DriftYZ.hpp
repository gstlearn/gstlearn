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

class GSTLEARN_EXPORT DriftYZ : public ADriftElem
{
public:
  DriftYZ(const CovContext& ctxt);
  DriftYZ(const DriftYZ &r);
  DriftYZ& operator= (const DriftYZ &r);
  virtual ~DriftYZ();

  IClonable* clone() const override;

  String getDriftSymbol() const override { return "yz"; }
  String getDriftName() const override { return "Drift YZ"; }
  int getOrderIRF() const override { return 2; }
  double eval(const Db* db, int iech) const override;
};

