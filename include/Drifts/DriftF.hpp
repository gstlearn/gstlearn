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

class GSTLEARN_EXPORT DriftF : public ADriftElem
{
public:
  DriftF(const CovContext& ctxt);
  DriftF(const DriftF &r);
  DriftF& operator= (const DriftF &r);
  virtual ~DriftF();

  IClonable* clone() const override;

  String getDriftSymbol() const override { return "f"; }
  String getDriftName() const override { return "External Drift"; }
  int getOrderIRF() const override { return 0; }
  double eval(const Db* db, int iech) const override;
};

