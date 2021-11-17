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

#include "Covariances/ACovFunc.hpp"

/* Be careful ! This is not a real covariance */

class CovContext;

class CovGC1 : public ACovFunc
{
public:
  CovGC1(const CovContext& ctx);
  CovGC1(const CovGC1 &r);
  CovGC1& operator= (const CovGC1 &r);
  virtual ~CovGC1();

  int    hasRange()    const override { return -1; }
  int    getMinOrder() const override { return 0; }
  String getCovName()  const override { return "Order-1 G.C."; }

protected:
  double _evaluateCov(double h) const override;
};

