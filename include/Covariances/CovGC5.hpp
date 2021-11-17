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

class CovGC5 : public ACovFunc
{
public:
  CovGC5(const CovContext& ctx);
  CovGC5(const CovGC5 &r);
  CovGC5& operator= (const CovGC5 &r);
  virtual ~CovGC5();

  int    hasRange() const override { return -1; }
  int    getMinOrder()  const override { return 2; }
  String getCovName() const override { return "Order-5 G.C."; }

protected:
  double _evaluateCov(double h) const override;
};

