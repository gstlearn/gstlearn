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

class CovGC3 : public ACovFunc
{
public:
  CovGC3(const CovContext& ctx);
  CovGC3(const CovGC3 &r);
  CovGC3& operator= (const CovGC3 &r);
  virtual ~CovGC3();

  int    hasRange() const override { return -1; }
  unsigned int getMinOrder()  const override { return 1; }
  String         getCovName() const override { return "Order-3 G.C."; }

protected:
  double _evaluateCov(double h) const override;
};

