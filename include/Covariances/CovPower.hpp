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

class CovPower : public ACovFunc
{
public:
  CovPower(const CovContext& ctx);
  CovPower(const CovPower &r);
  CovPower& operator= (const CovPower &r);
  virtual ~CovPower();

  int          hasRange()    const override { return -1; }
  bool         hasParam()    const override { return true; }
  double       getParMax()   const override { return 1.99; }
  unsigned int getMinOrder() const override { return 0; }
  String       getCovName()  const override { return "Power"; }

protected:
  double _evaluateCov(double h) const override;
};

