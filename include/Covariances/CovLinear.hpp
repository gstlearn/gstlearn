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
#include "Covariances/ACovFunc.hpp"

/* Be careful ! This is not a real covariance */

class CovContext;

class GSTLEARN_EXPORT CovLinear : public ACovFunc
{
public:
  CovLinear(const CovContext& ctx);
  CovLinear(const CovLinear &r);
  CovLinear& operator= (const CovLinear &r);
  virtual ~CovLinear();

  int    hasRange() const override { return -1; }
  int    getMinOrder()  const override { return 0; }
  String getCovName() const override { return "Linear"; }

protected:
  double _evaluateCov(double h) const override;
};

