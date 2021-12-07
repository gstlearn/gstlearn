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

class GSTLEARN_EXPORT CovGCspline2 : public ACovFunc
{
public:
  CovGCspline2(const CovContext& ctx);
  CovGCspline2(const CovGCspline2 &r);
  CovGCspline2& operator= (const CovGCspline2 &r);
  virtual ~CovGCspline2();

  int          getMinOrder()  const override { return 1; }
  unsigned int getMaxNDim()   const  override { return 3; }
  String       getCovName() const override { return "Spline-2 G.C."; }

protected:
  double _evaluateCov(double h) const override;
  double _evaluateCovDerivate(int degree, double h) const override;
};

