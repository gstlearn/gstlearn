/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Covariances/ACovFunc.hpp"

/* Be careful ! This is not a real covariance */

class CovContext;

class GSTLEARN_EXPORT CovGCspline : public ACovFunc
{
public:
  CovGCspline(const CovContext& ctx);
  CovGCspline(const CovGCspline &r);
  CovGCspline& operator= (const CovGCspline &r);
  virtual ~CovGCspline();

  int            hasRange() const override { return -1; }
  String         getCovName() const override { return "Spline G.C."; }
  int            getMinOrder() const override { return 1; }

protected:
  double _evaluateCov(double h) const override;
};

