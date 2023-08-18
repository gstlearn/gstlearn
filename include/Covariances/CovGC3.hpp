/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Covariances/ACovFunc.hpp"

/* Be careful ! This is not a real covariance */

class CovContext;

class GSTLEARN_EXPORT CovGC3 : public ACovFunc
{
public:
  CovGC3(const CovContext& ctx);
  CovGC3(const CovGC3 &r);
  CovGC3& operator= (const CovGC3 &r);
  virtual ~CovGC3();

  int    hasRange() const override { return -1; }
  int    getMinOrder()  const override { return 1; }
  String getCovName() const override { return "Order-3 G.C."; }

protected:
  double _evaluateCov(double h) const override;
};

