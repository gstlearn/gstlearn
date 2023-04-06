/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
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

