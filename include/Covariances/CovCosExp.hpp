/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Covariances/ACovFunc.hpp"

class CovContext;

class GSTLEARN_EXPORT CovCosExp : public ACovFunc
{
public:
  CovCosExp(const CovContext& ctx);
  CovCosExp(const CovCosExp &r);
  CovCosExp& operator= (const CovCosExp &r);
  virtual ~CovCosExp();

  double getParMax() const override { return TEST; }
  bool   hasParam()  const override { return true; }
  double getScadef() const override;

  String getCovName() const override { return "Cosexp"; }
  int    getMinOrder() const override { return -1; }

protected:
  double _evaluateCov(double h) const override;
};

