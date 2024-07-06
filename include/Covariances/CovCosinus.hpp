/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Covariances/ACovFunc.hpp"

class CovContext;

class GSTLEARN_EXPORT CovCosinus : public ACovFunc
{
public:
  CovCosinus(const CovContext& ctx);
  CovCosinus(const CovCosinus &r);
  CovCosinus& operator= (const CovCosinus &r);
  virtual ~CovCosinus();

  unsigned int getMaxNDim() const override { return 1; }

  String         getCovName() const override { return "Cosinus"; }
  int            getMinOrder() const override { return -1; }
  bool           getCompatibleSpaceR() const override { return true; }

protected:
  double _evaluateCov(double h) const override;
};

