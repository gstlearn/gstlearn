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

class GSTLEARN_EXPORT CovNugget : public ACovFunc
{
public:
  CovNugget(const CovContext& ctx);
  CovNugget(const CovNugget &r);
  CovNugget& operator= (const CovNugget &r);
  virtual ~CovNugget();

  virtual String getFormula() const override;
  String         getCovName() const override { return "Nugget Effect"; }
  int            getMinOrder() const override { return -1; }
  bool           getCompatibleSpaceR() const override { return true; }

  int    hasRange() const override { return 0; }

  bool isValidForTurningBand() const override { return true; }

protected:
  double _evaluateCov(double h) const override;
};

