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

namespace gstlrn
{
class CovContext;

class GSTLEARN_EXPORT CovNugget : public ACovFunc
{
public:
  CovNugget(const CovContext& ctx);
  CovNugget(const CovNugget &r);
  CovNugget& operator= (const CovNugget &r);
  virtual ~CovNugget();

  String getFormula() const override;
  String         getCovName() const override { return "Nugget Effect"; }
  Id            getMinOrder() const override { return -1; }
  bool           getCompatibleSpaceR() const override { return true; }

  Id    hasRange() const override { return 0; }

  bool isValidForTurningBand() const override { return true; }

protected:
  double _evaluateCov(double h) const override;
};

}