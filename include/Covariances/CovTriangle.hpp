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

class GSTLEARN_EXPORT CovTriangle : public ACovFunc
{
public:
  CovTriangle(const CovContext& ctx);
  CovTriangle(const CovTriangle &r);
  CovTriangle& operator= (const CovTriangle &r);
  virtual ~CovTriangle();

  unsigned int   getMaxNDim()   const  override { return 1; }

  String         getCovName() const override { return "Triangle"; }
  int            getMinOrder() const override { return -1; }
  bool           getCompatibleSpaceR() const override { return true; }

protected:
  double _evaluateCov(double h) const override;
};

