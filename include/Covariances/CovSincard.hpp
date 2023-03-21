/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Covariances/ACovFunc.hpp"

class CovContext;

class GSTLEARN_EXPORT CovSincard : public ACovFunc
{
public:
  CovSincard(const CovContext& ctx);
  CovSincard(const CovSincard &r);
  CovSincard& operator= (const CovSincard &r);
  virtual ~CovSincard();

  double getScadef() const override;
  virtual String getFormula() const override;
  String         getCovName() const override { return "Cardinal Sine"; }
  int            getMinOrder() const override { return -1; }

protected:
  double _evaluateCov(double h) const override;
};

