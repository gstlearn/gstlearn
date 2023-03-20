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

class CovContext;

class GSTLEARN_EXPORT CovBesselJ : public ACovFunc
{
public:
  CovBesselJ(const CovContext& ctx);
  CovBesselJ(const CovBesselJ &r);
  CovBesselJ& operator= (const CovBesselJ &r);
  virtual ~CovBesselJ();

  bool         hasParam() const override { return true; }
  double       getParMax() const override { return 2; }

  virtual String getFormula() const override;
  String         getCovName() const override { return "J-Bessel"; }
  int            getMinOrder() const override { return -1; }

protected:
  double _evaluateCov(double h) const override;
};

