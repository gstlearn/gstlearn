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

class GSTLEARN_EXPORT CovSpherical : public ACovFunc
{
public:
  CovSpherical(const CovContext& ctx);
  CovSpherical(const CovSpherical &r);
  CovSpherical& operator= (const CovSpherical &r);
  virtual ~CovSpherical();

  unsigned int   getMaxNDim() const override { return 3; }
  String         getFormula() const override;
  String         getCovName() const override { return "Spherical"; }
  int            getMinOrder() const override { return -1; }

protected:
  double _evaluateCov(double h) const override;
};

