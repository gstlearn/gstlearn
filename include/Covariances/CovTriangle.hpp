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

class GSTLEARN_EXPORT CovTriangle : public ACovFunc
{
public:
  CovTriangle(const CovContext& ctx);
  CovTriangle(const CovTriangle &r);
  CovTriangle& operator= (const CovTriangle &r);
  virtual ~CovTriangle();

  unsigned int   getMaxNDim()   const  override { return 1; }

  virtual String getFormula() const override { return String("Equation not yet implemented"); }
  String         getCovName() const override { return "Triangle"; }
  int            getMinOrder() const override { return -1; }

protected:
  double _evaluateCov(double h) const override;
};

