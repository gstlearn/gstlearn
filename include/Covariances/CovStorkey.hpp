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

class GSTLEARN_EXPORT CovStorkey : public ACovFunc
{
public:
  CovStorkey(const CovContext& ctx);
  CovStorkey(const CovStorkey &r);
  CovStorkey& operator= (const CovStorkey &r);
  virtual ~CovStorkey();

  unsigned int getMaxNDim()   const  override { return 1; }

  virtual String getFormula() const override { return String("Equation not yet implemented"); }
  String         getCovName() const override { return "Storkey"; }
  int            getMinOrder() const override { return -1; }

protected:
  double _evaluateCov(double h) const override;
};

