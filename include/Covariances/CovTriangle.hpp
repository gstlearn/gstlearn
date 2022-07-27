/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
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

