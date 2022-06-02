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

class GSTLEARN_EXPORT CovMarkov : public ACovFunc
{
public:
  CovMarkov(const CovContext& ctx);
  CovMarkov(const CovMarkov &r);
  CovMarkov& operator= (const CovMarkov &r);
  virtual ~CovMarkov();

  virtual String getFormula() const override;
  String         getCovName() const override { return "K-Bessel"; }


  bool   hasParam() const override { return true; }
  double getParMax() const override { return MAX_PARAM; }
  double getScadef() const override;
  bool   hasCovOnSphere() const override { return true; }
  bool   hasSpectrum() const { return true; }

protected:
  double _evaluateCov(double h) const override;
  double _evaluateCovOnSphere(double scale, int degree = 50) const override;
  double _evaluateSpectrum(double freq, double scale, int ndim) const override;
};
