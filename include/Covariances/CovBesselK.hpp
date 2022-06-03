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

class GSTLEARN_EXPORT CovBesselK : public ACovFunc
{
public:
  CovBesselK(const CovContext& ctx);
  CovBesselK(const CovBesselK &r);
  CovBesselK& operator= (const CovBesselK &r);
  virtual ~CovBesselK();

  virtual String getFormula() const override;
  String         getCovName() const override { return "K-Bessel"; }


  bool   hasParam() const override { return true; }
  double getParMax() const override { return MAX_PARAM; }
  double getScadef() const override;
  bool   hasCovOnSphere() const override { return true; }
  bool   hasSpectrum() const override { return true; }
  bool   hasMarkovCoeffs() const override { return false; } //TODO pas to True when _computeCoeffs will be implemented
  double evaluateSpectrum(double freq, double scale, int ndim) const override;
  void   setMarkovCoeffs(VectorDouble coeffs){ _markovCoeffs = coeffs;}
  VectorDouble getMarkovCoeffs() const;
protected:
  double _evaluateCov(double h) const override;
  double _evaluateCovOnSphere(double scale, int degree = 50) const override;
  void _computeCoeffs() const;

private:
  mutable VectorDouble _markovCoeffs;
  mutable bool _coeffsComputed;
};
