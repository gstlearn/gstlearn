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
  bool   hasMarkovCoeffs() const override { return true; }
  double evaluateSpectrum(double freq, int ndim) const override;
  void   setMarkovCoeffs(VectorDouble coeffs) override { _markovCoeffs = coeffs;}
  VectorDouble getMarkovCoeffs() const override;
  double getCorrec() const override { return _correc;}
  void   computeCorrec(int ndim);
  void   setCorrec(double val) override { _correc = val;}
  void   computeMarkovCoeffs(int dim) override;

protected:
  double _evaluateCov(double h) const override;
  double _evaluateCovOnSphere(double scale, int degree = 50) const override;


private:
  double _correc;
  VectorDouble _markovCoeffs;

};
