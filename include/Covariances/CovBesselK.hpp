/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
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
  int            getMinOrder() const override { return -1; }

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
