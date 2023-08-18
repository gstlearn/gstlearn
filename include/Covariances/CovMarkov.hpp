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

class GSTLEARN_EXPORT CovMarkov : public ACovFunc
{
public:
  CovMarkov(const CovContext& ctx);
  CovMarkov(const CovMarkov &r);
  CovMarkov& operator= (const CovMarkov &r);
  virtual ~CovMarkov();

  virtual String getFormula() const override;
  String         getCovName() const override { return "Markov"; }

  bool   hasParam() const override { return true; }
  double getParMax() const override { return MAX_PARAM; }
  double getScadef() const override;
  bool   hasCovOnSphere() const override { return true; }
  bool   hasSpectrum() const override { return true; }
  bool   hasMarkovCoeffs() const override { return true; }

  double evaluateSpectrum(double freq, int ndim) const override;
  VectorDouble getMarkovCoeffs() const override {return _markovCoeffs;}
  void   setMarkovCoeffs(VectorDouble coeffs) override { _markovCoeffs = coeffs;}
  double getCorrec() const override { return _correc; }
  void   setCorrec(double val) override { _correc = val;}

protected:

  double _evaluateCov(double h) const override;
  double _evaluateCovOnSphere(double scale, int degree = 50) const override;

private :
  VectorDouble _markovCoeffs;
  double _correc;
};
