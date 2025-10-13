/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Covariances/ACovFunc.hpp"

namespace gstlrn
{
  
class CovContext;

class GSTLEARN_EXPORT CovMarkov : public ACovFunc
{
public:
  CovMarkov(const CovContext& ctx);
  CovMarkov(const CovMarkov &r);
  CovMarkov& operator= (const CovMarkov &r);
  virtual ~CovMarkov();

  String getFormula() const override;
  String         getCovName() const override { return "Markov"; }

  bool   hasParam() const override { return true; }
  double getParMax() const override { return MAX_PARAM; }
  double getScadef() const override;
  bool   hasCovOnRn() const override { return false; }
  bool   hasSpectrumOnRn() const override { return true; }
  bool   hasSpectrumOnSphere() const override { return true; }
  bool   hasMarkovCoeffs() const override { return true; }
  double normalizeOnSphere(Id n = 50,double scale = 1.) const override;
  double evaluateSpectrum(double freq) const override;
  VectorDouble getMarkovCoeffs() const override { return _markovCoeffs; }
  void   setMarkovCoeffs(const VectorDouble& coeffs) override { _markovCoeffs = coeffs;}
  double getCorrec() const override { return _correc; }
  void   setCorrec(double val) override { _correc = val;}

protected:
  VectorDouble _evaluateSpectrumOnSphere(Id n,
                                         double scale = 1.) const override;

private:
  VectorDouble _evaluateSpectrumOnSphereWithoutNormalization(Id n, double scale = 1.) const;

private :
  VectorDouble _markovCoeffs;
  double _correc;
};
}