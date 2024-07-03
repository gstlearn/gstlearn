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
#include "Covariances/CovMarkov.hpp"
#include "Covariances/CovContext.hpp"
#include "Basic/MathFunc.hpp"

#include "math.h"

#define MAXTAB 100

CovMarkov::CovMarkov(const CovContext &ctxt)
    : ACovFunc(ECov::MARKOV, ctxt),
      _markovCoeffs(),
      _correc(1.)
{
  setParam(1);
  _markovCoeffs.push_back(1.);
}

CovMarkov::CovMarkov(const CovMarkov &r)
  : ACovFunc(r),
    _markovCoeffs(r._markovCoeffs),
    _correc(r._correc)
{
}

CovMarkov& CovMarkov::operator=(const CovMarkov &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
    _markovCoeffs = r._markovCoeffs;
    _correc = r._correc;
  }
  return *this;
}

CovMarkov::~CovMarkov()
{
}

double CovMarkov::getScadef() const
{
  return sqrt(12. * _markovCoeffs.size());
}

String CovMarkov::getFormula() const
{
  return "C(h)=\\int_{R^d} \\frac{e^{-i\\omega^t.h}}{P(||\\omega||^2)}d\\omega";
}

VectorDouble CovMarkov::_evaluateSpectrumOnSphere(double scale) const
{
  int degree = getDegree();
  int nm = (int) _markovCoeffs.size();
  VectorDouble sp(1+degree, 0.);

  for (int k = 0; k <= degree; k++)
  {
    double nnp1 = scale * scale * (double) k * ((double) k + 1.);

    double s = 0.;
    for (int i = 0; i < nm; i++)
      s += _markovCoeffs[i] * pow(nnp1,i);
    sp[k] = scale * scale * (2. * k + 1.) / (4. * GV_PI * s);
  }

  if (isFlagNormalizeSpectrum()) VH::normalize(sp, 1);

  return sp;
}

double CovMarkov::evaluateSpectrum(double freq, int /*ndim*/) const
{
  double s = 0.;
  int nm = (int)_markovCoeffs.size();
  if (nm == 0) return TEST;
  for (int i = 0; i < nm; i++)
    s += _markovCoeffs[i] * pow(freq,i);
  return 1. /  s;
}
