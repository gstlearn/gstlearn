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

VectorDouble CovMarkov::_evaluateSpectrumOnSphere(int n, double scale) const
{
  auto sp = _evaluateSpectrumOnSphereWithoutNormalization(n,scale);
  VH::normalize(sp,1);
  return sp;
}

VectorDouble CovMarkov::_evaluateSpectrumOnSphereWithoutNormalization(int n, double scale) const
{
  VectorDouble sp(1+n, 0.);

  for (int j = 0; j < (int)sp.size(); j++)
  {
    double nnp1 = scale * scale * (double) j * ((double) j + 1.);
    double s = 0.;
    for (int i = 0; i < (int)_markovCoeffs.size(); i++)
    {
      s += _markovCoeffs[i] * pow(nnp1,i);
    }
    sp[j] = (2. * j + 1.) / (4 * GV_PI * s);
  }
  return sp;

}

double CovMarkov::normalizeOnSphere(int n, double scale) const 
{ 
  auto sp = _evaluateSpectrumOnSphereWithoutNormalization(n,scale);
  double s = 0.;
  for (auto &e : sp)
  {
    s += e;
  }
  return s;
}


double CovMarkov::evaluateSpectrum(double freq) const
{
  double s = 0.;
  int n = (int)_markovCoeffs.size();
  if (n == 0)
  {
    return TEST;
  }
  for (int i = 0; i < n; i++)
  {
    s += _markovCoeffs[i] * pow(freq,i);
  }
  return 1. /  s;
}
