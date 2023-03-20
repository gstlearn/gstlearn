/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "Covariances/CovMarkov.hpp"
#include "Covariances/CovContext.hpp"
#include "Basic/MathFunc.hpp"

#include "math.h"

#define MAXTAB 100

CovMarkov::CovMarkov(const CovContext& ctxt)
: ACovFunc(ECov::MARKOV, ctxt)
  , _correc(1.)
{
  setParam(1);
  _markovCoeffs.push_back(1.);

}

CovMarkov::CovMarkov(const CovMarkov &r)
: ACovFunc(r)
{
  _markovCoeffs = r._markovCoeffs;
  _correc = r._correc;
}

CovMarkov& CovMarkov::operator=(const CovMarkov &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
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

double CovMarkov::_evaluateCov(double /*h*/) const
{
  return TEST;
}

String CovMarkov::getFormula() const
{
  return "C(h)=\\int_{R^d} \\frac{e^{-i\\omega^t.h}}{P(||\\omega||^2)}d\\omega";
}

double CovMarkov::_evaluateCovOnSphere(double scale, int degree) const
{
  double s = 0.;
  int n = (int)_markovCoeffs.size();
  double nnp1 = scale * scale * (double) degree * ((double) degree + 1.);
  for(int i = 0; i< n;i++)
  {
    s += _markovCoeffs[i] * pow(nnp1,i);
  }
  return scale * scale * (2. * degree + 1.) / (4 * GV_PI * s);
}

double CovMarkov::evaluateSpectrum(double freq, int /*ndim*/) const
{
  double s = 0.;
  int n = (int)_markovCoeffs.size();
  if (n == 0)
  {
    return TEST;
  }
  for(int i = 0; i< n;i++)
  {
      s += _markovCoeffs[i] * pow(freq,i);
  }
  return 1./  s;
}

