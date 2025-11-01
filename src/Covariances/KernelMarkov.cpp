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
#include "Covariances/KernelMarkov.hpp"
#include "Covariances/CovContext.hpp"

#include <cmath>

#define MAXTAB 100

namespace gstlrn
{
KernelMarkov::KernelMarkov(const CovContext& ctxt)
  : AKernel(ECov::MARKOV, ctxt)
  , _markovCoeffs()
  , _correc(1.)
{
  setParam(1);
  _markovCoeffs.push_back(1.);
}

KernelMarkov::KernelMarkov(const KernelMarkov& r)
  : AKernel(r)
  , _markovCoeffs(r._markovCoeffs)
  , _correc(r._correc)
{
}

KernelMarkov& KernelMarkov::operator=(const KernelMarkov& r)
{
  if (this != &r)
  {
    AKernel::operator=(r);
    _markovCoeffs = r._markovCoeffs;
    _correc       = r._correc;
  }
  return *this;
}

KernelMarkov::~KernelMarkov()
{
}

void KernelMarkov::setCorrec(double val)
{
  for(auto& e: _markovCoeffs)
  {
    e *= val;
  }
  _correc = 1;
}
double KernelMarkov::getScadef() const
{
  return sqrt(12. * _markovCoeffs.size());
}

String KernelMarkov::getFormula() const
{
  return "C(h)=\\int_{R^d} \\frac{e^{-i\\omega^t.h}}{P(||\\omega||^2)}d\\omega";
}

VectorDouble KernelMarkov::_evaluateSpectrumOnSphere(Id n, double scale) const
{
  auto sp = _evaluateSpectrumOnSphereWithoutNormalization(n, scale);
  VH::normalize(sp, 1);
  return sp;
}

VectorDouble KernelMarkov::_evaluateSpectrumOnSphereWithoutNormalization(Id n, double scale) const
{
  VectorDouble sp(1 + n, 0.);

  for (Id j = 0; j < static_cast<Id>(sp.size()); j++)
  {
    double nnp1 = scale * scale * static_cast<double>(j) * (static_cast<double>(j) + 1.);
    double s    = 0.;
    for (Id i = 0; i < static_cast<Id>(_markovCoeffs.size()); i++)
    {
      s += _markovCoeffs[i] * pow(nnp1, i);
    }
    sp[j] = (2. * j + 1.) / (4 * GV_PI * s);
  }
  return sp;
}

double KernelMarkov::normalizeOnSphere(Id n, double scale) const
{
  auto sp  = _evaluateSpectrumOnSphereWithoutNormalization(n, scale);
  double s = 0.;
  for (auto& e: sp)
  {
    s += e;
  }
  return s;
}

double KernelMarkov::evaluateSpectrum(double freq) const
{
  double s = 0.;
  size_t ndim = getContext().getNDim();

  Id n     = static_cast<Id>(_markovCoeffs.size());
  if (n == 0)
  {
    return TEST;
  }
  for (Id i = 0; i < n; i++)
  {
    s += _markovCoeffs[i] * pow(freq * freq, i);
  }
  return pow(2, ndim) / s;
}
} // namespace gstlrn