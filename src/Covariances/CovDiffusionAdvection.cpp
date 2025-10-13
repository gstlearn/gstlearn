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
#include "Covariances/CovDiffusionAdvection.hpp"
#include "Basic/AException.hpp"
#include "Basic/FFT.hpp"
#include "Covariances/CovAniso.hpp"
#include <cmath>
#include <complex>

namespace gstlrn
{
CovDiffusionAdvection::CovDiffusionAdvection()
  : _markovL(nullptr)
  , _markovR(nullptr)
  , _scaleTime(1.)
  , _vel(VectorDouble())
  , _sigma2(1.)
  , _globalCorrec(1.)
  , _spatialTrace(nullptr)
  , _ctxt(CovContext())
  , _destroyMarkovL(false)
  , _destroyMarkovR(false)
  , _markovRdefined(false)
  , _markovLdefined(false)
{
}

CovDiffusionAdvection::CovDiffusionAdvection(const CovDiffusionAdvection& r)
  : _markovL(r._markovL->clone())
  , _markovR(r._markovR->clone())
  , _scaleTime(r._scaleTime)
  , _vel(r._vel)
  , _sigma2(r._sigma2)
  , _globalCorrec(r._globalCorrec)
  , _spatialTrace(r._spatialTrace)
  , _destroyMarkovL(true) /// TODO : Why not using r flags?
  , _destroyMarkovR(true)
  , _markovRdefined(r._markovRdefined)
  , _markovLdefined(r._markovLdefined)
{
}

CovDiffusionAdvection& CovDiffusionAdvection::operator=(const CovDiffusionAdvection& r)
{
  if (this != &r)
  {
    _markovL        = r._markovL->clone();
    _markovR        = r._markovR->clone();
    _scaleTime      = r._scaleTime;
    _vel            = r._vel;
    _sigma2         = r._sigma2;
    _globalCorrec   = r._globalCorrec;
    _spatialTrace   = r._spatialTrace->clone();
    _destroyMarkovL = true; /// TODO : Why not using r flags?
    _destroyMarkovR = true;
    _markovRdefined = r._markovRdefined;
    _markovLdefined = r._markovLdefined;
  }
  return *this;
}

CovDiffusionAdvection::~CovDiffusionAdvection()
{
  delete _spatialTrace;
  if (_destroyMarkovL)
    delete _markovL;
  if (_destroyMarkovR)
    delete _markovR;
}

CovDiffusionAdvection* CovDiffusionAdvection::create(CovAniso* markovL,
                                                     CovAniso* markovR,
                                                     double scaleTime,
                                                     const VectorDouble& vel,
                                                     double sigma2)
{
  auto* cov = new CovDiffusionAdvection();
  cov->setMarkovL(markovL);
  cov->setMarkovR(markovR);
  cov->setScaleTime(scaleTime);
  cov->setVel(vel);
  cov->setSigma2(sigma2);
  cov->_init();

  return cov;
}

void CovDiffusionAdvection::_init()
{
  if (_markovL == nullptr && _markovR == nullptr)
  {
    my_throw("At least one of the covariances has to be defined to make a valid advection diffusion equation!");
  }

  const CovAniso* cova = _markovL == nullptr ? _markovR : _markovL;

  double correcR = 1.;
  double correcL = 1.;

  _ctxt   = cova->getContext();
  Id ndim = cova->getNDim();

  VectorDouble temp(ndim, 1.);
  if (_markovL == nullptr)
  {

    _markovL        = CovAniso::createAnisotropic(_ctxt, ECov::MARKOV, temp, 1., 1., temp, false);
    _destroyMarkovL = true;
    _markovLdefined = false;
    correcL         = 1.;
  }
  else
  {
    _markovLdefined = true;
    correcL         = _markovL->getCorrec();
  }
  if (_markovR == nullptr)
  {
    _markovR        = CovAniso::createAnisotropic(_ctxt, ECov::MARKOV, temp, 1., 1., temp, false);
    _destroyMarkovR = true;
    _markovRdefined = false;
    correcR         = 1.;
  }
  else
  {
    _markovRdefined = true;
    correcR         = _markovR->getCorrec();
  }

  _computeSpatialTrace();

  _globalCorrec = _spatialTrace->getFullCorrec() / (correcR * correcL);
}

void CovDiffusionAdvection::_computeSpatialTrace()
{
  delete _spatialTrace;

  VectorDouble scales;
  VectorDouble angles;
  if (_markovLdefined)
  {
    scales = _markovL->getScales();
    angles = _markovL->getAnisoAngles();
  }
  else
  {
    scales = _markovR->getScales();
    angles = _markovR->getAnisoAngles();
  }

  _spatialTrace = CovAniso::createAnisotropic(_ctxt, ECov::MARKOV, scales, _sigma2, 1., angles, false);

  VectorDouble coeffsL = _markovL->getMarkovCoeffs();
  VectorDouble coeffsR = _markovR->getMarkovCoeffs();
  Id degree            = (static_cast<Id>(coeffsL.size()) + static_cast<Id>(coeffsR.size()) - 2);
  VectorDouble coeffs;
  coeffs.resize(degree + 1, 0.);

  for (Id i = 0; i < static_cast<Id>(coeffsR.size()); i++)
    for (Id j = 0; j < static_cast<Id>(coeffsL.size()); j++)
    {
      coeffs[i + j] += coeffsR[i] * coeffsL[j];
    }

  _spatialTrace->setMarkovCoeffs(coeffs);
}

std::complex<double> CovDiffusionAdvection::evalSpatialSpectrum(VectorDouble freq, double time) const
{

  double velinner = 0.;

  for (Id i = 0; i < static_cast<Id>(freq.size()); i++)
  {
    velinner += _vel[i] * freq[i];
  }

  double s1 = 1.;

  if (_markovLdefined)
    s1 = 1. / (_markovL->evalSpectrum(freq));

  double s2 = 1.;

  if (_markovRdefined)
    s2 = 1. / (_markovR->evalSpectrum(freq));

  // std::complex<double> temp = _scaleTime * (-1i * velinner * time - abs(time * s1));
  std::complex<double> a(-_scaleTime * abs(time * s1), -_scaleTime * velinner * time);
  std::complex<double> temp = a;

  double ratio = _sigma2 / (_globalCorrec * s1 * s2);
  return ratio * exp(temp);
}

Array CovDiffusionAdvection::evalCovFFT(const VectorDouble& hmax, double time, Id N) const
{
  std::function<std::complex<double>(VectorDouble, double)> funcSpectrum;
  funcSpectrum = [this](VectorDouble freq, double time)
  { return evalSpatialSpectrum(freq, time); };

  return evalCovFFTTimeSlice(hmax, time, N, funcSpectrum);
}
} // namespace gstlrn