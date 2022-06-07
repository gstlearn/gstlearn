#include "Covariances/CovDiffusionAdvection.hpp"
#include "Covariances/CovAniso.hpp"
#include "Basic/FFT.hpp"

#include <cmath>

CovDiffusionAdvection::CovDiffusionAdvection()
    : _markovL(nullptr),
      _markovR(nullptr),
      _scaleTime(1.),
      _vel(VectorDouble()),
      _sigma2(1.),
      _ndim(2),
      _globalCorrec(1.)
{

}


CovDiffusionAdvection::CovDiffusionAdvection(const CovDiffusionAdvection& r)
:   _markovL((CovAniso*)r._markovL->clone())
  , _markovR((CovAniso*)r._markovR->clone())
  , _scaleTime(r._scaleTime)
  , _vel(r._vel)
  , _sigma2(r._sigma2)
  , _ndim(r._ndim)
  , _globalCorrec(r._globalCorrec)
{

}

CovDiffusionAdvection& CovDiffusionAdvection::operator=(const CovDiffusionAdvection& r)
{
  if (this != &r)
  {
    _markovL      = (CovAniso*)r._markovL->clone();
    _markovR      = (CovAniso*)r._markovR->clone();
    _scaleTime    = r._scaleTime;
    _vel          = r._vel;
    _sigma2       = r._sigma2;
    _ndim         = r._ndim;
    _globalCorrec = r._globalCorrec;
   }
   return *this;
}

CovDiffusionAdvection::~CovDiffusionAdvection()
{

}



CovDiffusionAdvection* CovDiffusionAdvection::create(CovAniso* markovL,
                                                     CovAniso* markovR ,
                                                     double scaleTime ,
                                                     VectorDouble vel,
                                                     double sigma2)
{
  CovDiffusionAdvection* cov = new CovDiffusionAdvection();
  cov->setMarkovL(markovL);
  cov->setMarkovR(markovR);
  cov->setScaleTime(scaleTime);
  cov->setVel(vel);
  cov->setSigma2(sigma2);

  if (!cov->isNoneMarkovL() && !cov->isNoneMarkovR())
  {
    CovContext ctxt = cov->getMarkovL()->getContext();
    CovAniso* covtemp = CovAniso::createIsotropic(ctxt,ECov::MARKOV,1.,1.,1);
    VectorDouble coeffsL = markovL->getMarkovCoeffs();
    VectorDouble coeffsR = markovR->getMarkovCoeffs();
    int degree = (coeffsL.size() + coeffsR.size() - 2);
    VectorDouble coeffs;
    coeffs.resize(degree + 1,0.);

    for(int i = 0; i<(int)coeffsL.size();i++)
      for(int j = 0; j<(int)coeffsL.size();j++)
      {
        coeffs[i+j] += coeffsL[i] * coeffsR[j];
      }

    covtemp->setMarkovCoeffs(coeffs);
    cov->_globalCorrec = markovL->getCorrec() * markovR->getCorrec() / covtemp->getCorrec();
    delete covtemp;
  }

  if (!cov->isNoneMarkovL())
  {
    cov->_ndim = markovL->getNDim();
  }
  else
  {
    cov->_ndim = markovR->getNDim();
  }


  return cov;
}

std::complex<double> CovDiffusionAdvection::evalSpatialSpectrum(VectorDouble freq, double time) const
{

  double velinner = 0.;

  for(int i = 0; i<(int) freq.size();i++)
  {
    velinner += _vel[i] * freq[i];
  }

  double s1 = 1.;
  if (!isNoneMarkovL())
  {
    s1 = 1./(_markovL->evalSpectrum(freq));
  }

  double s2 = 1.;
  if (!isNoneMarkovR())
  {
    s2 = 1./(_markovR->evalSpectrum(freq));
  }

  std::complex<double> temp = _scaleTime * (-1i * velinner * time - abs(time * s1));

  double ratio = 0.5 * _globalCorrec / (s1 * s2);
  return ratio * exp(temp);
}

Array CovDiffusionAdvection::evalCovFFT(const VectorDouble& hmax, double time, int N) const
{
  std::function<std::complex<double>(VectorDouble, double)> funcSpectrum;
  funcSpectrum = [this](VectorDouble freq, double time)
      { return evalSpatialSpectrum(freq, time);};

 return evalCovFFTTimeSlice(hmax, time, N, funcSpectrum);
}
