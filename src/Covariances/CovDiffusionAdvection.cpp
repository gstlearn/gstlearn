#include "Covariances/CovDiffusionAdvection.hpp"
#include <cmath>
#include "Covariances/CovAniso.hpp"
#include "Basic/FFT.hpp"


CovDiffusionAdvection::CovDiffusionAdvection()
:   _markovL(nullptr)
  , _markovR(nullptr)
  , _scaleTime(1.)
  , _vel(VectorDouble())
  , _sigma2(1.)
  , _ndim(2)
  , _globalCorrec(1.)
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

    ut_vector_display("L",coeffsL);
    ut_vector_display("R",coeffsR);
    ut_vector_display("prod",coeffs);

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

  double normfreq = 0.;
  double velinner = 0.;

  for(int i = 0; i<(int) freq.size();i++)
  {
    velinner += _vel[i] * freq[i];
    normfreq += freq[i] * freq[i];
  }

  double s1 = 1.;
  double s2 = 1.;

  if (!isNoneMarkovL())
  {
    s1 = 1./(_markovL->evalSpectrum(normfreq));
  }

  if (!isNoneMarkovR())
  {
    s2 = 1./(_markovR->evalSpectrum(normfreq));
  }

  std::complex<double> temp = _scaleTime * (-1i * velinner * time);

  double ratio = 0.5 * _globalCorrec / (s1 * s2);
  return ratio * exp(temp) * exp(- _scaleTime * abs(time * s1));;

}

Array CovDiffusionAdvection::evalCovFFT(const VectorDouble& hmax,double time, int N) const
{
  N *= 2;
  VectorInt nxs(_ndim);
  for (int idim = 0; idim < _ndim; idim++)
    nxs[idim] = N ;
  Array array(nxs);


  int ntotal = pow(N, _ndim);
  VectorDouble a(_ndim);
  double coeff = 0;
  double prod = 1.;

  for(int idim = 0; idim < _ndim; idim++)
  {
    coeff = 1. / (2. * hmax[idim]);
    a[idim]=    GV_PI * (N-1) / ( hmax[idim]);
    prod *= coeff;
  }

  VectorDouble Re(ntotal);
  VectorDouble Im(ntotal);
  VectorInt indices(_ndim);
  VectorDouble temp(_ndim);

  for (int iad = 0; iad < ntotal; iad++)
  {
    array.rankToIndice(iad,indices);

    for (int idim = 0; idim < _ndim; idim++)
    {
      temp[idim] = a[idim] * ((double)indices[idim] / (N - 1) - 0.5);
    }

    std::complex<double> fourier = evalSpatialSpectrum(temp,time);
    Re[iad] = prod * fourier.real();
    Im[iad] = prod * fourier.imag();
    //array.setValue(indices,Re[iad]);
    //arrayImag.setValue(indices,Im[iad]);

  }

  FFTn(_ndim, nxs, Re, Im);


  // Retrieve information from the Re array and load them back in the array result.

  VectorInt nxs2(_ndim);
  for (int idim = 0; idim < _ndim; idim++)
    nxs2[idim] = N/2 ;
  Array result(nxs2);
  VectorInt newIndices(_ndim);

  for (int iad = 0; iad < ntotal; iad++)
  {
    array.rankToIndice(iad,indices);
    for (int idim = 0;  idim < _ndim; idim++)
    {
      int odd = indices[idim] % 2;
      int s = 1 - 2 * odd;
      newIndices[idim] = nxs[idim]/2 + s * (indices[idim]/2 + odd);
      Re[iad] *= s;

    }
    array.setValue(newIndices,Re[iad]);
  }

  bool cont;
  int iadr = 0;
  for (int iad = 0; iad < ntotal; iad++)
  {
    array.rankToIndice(iad,indices);
    cont = true;
    for (int idim = 0;  idim < _ndim; idim++)
    {
      if (indices[idim]<(nxs2[idim]/2) || indices[idim] >= (3 * nxs2[idim]/2 ))
      {
        cont = false;
        continue;
      }
    }
    if (cont)
    {
      result.rankToIndice(iadr++,newIndices);
      result.setValue(newIndices,  array.getValue(indices));
    }
  }

  return result;
}
