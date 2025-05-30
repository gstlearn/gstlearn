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
#include "Covariances/ACovFunc.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"
#include "Basic/FFT.hpp"

#include <math.h>

ACovFunc::ACovFunc(const ECov& type, const CovContext& ctxt)
: AStringable(),
  _type(type),
  _ctxt(ctxt),
  _param(TEST)
{
  if (!isConsistent())
    my_throw ("Cannot create such covariance function in that context");
}

ACovFunc::ACovFunc(const ACovFunc &r)
: AStringable(r),
  _type(r._type),
  _ctxt(r._ctxt),
  _param(r._param)
{
}

ACovFunc& ACovFunc::operator=(const ACovFunc &r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _type = r._type;
    _ctxt = r._ctxt;
    _param = r._param;
  }
  return *this;
}

ACovFunc::~ACovFunc()
{
}

void ACovFunc::setParam(double param)
{
  /// TODO : Do not throw in setter. Check range and build the error message here.
  if (! hasParam()) return;
  double max = getParMax();
  if (param < 0. || (!FFFF(max) && param > max))
    my_throw("Wrong third parameter value");
  _param = param;
}

void ACovFunc::setField(double field)
{
  if (isZero(field))
    my_throw("Cannot scale with zero");
  _ctxt.setField(field);
}

double ACovFunc::evalCorFunc(double h) const
{
  return _evaluateCov(h);
}
double ACovFunc::evalCovDerivative(int degree, double h) const
{
  return _evaluateCovDerivative(degree, h);
}

double ACovFunc::evalCovOnSphere(double alpha,
                                 double scale,
                                 int degree) const
{
  return _evaluateCovOnSphere(alpha, scale, degree);
}

VectorDouble ACovFunc::evalSpectrumOnSphere(int n, double scale) const
{
  return _evaluateSpectrumOnSphere(n, scale);
}

VectorDouble ACovFunc::evalCovVec(const VectorDouble& vech) const
{
  VectorDouble vec;
  for (const auto& h : vech)
    vec.push_back(evalCorFunc(h));
  return vec;
}
VectorDouble ACovFunc::evalCovDerivativeVec(int degree,
                                            const VectorDouble& vech) const
{
  VectorDouble vec;
  for (const auto& i : vech)
    vec.push_back(evalCovDerivative(degree, i));
  return vec;
}
String ACovFunc::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  sstr << getCovName();
  if (hasParam())
    sstr << " (Third Parameter = " << getParam() << ")";
  sstr << std::endl;
  return sstr.str();
}

bool ACovFunc::hasCovOnSphere() const
{
  // If a spectrum is available, the covariance can be calculated
  return hasSpectrumOnSphere();
}

/// Test consistency with the current context
bool ACovFunc::isConsistent() const
{
  unsigned int maxndim = getMaxNDim();
  if (maxndim <= 0.) return true;
  if (maxndim >= _ctxt.getNDim()) return true;
    /// TODO : Test irfDegree vs getMinOrder in CovElem because zonal anisotropies
  return false;
}

bool ACovFunc::hasInt1D() const
{
  return (getMaxNDim() >= 1 && getMinOrder() <= 0);
}
bool ACovFunc::hasInt2D() const
{
  return (getMaxNDim() >= 2 && getMinOrder() <= 0);
}
 /**
  * Calculate covariance derivatives, i.e.
  * - Degree 1: C^1(r) / r
  * - degree 2: C^2(r)
  * - Degree 3: C^3(r)
  * - Degree 4: C^4(r)
  * @param degree Level of derivation
  * @param h Normalized distance
  * @return
  */
double ACovFunc::_evaluateCovDerivative(int degree, double h) const
{
  DECLARE_UNUSED(degree);
  DECLARE_UNUSED(h);
  if (! hasCovDerivative())
  {
    messerr("This covariance does not allow Derivative calculations");
    return TEST;
  }
  messerr("This covariance should have Derivative calculations");
  messerr("But _evaluateCovDerivative has not been coded");
  my_throw("This should never happen");
}

void ACovFunc::setMarkovCoeffs(const VectorDouble& coeffs)
{
  DECLARE_UNUSED(coeffs);
  if (! hasMarkovCoeffs())
  {
    messerr("This covariance is not known to be Markovian");
  }
  messerr("This covariance should have a method giving the Markov coefficients");
  messerr("But getMarkovCoeffs has not been coded");
  my_throw("This should never happen");
}

VectorDouble ACovFunc::getMarkovCoeffs() const
{
  if (! hasMarkovCoeffs())
  {
      messerr("This covariance is not known to be Markovian");
      return VectorDouble();
  }
  messerr("This covariance should have a method giving the Markov coefficients");
  messerr("But getMarkovCoeffs has not been coded");
  my_throw("This should never happen");
}

double ACovFunc::evaluateSpectrum(double freq) const
{
  DECLARE_UNUSED(freq);
  if (! hasSpectrumOnRn())
  {
      messerr("This covariance does not allow spectrum calculations");
      return TEST;
  }
  messerr("This covariance should have a method giving the spectrum");
  messerr("But evaluateSpectrum has not been coded");
  my_throw("This should never happen");
}

double ACovFunc::_evaluateCovOnSphere(double alpha,
                                      double scale,
                                      int degree) const
{
  double s = 0.;

  VectorDouble spectrum = _evaluateSpectrumOnSphere(degree, scale);

  if (isZero(alpha))
  {
    for (int i = 0; i < degree; i++)
    {
      s += spectrum[i];
    }
  }
  else
  {
    double calpha = cos(alpha);
    double u0 = 1.;
    double u2 = 0.;
    double u1 = calpha;
    for (int i = 1; i < (degree + 2); i++)
    {
      u2 = 1. / (i + 1) * ((2 * i + 1) * calpha * u1 - i * u0);
      s += u0 * spectrum[i-1];
      u0 = u1;
      u1 = u2;
    }
  }
  return s;
}

VectorDouble ACovFunc::_evaluateSpectrumOnSphere(int n, double scale) const
{
  DECLARE_UNUSED(n);
  DECLARE_UNUSED(scale);
  if (!hasSpectrumOnSphere())
  {
    messerr("This covariance does not allow On Sphere calculations");
    return VectorDouble();
  }
  messerr("This covariance should have On Sphere calculations");
  messerr("But '_evaluateSpectrumOnSphere()' has not been coded");
  my_throw("This should never happen");
}

Array ACovFunc::_evalCovFFT(const VectorDouble& hmax, int N) const
{
  N *= 2;
  int ndim = (int) hmax.size();
  VectorInt nxs(ndim);
  for (int idim = 0; idim < ndim; idim++)
    nxs[idim] = N ;
  Array array(nxs);

  int ntotal = (int) pow(N, ndim);
  VectorDouble a(ndim);
  double coeff = 0;
  double prod = 1.;

  for(int idim = 0; idim < ndim; idim++)
  {
    coeff = 1. / (2. * hmax[idim]);
    a[idim]=    GV_PI * (N-1) / (hmax[idim]);
    prod *= coeff;
  }

  VectorDouble Re(ntotal,0.);
  VectorDouble Im(ntotal,0.);
  VectorInt    indices(ndim);

  for (int iad = 0; iad < ntotal; iad++)
  {
    array.rankToIndice(iad,indices);

    double s = 0.;
    for (int idim = 0; idim < ndim; idim++)
    {
      double temp = a[idim] * ((double) indices[idim] / (N - 1) - 0.5);
      s += temp * temp;
    }
    Re[iad] = prod * evaluateSpectrum(s);
    array.setValue(indices,Re[iad]);
  }

  FFTn(ndim, nxs, Re, Im);

  // Retrieve information from the Re array and load them back in the array result.

  VectorInt nxs2(ndim);
  for (int idim = 0; idim < ndim; idim++)
    nxs2[idim] = N/2 ;
  Array result(nxs2);
  VectorInt newIndices(ndim);

  for (int iad = 0; iad < ntotal; iad++)
  {
    array.rankToIndice(iad,indices);
    for (int idim = 0;  idim < ndim; idim++)
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
    for (int idim = 0;  idim < ndim; idim++)
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
      result.setValue(newIndices,array.getValue(indices));
    }
  }
  return result;
}

void ACovFunc::computeCorrec(int ndim) 
{
  if (! hasSpectrumOnRn()) return;
  int N = (int) pow(2,8);
  VectorInt Nv(ndim);
  VectorDouble hmax(ndim);
  for (int idim = 0; idim<ndim; idim++)
  {
    hmax[idim] = 3 * getScadef();
    Nv[idim] = N/2;
  }
  Array res = _evalCovFFT(hmax,N);
  double correc = res.getValue(Nv);
  setCorrec(correc);
}

