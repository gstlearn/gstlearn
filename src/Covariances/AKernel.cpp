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
#include "Covariances/AKernel.hpp"
#include "Basic/AException.hpp"
#include "Basic/FFT.hpp"
#include "Basic/Utilities.hpp"

#include <cmath>

namespace gstlrn
{
AKernel::AKernel(const ECov& type, const CovContext& ctxt)
  : AStringable()
  , _type(type)
  , _ctxt(ctxt)
  , _params()
{
  if (!isConsistent())
    my_throw("Cannot create such covariance function in that context");
}

AKernel::AKernel(const AKernel& r)
  : AStringable(r)
  , _type(r._type)
  , _ctxt(r._ctxt)
  , _params(r._params)
{
}

AKernel& AKernel::operator=(const AKernel& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _type   = r._type;
    _ctxt   = r._ctxt;
    _params = r._params;
  }
  return *this;
}

AKernel::~AKernel()
{
}

void AKernel::setParam(double param, Id ipar)
{
  /// TODO : Do not throw in setter. Check range and build the error message here.
  if (!hasParam()) return;
  double max = getParMax();
  if (param < 0. || (!FFFF(max) && param > max))
    my_throw("Wrong third parameter value");
  
  // Ensure the vector is large enough
  if (ipar >= static_cast<Id>(_params.size()))
    _params.resize(ipar + 1, TEST);
  
  _params[ipar] = param;
}

void AKernel::setParams(const VectorDouble& params)
{
  if (!hasParam()) return;
  double max = getParMax();
  
  // Validate all parameters
  for (const auto& param : params)
  {
    if (param < 0. || (!FFFF(max) && param > max))
      my_throw("Wrong parameter value");
  }
  
  _params = params;
}

double AKernel::getParam(Id ipar) const
{
  if (ipar >= static_cast<Id>(_params.size()))
    return TEST;
  return _params[ipar];
}

void AKernel::setField(double field)
{
  if (isZero(field))
    my_throw("Cannot scale with zero");
  _ctxt.setField(field);
}

double AKernel::evalCorFunc(double h) const
{
  return _evaluateCov(h);
}
double AKernel::evalCovDerivative(Id degree, double h) const
{
  return _evaluateCovDerivative(degree, h);
}
double AKernel::evalCovFirstDerivativeOverH(double h) const
{
  return _evaluateCovFirstDerivativeOverH(h);
}
double AKernel::evalCovOnSphere(double alpha,
                                 double scale,
                                 Id degree) const
{
  return _evaluateCovOnSphere(alpha, scale, degree);
}

VectorDouble AKernel::evalSpectrumOnSphere(Id n, double scale) const
{
  return _evaluateSpectrumOnSphere(n, scale);
}

String AKernel::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  sstr << getCovName();
  if (hasParam())
  {
    if (_params.size() == 1)
    {
      sstr << " (Third Parameter = " << getParam() << ")";
    }
    else if (_params.size() > 1)
    {
      sstr << " (Parameters =";
      for (Id i = 0; i < static_cast<Id>(_params.size()); i++)
      {
        sstr << " " << _params[i];
      }
      sstr << ")";
    }
  }
  sstr << std::endl;
  return sstr.str();
}

bool AKernel::hasCovOnSphere() const
{
  // If a spectrum is available, the covariance can be calculated
  return hasSpectrumOnSphere();
}

/// Test consistency with the current context
bool AKernel::isConsistent() const
{
  auto maxndim = getMaxNDim();
  if (maxndim <= 0.) return true;
  if (maxndim >= _ctxt.getNDim()) return true;
  /// TODO : Test irfDegree vs getMinOrder in CovElem because zonal anisotropies
  return false;
}

bool AKernel::hasInt1D() const
{
  return (getMaxNDim() >= 1 && getMinOrder() <= 0);
}
bool AKernel::hasInt2D() const
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
double AKernel::_evaluateCovDerivative(Id degree, double h) const
{
  DECLARE_UNUSED(degree);
  DECLARE_UNUSED(h);
  if (!hasCovDerivative())
  {
    messerr("This covariance does not allow Derivative calculations");
    return TEST;
  }
  messerr("This covariance should have Derivative calculations");
  messerr("But _evaluateCovDerivative has not been coded");
  my_throw("This should never happen");
}

void AKernel::setMarkovCoeffs(const VectorDouble& coeffs)
{
  DECLARE_UNUSED(coeffs);
  if (!hasMarkovCoeffs())
  {
    messerr("This covariance is not known to be Markovian");
  }
  messerr("This covariance should have a method giving the Markov coefficients");
  messerr("But getMarkovCoeffs has not been coded");
  my_throw("This should never happen");
}

VectorDouble AKernel::getMarkovCoeffs() const
{
  if (!hasMarkovCoeffs())
  {
    messerr("This covariance is not known to be Markovian");
    return VectorDouble();
  }
  messerr("This covariance should have a method giving the Markov coefficients");
  messerr("But getMarkovCoeffs has not been coded");
  my_throw("This should never happen");
}

double AKernel::evaluateSpectrum(double freq) const
{
  DECLARE_UNUSED(freq);
  if (!hasSpectrumOnRn())
  {
    messerr("This covariance does not allow spectrum calculations");
    return TEST;
  }
  messerr("This covariance should have a method giving the spectrum");
  messerr("But evaluateSpectrum has not been coded");
  my_throw("This should never happen");
}

double AKernel::_evaluateCovOnSphere(double alpha,
                                      double scale,
                                      Id degree) const
{
  double s = 0.;

  VectorDouble spectrum = _evaluateSpectrumOnSphere(degree, scale);

  if (isZero(alpha))
  {
    for (Id i = 0; i < degree; i++)
    {
      s += spectrum[i];
    }
  }
  else
  {
    double calpha = cos(alpha);
    double u0     = 1.;
    double u2     = 0.;
    double u1     = calpha;
    for (Id i = 1; i < (degree + 2); i++)
    {
      u2 = 1. / (i + 1) * ((2 * i + 1) * calpha * u1 - i * u0);
      s += u0 * spectrum[i - 1];
      u0 = u1;
      u1 = u2;
    }
  }
  return s;
}

VectorDouble AKernel::_evaluateSpectrumOnSphere(Id n, double scale) const
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

Array AKernel::_evalCovFFT(const VectorDouble& hmax, Id N) const
{
  N *= 2;
  Id ndim = static_cast<Id>(hmax.size());
  VectorInt nxs(ndim);
  for (Id idim = 0; idim < ndim; idim++)
    nxs[idim] = N;
  Array array(nxs);

  Id ntotal = static_cast<Id>(pow(N, ndim));
  VectorDouble a(ndim);
  double coeff = 0;
  double prod  = 1.;

  for (Id idim = 0; idim < ndim; idim++)
  {
    coeff   = 1. / (2. * hmax[idim]);
    a[idim] = GV_PI * (N - 1) / (hmax[idim]);
    prod *= coeff;
  }

  VectorDouble Re(ntotal, 0.);
  VectorDouble Im(ntotal, 0.);
  VectorInt indices(ndim);

  for (Id iad = 0; iad < ntotal; iad++)
  {
    array.rankToIndice(iad, indices);

    double s = 0.;
    for (Id idim = 0; idim < ndim; idim++)
    {
      double temp = a[idim] * (static_cast<double>(indices[idim]) / (N - 1) - 0.5);
      s += temp * temp;
    }
    Re[iad] = prod * evaluateSpectrum(sqrt(s));
    array.setValue(indices, Re[iad]);
  }

  FFTn(ndim, nxs, Re, Im);

  // Retrieve information from the Re array and load them back in the array result.

  VectorInt nxs2(ndim);
  for (Id idim = 0; idim < ndim; idim++)
    nxs2[idim] = N / 2;
  Array result(nxs2);
  VectorInt newIndices(ndim);

  for (Id iad = 0; iad < ntotal; iad++)
  {
    array.rankToIndice(iad, indices);
    for (Id idim = 0; idim < ndim; idim++)
    {
      Id odd           = indices[idim] % 2;
      Id s             = 1 - (2 * odd);
      newIndices[idim] = nxs[idim] / 2 + s * (indices[idim] / 2 + odd);
      Re[iad] *= s;
    }
    Re[iad] /= pow(2., ndim);
    array.setValue(newIndices, Re[iad]);
  }

  bool cont;
  Id iadr = 0;
  for (Id iad = 0; iad < ntotal; iad++)
  {
    array.rankToIndice(iad, indices);
    cont = true;
    for (Id idim = 0; idim < ndim; idim++)
    {
      if (indices[idim] < (nxs2[idim] / 2) || indices[idim] >= (3 * nxs2[idim] / 2))
      {
        cont = false;
        continue;
      }
    }
    if (cont)
    {
      result.rankToIndice(iadr++, newIndices);
      result.setValue(newIndices, array.getValue(indices));
    }
  }
  return result;
}

void AKernel::computeCorrec(Id ndim)
{
  if (!hasSpectrumOnRn()) return;
  Id N = static_cast<Id>(pow(2, 8));
  VectorInt Nv(ndim);
  VectorDouble hmax(ndim);
  for (Id idim = 0; idim < ndim; idim++)
  {
    hmax[idim] = 3 * getScadef();
    Nv[idim]   = N / 2;
  }
  Array res     = _evalCovFFT(hmax, N);
  double correc = res.getValue(Nv);
  setCorrec(correc);
}

} // namespace gstlrn
