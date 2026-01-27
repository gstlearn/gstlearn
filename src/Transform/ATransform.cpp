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

#include "Transform/ATransform.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Law.hpp"
#include "Basic/VectorNumT.hpp"
#include "geoslib_define.h"

#include <cmath>
namespace gstlrn
{

double ATransform::inverseTransform(double y) const
{
  double x           = y; // initialisation naturelle (si la fonction est proche de l’identité)
  const int max_iter = 2000;
  const double tol   = 1e-5;

  for (int i = 0; i < max_iter; ++i)
  {
    double fx  = transform(x) - y;
    double dfx = evalJacobian(x);
    double dx  = fx / dfx;

    x -= dx;
    if (std::fabs(dx) < tol)
      return x;
  }
  messerr("ATransform::inverse: no convergence");
  return x;
}

String ATransform::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  sstr << "Type: " << getName() << "\n";
  if (hasParameters())
  {
    _printParams(sstr, strfmt);
  }
  return sstr.str();
}
double ATransform::condExp(double mu, double sigma, Id power) const
{
  Id oldseed = law_get_random_seed();
  // Monte Carlo estimation of E[transform(X)] where X ~ N(mu, sigma^2)
  double sum = 0.0;

  for (Id i = 0; i < _nMonteCarlo; ++i)
  {
    double x = mu + (sigma * law_gaussian());
    sum += pow(transform(x), power);
  }
  law_set_random_seed(oldseed);
  return sum / _nMonteCarlo;

}

void ATransform::condExpVec(constvect mu, constvect sigma, vect out, Id power) const
{
  Id n = static_cast<Id>(mu.size());
  for (Id i = 0; i < n; ++i)
  {
    out[i] = condExp(mu[i], sigma[i], power);
  }
}

VectorDouble ATransform::condExpVec(const VectorDouble& mu, const VectorDouble& sigma, Id power) const
{
  Id n = static_cast<Id>(mu.size());
  VectorDouble out(n);
  condExpVec(mu, sigma, out, power);
  return out;
}

VectorDouble ATransform::transformVec(const VectorDouble& in) const
{
  Id n = static_cast<Id>(in.size());
  VectorDouble out(n);
  transformVec(in, out);
  return out;
}

VectorDouble ATransform::inverseTransformVec(const VectorDouble& in) const
{
  Id n = static_cast<Id>(in.size());
  VectorDouble out(n);
  inverseTransformVec(in, out);
  return out;
}
double ATransform::evalJacobian(double x) const
{
  double eps = EPSILON4;
  return (transform(x + eps) - transform(x - eps)) / (2. * eps);
}

void ATransform::transformVec(constvect in, vect out) const
{
  Id n = static_cast<Id>(in.size());
  for (Id i = 0; i < n; i++)
    out[i] = transform(in[i]);
}

void ATransform::inverseTransformVec(constvect in, vect out) const
{
  Id n = static_cast<Id>(in.size());
  for (Id i = 0; i < n; i++)
    out[i] = inverseTransform(in[i]);
}

double ATransform::evalLogJacobianVec(constvect in) const
{
  double logJacobian = 0.0;
  Id n               = in.size();
  for (Id i = 0; i < n; i++)
    logJacobian += std::log(std::fabs(evalJacobian(in[i])));
  return logJacobian;
}

} // namespace gstlrn