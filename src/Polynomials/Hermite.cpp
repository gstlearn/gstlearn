/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "Polynomials/Hermite.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "geoslib_f.h"
#include "geoslib_old_f.h"

#include <math.h>

double _convert2u(double yc, double krigest, double krigstd)
{
  if (ABS(krigstd) < EPSILON6)
  {
    return ((yc >= krigest) ? +10. :
                              -10.);
  }
  else
  {
    return ((yc - krigest) / krigstd);
  }
}

/**
 * Calculate In(u,v) = int_u^v H_n(yk + sk * t) g(t) dt (with v = +Inf)
 * @param In Returned vector (in place)
 * @param yk Kriged value
 * @param sk Standard deviation of estimation error
 * @param u  Lower bound of integral (-inf if set to TEST)
 * @param hnYc Vector of Hermite polynomials at cutoff
 */
void _calculateIn(VectorDouble &In,
                  double yk,
                  double sk,
                  double u,
                  const VectorDouble &hnYc)
{
  double Gcomp, gusk;
  bool flag_u = !FFFF(u);
  double r2 = 1 - sk * sk;
  int nbpoly = static_cast<int>(In.size());
  if (flag_u)
  {
    Gcomp = 1 - law_cdf_gaussian(u);
    gusk  = sk * law_df_gaussian(u);
  }
  else
  {
    Gcomp = 1.;
    gusk  = 0.;
  }

  In.resize(nbpoly, 0.);
  In[0] = Gcomp;
  In[1] = -(yk * In[0] + gusk);

  double cutval = 0.;
  for (int ih = 2; ih < nbpoly; ih++)
  {
    if (flag_u) cutval = gusk * hnYc[ih - 1];
    double sqh  = sqrt((double) ih);
    double sqh1 = sqrt((double) (ih - 1.));
    In[ih] = -(yk * In[ih - 1] + r2 * sqh1 * In[ih - 2] + cutval) / sqh;
  }
}

/**
 * Calculate matrix JJ(n,p) = int_u_v H_n(y+st) H_p(y+st) g(t) dt (with v = +Inf)
 * @param JJ Output matrix
 * @param In Preliminary calculations (see _calculateII)
 * @param yk Kriged value
 * @param sk Standard deviation of Krging error
 * @param u  Lower bound for integration
 * @param hnYc Vector of Hermite polynomials at cutoff
 * @param phi
 */
void _calculateJJ(MatrixSquareGeneral &JJ,
                  VectorDouble &In,
                  double yk,
                  double sk,
                  double u,
                  const VectorDouble &hnYc,
                  const VectorDouble &phi)
{
  int nbpoly = static_cast<int>(phi.size());

  bool flag_u = !FFFF(u);
  double s2 = sk * sk;
  double r2 = 1 - s2;
  double gusk = (flag_u) ? sk * law_df_gaussian(u) :  0.;

  _calculateIn(In, yk, sk, u, hnYc);

  double cutval = 0.;
  for (int n = 0; n < nbpoly; n++)
  {
    JJ.setValue(n, 0, In[n]);
    JJ.setValue(0, n, In[n]);
  }
  for (int n = 1; n < nbpoly; n++)
  {
    if (flag_u) cutval = gusk * hnYc[n];
    double sqn = sqrt((double) n);
    double value = -yk * JJ.getValue(n, 0) + s2 * sqn * JJ.getValue(n - 1, 0) - cutval;

    JJ.setValue(n, 1, value);
    JJ.setValue(1, n, value);
  }
  for (int n = 1; n < nbpoly; n++)
  {
    double sqn  = sqrt((double) n);
    double sqn1 = sqrt((double) (n + 1.));
    for (int p = n + 1; p < nbpoly; p++)
    {
      if (flag_u) cutval = gusk * hnYc[n] * hnYc[p];
      double sqp = sqrt((double) p);
      double value = -(yk * JJ.getValue(n, p) + sqn * r2 * JJ.getValue(n - 1, p)
          - sqp * s2* JJ.getValue(n, p - 1) + cutval) / sqn1;
      JJ.setValue(n + 1, p, value);
      JJ.setValue(p, n + 1, value);
    }
  }
}

/**
 * Calculation of the Hermite Polynomials for a given value
 *
 * @param y Gaussian value for which the Hermite polynomials are calculated
 * @param r Change of support coefficient
 * @param nbpoly Number of Hermite polynomials
 * @return The vector of polynomials (Dimension: nbpoly)
 */
VectorDouble hermitePolynomials(double y, double r, int nbpoly)
{
  VectorDouble poly(nbpoly);
  if (nbpoly < 1) return poly;

  poly[0] = 1.;
  if (nbpoly > 1)
  {
    poly[1] = -y;
    if (nbpoly > 2)
    {
      for (int ih = 2; ih < nbpoly; ih++)
      {
        double sqh   = sqrt((double) ih);
        double sqhm1 = sqrt((double) (ih - 1.));
        poly[ih] = -(y * poly[ih - 1] + sqhm1 * poly[ih - 2]) / sqh;
      }
    }
  }

  if (r != 1)
  {
    double rk = 1;
    for (int ih = 0; ih < nbpoly; ih++)
    {
      poly[ih] *= rk;
      rk *= r;
    }
  }
  return poly;
}

/**
 * Returns the vector of Hermite Polynomials selected by ranks
 * @param y Target variable
 * @param r Change of support coefficient
 * @param ifacs Vector of ranks (staring from 0)
 * @return The vector of Hi(y) where 'i' is in 'ifacs'
 */
VectorDouble hermitePolynomials(double y, double r, const VectorInt& ifacs)
{
  int nfact = (int) ifacs.size();
  VectorDouble vec(nfact);

  int nbpoly = ut_ivector_max(ifacs);
  VectorDouble poly = hermitePolynomials(y, r, nbpoly);

  for (int ifac = 0; ifac < nfact; ifac++)
    vec[ifac] = poly[ifacs[ifac]];

  return poly;
}

/**
 * Calculate the Conditional Expectation:
 *    E[Z | Z1=z1, Z2=z2, ..., Zn=zn] = int Phi(y_kk + s_k u) g(u) du
 *
 * @param krigest Vector of Kriging estimates
 * @param krigstd Vector of Kriging standard deviations
 * @param phi Array of Hermite coefficients
 * @return Conditional Expectation
 */
VectorDouble hermiteCondExp(VectorDouble krigest,
                            VectorDouble krigstd,
                            const VectorDouble &phi)
{
  VectorDouble condexp;

  int nech = static_cast<int>(krigest.size());
  condexp.resize(nech);

  for (int iech = 0; iech < nech; iech++)
  {
    condexp[iech] = hermiteCondExpElement(krigest[iech], krigstd[iech], phi);
  }
  return condexp;
}

double hermiteCondExpElement(double krigest,
                             double krigstd,
                             const VectorDouble &phi)
{
  int nbpoly = static_cast<int>(phi.size());
  VectorDouble In(nbpoly);
  _calculateIn(In, krigest, krigstd, TEST, VectorDouble());

  double condexp = 0.;
  for (int ih = 0; ih < nbpoly; ih++)
    condexp += phi[ih] * In[ih];
  return condexp;
}

/**
 * Vector of conditional variances (same dimension as krigest and krigstd)
 * @param krigest Vector of Kriging estimate
 * @param krigstd Vector of Kriging standard deviations
 * @param phi Array of Hermite coefficients
 * @return
 */
VectorDouble hermiteCondStd(VectorDouble krigest,
                            VectorDouble krigstd,
                            const VectorDouble &phi)
{
  int nech = static_cast<int>(krigest.size());
  VectorDouble condstd(nech, 0);

  /* Loop on the samples */

  for (int iech = 0; iech < nech; iech++)
    condstd[iech] = hermiteCondStdElement(krigest[iech], krigstd[iech], phi);

  return condstd;
}

double hermiteCondStdElement(double krigest,
                             double krigstd,
                             const VectorDouble &phi)
{
  MatrixSquareGeneral JJ;
  int nbpoly = static_cast<int>(phi.size());
  VectorDouble In(nbpoly);
  JJ.reset(nbpoly, nbpoly, TEST);
  _calculateJJ(JJ, In, krigest, krigstd, TEST, VectorDouble(), phi);

  double constd = 0.;
  for (int ih = 0; ih < nbpoly; ih++)
    for (int jh = 0; jh < nbpoly; jh++)
      constd += JJ.getValue(ih, jh) * phi[ih] * phi[jh];

  double condexp = hermiteCondExpElement(krigest, krigstd, phi);
  constd -= condexp * condexp;
  constd = (constd > 0) ? sqrt(constd) : 0.;

  return constd;
}

/**
 *
 * @param yc Cutoff Value
 * @param krigest Estimation
 * @param krigstd Standard deviation of estimation error
 * @return The indicator above Cutoff
 */
VectorDouble hermiteIndicator(double yc,
                              VectorDouble krigest,
                              VectorDouble krigstd)
{
  int nech = static_cast<int>(krigest.size());
  VectorDouble proba(nech);

  for (int iech = 0; iech < nech; iech++)
  {
    proba[iech] = hermiteIndicatorElement(yc, krigest[iech], krigstd[iech]);
  }
  return proba;
}

double hermiteIndicatorElement(double yc, double krigest, double krigstd)
{
  double proba;

  double std = krigstd;
  if (ABS(std) < EPSILON6) std = EPSILON6;
  proba = 1. - law_cdf_gaussian((yc - krigest) / std);
  return proba;
}

VectorDouble hermiteIndicatorStd(double yc,
                                 VectorDouble krigest,
                                 VectorDouble krigstd)
{
  int nech = static_cast<int>(krigest.size());
  VectorDouble probstd(nech);

  for (int iech = 0; iech < nech; iech++)
    probstd[iech] = hermiteIndicatorStdElement(yc, krigest[iech],
                                               krigstd[iech]);

  return probstd;
}

double hermiteIndicatorStdElement(double yc, double krigest, double krigstd)
{
  double proba = hermiteIndicatorElement(yc, krigest, krigstd);
  double probstd = sqrt(proba * (1. - proba));
  return probstd;
}

/**
 *
 * @param yc Cutoff Value
 * @param krigest Estimation
 * @param krigstd Standard deviation of estimation error
 * @param phi  Hermite coefficients
 * @return The Metal
 */
VectorDouble hermiteMetal(double yc,
                          VectorDouble krigest,
                          VectorDouble krigstd,
                          const VectorDouble &phi)
{
  int nech = static_cast<int>(krigest.size());
  int nbpoly = static_cast<int>(phi.size());
  VectorDouble In(nbpoly);
  VectorDouble metal(nech);
  VectorDouble hnYc = hermitePolynomials(yc, 1., nbpoly);

  for (int iech = 0; iech < nech; iech++)
  {
    double u = _convert2u(yc, krigest[iech], krigstd[iech]);
    _calculateIn(In, krigest[iech], krigstd[iech], u, hnYc);

    double result = 0.;
    for (int ih = 0; ih < nbpoly; ih++)
      result += phi[ih] * In[ih];
    metal[iech] = result;
  }
  return metal;
}

double hermiteMetalElement(double yc,
                           double krigest,
                           double krigstd,
                           const VectorDouble &phi)
{
  int nbpoly = static_cast<int>(phi.size());
  VectorDouble In(nbpoly);
  VectorDouble hnYc = hermitePolynomials(yc, 1., nbpoly);

  double u = _convert2u(yc, krigest, krigstd);
  _calculateIn(In, krigest, krigstd, u, hnYc);

  double result = 0.;
  for (int ih = 0; ih < nbpoly; ih++)
    result += phi[ih] * In[ih];

  return result;
}

VectorDouble hermiteMetalStd(double yc,
                             VectorDouble krigest,
                             VectorDouble krigstd,
                             const VectorDouble &phi)
{
  MatrixSquareGeneral JJ;

  int nech = static_cast<int>(krigest.size());
  int nbpoly = static_cast<int>(phi.size());
  VectorDouble In(nbpoly);
  JJ.reset(nbpoly, nbpoly, TEST);

  VectorDouble metstd(nech, 0.);
  VectorDouble hnYc = hermitePolynomials(yc, 1., nbpoly);
  VectorDouble metal = hermiteMetal(yc, krigest, krigstd, phi);

  for (int iech = 0; iech < nech; iech++)
  {
    double u = _convert2u(yc, krigest[iech], krigstd[iech]);
    _calculateJJ(JJ, In, krigest[iech], krigstd[iech], u, hnYc, phi);

    double result = 0.;
    for (int ih = 0; ih < nbpoly; ih++)
      for (int jh = 0; jh < nbpoly; jh++)
        result += JJ.getValue(ih, jh) * phi[ih] * phi[jh];
    result -= metal[iech] * metal[iech];
    if (result > 0) metstd[iech] = sqrt(result);
  }
  return metstd;
}

double hermiteMetalStdElement(double yc,
                              double krigest,
                              double krigstd,
                              const VectorDouble &phi)
{
  MatrixSquareGeneral JJ;
  int nbpoly = static_cast<int>(phi.size());
  VectorDouble In(nbpoly);
  JJ.reset(nbpoly, nbpoly, TEST);
  VectorDouble hnYc = hermitePolynomials(yc, 1., nbpoly);

  double u = _convert2u(yc, krigest, krigstd);
  _calculateJJ(JJ, In, krigest, krigstd, u, hnYc, phi);

  double result = 0.;
  for (int ih = 0; ih < nbpoly; ih++)
    for (int jh = 0; jh < nbpoly; jh++)
      result += JJ.getValue(ih, jh) * phi[ih] * phi[jh];

  double metal = hermiteMetalElement(yc, krigest, krigstd, phi);
  result -= metal * metal;

  double metstd = (metal > 0) ? sqrt(result) : 0.;

  return metstd;
}

/**
 *
 * @param yc Cutoff Value
 * @param nbpoly Number of Hermite polynomials
 * @return The vector of coefficients of the Indicator
 */
VectorDouble hermiteCoefIndicator(double yc, int nbpoly)
{
  VectorDouble hn = hermitePolynomials(yc, 1., nbpoly);
  VectorDouble an(nbpoly);
  double gyc = law_df_gaussian(yc);
  an[0] = 1. - law_cdf_gaussian(yc);
  for (int n = 1; n < nbpoly; n++)
    an[n] = -gyc * hn[n - 1] / sqrt(n);

  return an;
}

/**
 *
 * @param yc Cutoff Value
 * @param phi Coefficients of Hermite polynomial
 * @return The vector of coefficients of the Metal Quantity
 */
VectorDouble hermiteCoefMetal(double yc, const VectorDouble &phi)
{
  int nbpoly = static_cast<int>(phi.size());
  VectorDouble vect(nbpoly);
  MatrixSquareGeneral TAU = hermiteIncompleteIntegral(yc, nbpoly);
  TAU.prodVector(phi, vect);
  return vect;
}

/**
 *
 * @param yc Cutoff Value
 * @param nbpoly Number of Hermite polynomials
 * @return The matrix of Incomplete Integral (Dimension: nbpoly * nbpoly)
 */
MatrixSquareGeneral hermiteIncompleteIntegral(double yc, int nbpoly)
{
  MatrixSquareGeneral TAU;

  TAU.reset(nbpoly, nbpoly, 0.);
  VectorDouble hn = hermitePolynomials(yc, 1., nbpoly);
  double gy = law_df_gaussian(yc);

  /* Calculation of S_0n */

  TAU.setValue(0, 0, law_cdf_gaussian(yc));
  for (int ip = 1; ip < nbpoly; ip++)
  {
    double aa = hn[ip - 1] / sqrt(ip) * gy;
    TAU.setValue(ip, 0, aa);
    TAU.setValue(0, ip, aa);
  }

  /* Calculation of diagonals and symmetrization */

  for (int n = 0; n < nbpoly - 1; n++)
    for (int m = 1; m < nbpoly - n; m++)
    {
      double aa = sqrt((double) m / (double) (m + n))
          * TAU.getValue(m - 1, n + m - 1)
                  + gy * hn[m] * hn[m + n - 1] / sqrt((double) (m + n));
      TAU.setValue(m, m + n, aa);
      TAU.setValue(m + n, m, aa);
    }

  for (int n = 0; n < nbpoly; n++)
    for (int m = 0; m < nbpoly; m++)
      TAU.setValue(m, n, (m == n) ? 1. - TAU.getValue(m, n) : -TAU.getValue(m, n));

  return TAU;
}

/**
 * Hermite coefficient for a lognormal transform
 *        mean * exp(sigma * Y + 1/2 * sigma^2)
 * @param mean Mean value
 * @param sigma Standard deviation
 * @param nbpoly Number of Hermite polynomials
 * @return The array of coefficients
 */
VectorDouble hermiteLognormal(double mean, double sigma, int nbpoly)
{
  VectorDouble hn(nbpoly);

  double fact = 1.;
  hn[0] = mean;
  for (int i = 1; i < nbpoly; i++)
  {
    fact *= (double) i;
    hn[i] = mean * pow(-sigma, (double) i) / sqrt(fact);
  }
  return hn;
}

/**
 * Evaluate the Hermite expansion
 * @param an Series of coefficients of the Hermite polynomials
 * @param hn Hermite polynomial values
 * @return The result of the expansion
 */
double hermiteSeries(const VectorDouble &an, const VectorDouble &hn)
{
  double value = 0.;
  for (int ih = 0; ih < (int) hn.size(); ih++)
  {
    value += an[ih] * hn[ih];
  }
  return value;
}

/**
 * Returns the vector of Hermite coefficients of the gaussian floored at 'y'
 * @param y Floor value
 * @param nbpoly Number of Polynomial functions
 * @return Hermite Coefficients
 */
VectorDouble hermiteCoefLower(double y, int nbpoly)
{
  VectorDouble hn = hermitePolynomials(y, 1., nbpoly);
  VectorDouble coeff(nbpoly);

  /* Calculate the coefficients */

  double dg = law_df_gaussian(y);
  double pg = law_cdf_gaussian(y);
  coeff[0] = dg + y * pg;
  coeff[1] = pg - 1.;
  for (int n = 2; n < nbpoly; n++)
  {
    double sqn = sqrt((double) n);
    double sqnnm1 = sqrt((double) n * (n - 1.));
    coeff[n] = dg * (y * hn[n - 1] / sqn + hn[n - 2] / sqnnm1);
  }
  return coeff;
}
