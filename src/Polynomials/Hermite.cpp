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

void _calculateIn(VectorDouble &In,
                  double yk,
                  double sk,
                  double u,
                  const VectorDouble &hnYc)
{
  double Gcomp, gusk;
  bool flag_u = !FFFF(u);
  double s2 = 1 - sk * sk;
  int nbpoly = static_cast<int>(In.size());
  if (flag_u)
  {
    Gcomp = 1 - law_cdf_gaussian(u);
    gusk = sk * law_df_gaussian(u);
  }
  else
  {
    Gcomp = 1.;
    gusk = 0.;
  }

  In.resize(nbpoly, 0.);
  In[0] = Gcomp;
  In[1] = -(yk * In[0] + gusk);

  double cutval = 0.;
  for (int ih = 2; ih < nbpoly; ih++)
  {
    if (flag_u) cutval = gusk * hnYc[ih - 1];
    In[ih] = -(yk * In[ih - 1] + s2 * sqrt((double) ih - 1) * In[ih - 2]
               + cutval)
             / sqrt((double) ih);
  }
}

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
  double gusk = (flag_u) ? sk * law_df_gaussian(u) :
                           0.;

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
    double sq2n = s2 * sqrt((double) n);
    double value = -yk * JJ.getValue(n, 0) + sq2n * JJ.getValue(n - 1, 0)
        - cutval;

    JJ.setValue(n, 1, value);
    JJ.setValue(1, n, value);
  }
  for (int n = 1; n < nbpoly; n++)
  {
    double sq2n = (1. - s2) * sqrt((double) n);
    double sqn1 = sqrt((double) (n + 1.));
    for (int p = n + 1; p < nbpoly; p++)
    {
      if (flag_u) cutval = gusk * hnYc[n] * hnYc[p];
      double sq2p = s2 * sqrt((double) p);
      double value = -(yk * JJ.getValue(n, p) + sq2n * JJ.getValue(n - 1, p)
          - sq2p * JJ.getValue(n, p - 1)
                       + cutval)
                     / sqn1;
      JJ.setValue(n + 1, p, value);
      JJ.setValue(p, n + 1, value);
    }
  }
}

/**
 *
 * @param y Gaussian value for which the Hermite polynomials are calculated
 * @param r Change of support coefficient
 * @param nbpoly Number of Hermite polynomials
 * @return The vector of polynomials (Dimension: nbpoly)
 */
GSTLEARN_EXPORT VectorDouble hermitePolynomials(double y, double r, int nbpoly)
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
        poly[ih] = -(y * poly[ih - 1] + sqrt((double) (ih - 1)) * poly[ih - 2])
            / sqrt((double) ih);
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
 * Calculate the Conditional Expectation:
 *    E[Z | Z1=z1, Z2=z2, ..., Zn=zn] = int Phi(y_kk + s_k u) g(u) du
 *
 * @param krigest Vector of Kriging estimated
 * @param krigstd Vector of Kriging standard deviations
 * @param phi Array of Hermite coefficients
 * @return Conditional Expectation
 */
GSTLEARN_EXPORT VectorDouble hermiteCondExp(VectorDouble krigest,
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

GSTLEARN_EXPORT double hermiteCondExpElement(double krigest,
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
GSTLEARN_EXPORT VectorDouble hermiteCondStd(VectorDouble krigest,
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

GSTLEARN_EXPORT double hermiteCondStdElement(double krigest,
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
  constd = (constd > 0) ? sqrt(constd) :
                          0.;

  return constd;
}

/**
 *
 * @param yc Cutoff Value
 * @param krigest Estimation
 * @param krigstd Standard deviation of estimation error
 * @return The indicator above Cutoff
 */
GSTLEARN_EXPORT VectorDouble hermiteIndicator(double yc,
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

GSTLEARN_EXPORT double hermiteIndicatorElement(double yc,
                                               double krigest,
                                               double krigstd)
{
  double proba;

  double std = krigstd;
  if (ABS(std) < EPSILON6) std = EPSILON6;
  proba = 1. - law_cdf_gaussian((yc - krigest) / std);
  return proba;
}

GSTLEARN_EXPORT VectorDouble hermiteIndicatorStd(double yc,
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

GSTLEARN_EXPORT double hermiteIndicatorStdElement(double yc,
                                                  double krigest,
                                                  double krigstd)
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
GSTLEARN_EXPORT VectorDouble hermiteMetal(double yc,
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

GSTLEARN_EXPORT double hermiteMetalElement(double yc,
                                           double krigest,
                                           double krigstd,
                                           const VectorDouble &phi)
{
  int nbpoly = static_cast<int>(phi.size());
  VectorDouble In(nbpoly);
  VectorDouble hnYc = hermitePolynomials(yc, 1., nbpoly);

  double u = _convert2u(yc, krigest, krigstd);
  _calculateIn(In, krigest, krigstd, u, hnYc);

  double metal = 0.;
  for (int ih = 0; ih < nbpoly; ih++)
    metal += phi[ih] * In[ih];

  return metal;
}

GSTLEARN_EXPORT VectorDouble hermiteMetalStd(double yc,
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

GSTLEARN_EXPORT double hermiteMetalStdElement(double yc,
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

  double metstd = (metal > 0) ? sqrt(result) :
                                0.;

  return metstd;
}

/**
 *
 * @param yc Cutoff Value
 * @param nbpoly Number of Hermite polynomials
 * @return The vector of coefficients of the Indicator
 */
GSTLEARN_EXPORT VectorDouble hermiteCoefIndicator(double yc, int nbpoly)
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
GSTLEARN_EXPORT VectorDouble hermiteCoefMetal(double yc,
                                              const VectorDouble &phi)
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
GSTLEARN_EXPORT MatrixSquareGeneral hermiteIncompleteIntegral(double yc,
                                                              int nbpoly)
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
      TAU.setValue(m, n, (m == n) ? 1. - TAU.getValue(m, n) :
                                    -TAU.getValue(m, n));

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
GSTLEARN_EXPORT VectorDouble hermiteLognormal(double mean,
                                              double sigma,
                                              int nbpoly)
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
GSTLEARN_EXPORT double hermiteSeries(const VectorDouble &an,
                                     const VectorDouble &hn)
{
  double value = 0.;
  for (int ih = 0; ih < (int) hn.size(); ih++)
  {
    value += an[ih] * hn[ih];
  }
  return value;
}

/**
 *
 * @param y value where the Hermite polynomials are calculated
 * @param nbpoly Number of Polynomial functions
 * @return Coefficients of f(Y)=Yp in Hermite polynomials
 */
GSTLEARN_EXPORT VectorDouble hermiteFunction(double y, int nbpoly)
{
// TODO a corriger car le texte actuel ne fait pas sens.

  /* Calculate the Hermite polynomials for value 'y' */

  VectorDouble hn = hermitePolynomials(y, 1., nbpoly);
  VectorDouble coeff(nbpoly);

  /* Calculate the coefficients */

  double dg = law_df_gaussian(y);
  double pg = law_cdf_gaussian(y);
  hn[0] = dg + y * pg;
  hn[1] = pg - 1.;
  for (int i = 2; i < nbpoly; i++)
    hn[i] = dg * (y * hn[i - 1] / sqrt(i) + hn[i - 2] / sqrt(i * (i - 1)));
  return coeff;
}

/**
 *
 * @param yk Gaussian value
 * @param sk Standard Deviation coefficient
 * @param phi Array of Hermite coefficient
 * @return Calculate: E[Z^2| z1,...,zn]
 */
GSTLEARN_EXPORT VectorDouble hermiteEvaluateZ2(VectorDouble yk,
                                               VectorDouble sk,
                                               const VectorDouble &phi)
{
  int nech = static_cast<int>(yk.size());
  int nbpoly = static_cast<int>(phi.size());
  double log2 = log(2.);
  VectorDouble tab(nech, 0);
  VectorDouble fact(nbpoly);
  ut_log_factorial(nbpoly, fact.data());

  for (int iech = 0; iech < nech; iech++)
  {
    double var = sk[iech] * sk[iech];
    double logs = 0.5 * log(1 - var);

    /* Calculate the Hermite and log-factorial terms */

    VectorDouble hn = hermitePolynomials(yk[iech], 1., nbpoly);

    /* Loop on the different terms */

    double z = 0.;
    for (int n = 0; n < nbpoly; n++)
    {
      double phin = phi[n];
      if (ABS(phin) < EPSILON10) continue;
      for (int m = 0; m < nbpoly; m++)
      {
        double phim = phi[m];
        if (ABS(phim) < EPSILON10) continue;

        double coef = 0.;
        double fkmn = fact[m] + fact[n];
        for (int k1 = 0; 2 * k1 <= n; k1++)
          for (int k2 = 0; 2 * k2 <= m; k2++)
          {
            int k1p2 = k1 + k2;
            double fk12 = fact[k1] + fact[k2];
            double c1 = k1p2 * (2. * logs - log2);
            for (int p1 = 0; p1 <= n - 2 * k1; p1++)
            {
              double hnp1 = hn[p1];
              if (ABS(hnp1) < EPSILON10) continue;
              for (int p2 = 0; p2 <= m - 2 * k2; p2++)
              {
                double hnp2 = hn[p2];
                if (ABS(hnp2) < EPSILON10) continue;
                double fp12 = fact[p1] + fact[p2];
                for (int q = 0; q <= n - 2 * k1 - p1 && q <= m - 2 * k2 - p2;
                    q++)
                {
                  double a = c1 + 2. * q * logs + fkmn - fk12 - fp12 - fact[q];
                  coef += exp(a) * hnp1 * hnp2;
                }
              }
            }
          }
        z += coef * phin * phim;
      }
      tab[iech] = z;
    }
  }
  return tab;
}

GSTLEARN_EXPORT double hermiteEvaluateZ2(double yk,
                                         double sk,
                                         const VectorDouble &phi)
{
  int nbpoly = static_cast<int>(phi.size());
  double log2 = log(2.);
  VectorDouble fact(nbpoly);
  ut_log_factorial(nbpoly, fact.data());

  double tab = 0.;
  double var = sk * sk;
  double logs = 0.5 * log(1 - var);

  /* Calculate the Hermite and log-factorial terms */

  VectorDouble hn = hermitePolynomials(yk, 1., nbpoly);

  /* Loop on the different terms */

  double z = 0.;
  for (int n = 0; n < nbpoly; n++)
  {
    double phin = phi[n];
    if (ABS(phin) < EPSILON10) continue;
    for (int m = 0; m < nbpoly; m++)
    {
      double phim = phi[m];
      if (ABS(phim) < EPSILON10) continue;

      double coef = 0.;
      double fkmn = fact[m] + fact[n];
      for (int k1 = 0; 2 * k1 <= n; k1++)
        for (int k2 = 0; 2 * k2 <= m; k2++)
        {
          int k1p2 = k1 + k2;
          double fk12 = fact[k1] + fact[k2];
          double c1 = k1p2 * (2. * logs - log2);
          for (int p1 = 0; p1 <= n - 2 * k1; p1++)
          {
            double hnp1 = hn[p1];
            if (ABS(hnp1) < EPSILON10) continue;
            for (int p2 = 0; p2 <= m - 2 * k2; p2++)
            {
              double hnp2 = hn[p2];
              if (ABS(hnp2) < EPSILON10) continue;
              double fp12 = fact[p1] + fact[p2];
              for (int q = 0; q <= n - 2 * k1 - p1 && q <= m - 2 * k2 - p2; q++)
              {
                double a = c1 + 2. * q * logs + fkmn - fk12 - fp12 - fact[q];
                coef += exp(a) * hnp1 * hnp2;
              }
            }
          }
        }
      z += coef * phin * phim;
    }
    tab = z;
  }
  return tab;
}

