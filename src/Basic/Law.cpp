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
#include "geoslib_old_f.h"
#include "Basic/Law.hpp"

#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/MathFunc.hpp"
#include "Basic/VectorHelper.hpp"

#include <math.h>
#include <random>

static int Random_factor     = 105;
static int Random_congruent  = 20000159;
static int Random_value      = 43241421;
static bool Random_Old_Style = true;
std::mt19937 Random_gen;

/*! \cond */
#define TABIN_BY_COL(iech,ivar)       (tabin [(ivar) * nechin  + (iech)])
#define TABIN_BY_LINE(iech,ivar)      (tabin [(iech) * nvarin  + (ivar)])
#define TABOUT_BY_COL(iech,ivar)      (tabout[(ivar) * nechout + (iech)])
#define TABOUT_BY_LINE(iech,ivar)     (tabout[(iech) * nvarout + (ivar)])
#define CONSTS(iconst,ivar1)          (consts[(iconst) * nvar1 + (ivar1)])
/*! \endcond */

/**
 * Set the type of Usage for Random Number Generation
 * @param style true for using Old Style; false for using New Style
 */
void law_set_old_style(bool style)
{
  Random_Old_Style = style;
}

/*****************************************************************************/
/*!
 **  read the seed for the random number generator
 **
 ** \return  The current value of the seed (integer)
 **
 *****************************************************************************/
int law_get_random_seed(void)

{
  return Random_value;
}

/*****************************************************************************/
/*!
 **  Sets the seed for the random number generator
 **
 ** \param[in]  seed the new value given to the seed
 **
 *****************************************************************************/
void law_set_random_seed(int seed)

{
  if (seed > 0)
  {
    Random_value = seed;
    if (! Random_Old_Style)
      Random_gen.seed((unsigned) seed);
  }
  return;
}

/*****************************************************************************/
/*!
 **  Draw a random number according to a uniform distribution
 **
 ** \return  Uniform random value within an interval
 **
 ** \param[in]  mini  minimum value
 ** \param[in]  maxi  maximum value
 **
 *****************************************************************************/
double law_uniform(double mini, double maxi)

{
  double value = 0.;

  if (Random_Old_Style)
  {
    unsigned int random_product;
    random_product = Random_factor * Random_value;
    Random_value = random_product % Random_congruent;
    value = (double) Random_value / (double) Random_congruent;
    value = mini + value * (maxi - mini);
  }
  else
  {
    std::uniform_real_distribution<double> d{mini,maxi};
    value = d(Random_gen);
  }
  return (value);
}

/*****************************************************************************/
/*!
 **  Draw an integer random number according to a uniform distribution
 **
 ** \return  Integer Uniform random value within an interval
 **
 ** \param[in]  mini  minimum value
 ** \param[in]  maxi  maximum value
 **
 *****************************************************************************/
int law_int_uniform(int mini, int maxi)
{
  double rndval;
  int number, rank;

  number = maxi - mini + 1;
  rndval = law_uniform(0., (double) number);
  rank = (int) floor(rndval);
  return (rank + mini);
}

/*****************************************************************************/
/*!
 **  Generate random numbers according to a gaussian distribution
 **
 ** \return  Gaussian random value
 **
 ** \param[in]  mean Mean of the Normal Distribution
 ** \param[in]  sigma Standard deviation of the Normal Distribution
 **
 *****************************************************************************/
double law_gaussian(double mean, double sigma)

{
  double value = 0.;

  if (Random_Old_Style)
  {
    double random1 = law_uniform();
    double random2 = law_uniform(0., 2. * GV_PI);
    value = sqrt(-2. * log(random1)) * cos(random2);
    value = value * sigma + mean;
  }
  else
  {
    std::normal_distribution<double> d{mean,sigma};
    value = d(Random_gen);
  }
  return value;
}

/*****************************************************************************/
/*!
 **  Generate random numbers according to exponential distribution
 **
 ** \return  Exponential random value
 **
 ** \param[in]  lambda Parameter of exponential distribution
 **
 *****************************************************************************/
double law_exponential(double lambda)

{
  double value = 0;

  if (Random_Old_Style)
  {
    value = law_uniform(0., 1.);
    value = -log(value) / lambda;
  }
  else
  {
    std::exponential_distribution<double> d(lambda);
    value = d(Random_gen);
  }
  return value;
}

/*****************************************************************************/
/*!
 **  Generate random numbers according to a gamma distribution
 **
 ** \return  Gamma random value
 **
 ** \param[in]  alpha parameter of the gamma distribution
 ** \param[in]  beta  Second parameter of the Gamma distribution
 **
 *****************************************************************************/
double law_gamma(double alpha, double beta)

{
  double value = 0.;
  if (alpha <= 0.) return (TEST);

  if (Random_Old_Style)
  {

    if (fabs((double) alpha - 1.) < 0.00001)
      return (-log(law_uniform(0., 1.)));
    else if (alpha > 1.)
    {
      double c1 = alpha - 1.;
      double c2 = alpha + c1;
      double c3 = sqrt(c2);
      double t;
      do
      {
        t = c3 * tan(GV_PI * (law_uniform(-0.5, 0.5)));
        value = c1 + t;
      }
      while (value < 0
          || law_uniform(0., 1.) > exp(c1 * log(value / c1) - t + log(1 + t * t / c2)));
      return (value);
    }
    else
    {
      double c1 = 1. + alpha / GV_EE;
      double c2 = 1 / alpha;
      double c3 = alpha - 1;
      int test;
      do
      {
        double v = law_uniform(0., 1.);
        value = c1 * law_uniform(0., 1.);
        if (value <= 1)
        {
          value = pow(value, c2);
          test = (v >= exp(-value));
        }
        else
        {
          value = -log((c1 - value) * c2);
          test = (log(v) > c3 * log(value));
        }
      }
      while (test);
      return (value);
    }
  }
  else
  {
    std::gamma_distribution<double> d(alpha, beta);
    value = d(Random_gen);
    return value;
  }
}

/*****************************************************************************/
/*!
 **  Generate random numbers according to a standard stable distribution
 **
 ** \return  Stable value with std. parameters (beta=gamma=1,delta=0,alpha!=1)
 **
 ** \param[in]  alpha  value of the alpha parameter
 **
 *****************************************************************************/
double law_stable_standard_abgd(double alpha)
{
  double unif, expo, temp, ialpha, b, res;

  b = GV_PI / 2;
  unif = law_uniform(-b, b);
  expo = law_exponential(1.);
  ialpha = 1. / alpha;
  if (alpha > 1) b *= (1 - 2. / alpha);
  temp = alpha * (unif + b);
  res = sin(temp) / pow(cos(unif), ialpha);
  res *= pow(cos(unif - temp) / expo, ialpha - 1.);
  res = (!FFFF(unif) && !FFFF(expo)) ? res : TEST;
  return (res);
}

/*****************************************************************************/
/*!
 **  Generate random numbers according to a standard stable distribution
 **
 ** \return  Stable value with standard parameters (gamma=1,delta=0,alpha!=1)
 **
 ** \param[in]  alpha  value of the alpha parameter
 ** \param[in]  beta   value of the beta parameter
 **
 *****************************************************************************/
double law_stable_standard_agd(double alpha, double beta)
{
  double unif, expo, unif_norm, ialpha, temp, temp1, b, temp2, temp3, res;

  temp = alpha * GV_PI / 2;
  ialpha = 1. / alpha;
  unif = law_uniform(-temp, temp);
  expo = law_exponential(1.);
  unif_norm = unif * ialpha;
  temp1 = beta * tan(temp);
  b = atan(temp1);
  temp2 = pow(1 + temp1 * temp1, ialpha / 2);
  temp3 = unif + b;
  res = temp2 * sin(temp3) / pow(cos(unif_norm), ialpha);
  res *= pow(cos(unif_norm - temp3) / expo, ialpha - 1);
  res = (!FFFF(unif) && !FFFF(expo)) ? res :
                                       TEST;
  return (res);
}

/*****************************************************************************/
/*!
 ** Generate random numbers according to a standard stable distribution (alpha=1)
 **
 ** \return  Stable value with standard parameters (alpha=gamma=1,delta=0)
 **
 ** \param[in]  beta   value of the beta parameter
 *****************************************************************************/
double law_stable_standard_a1gd(double beta)
{
  double unif, expo, temp, temp2, res;

  temp = GV_PI / 2;
  unif = law_uniform(-temp, temp);
  expo = law_exponential(1.);
  temp2 = temp + beta * unif;
  res = temp2 * tan(unif);
  res -= beta * log(expo * cos(unif) / temp2);
  res /= temp;
  res = (!FFFF(unif) && !FFFF(expo)) ? res :
                                       TEST;
  return (res);
}

/*****************************************************************************/
/*!
 **  Generate random numbers according to a stable distribution (alpha != 1)
 **
 ** \return  Stable value with unit parameters

 ** \param[in]  alpha  value of the alpha parameter
 ** \param[in]  beta   value of the beta  parameter
 ** \param[in]  gamma  value of the gamma parameter
 ** \param[in]  delta  value of the delta parameter
 *****************************************************************************/
double law_stable_a(double alpha, double beta, double gamma, double delta)
{
  double stable, res;
  stable = law_stable_standard_agd(alpha, beta);
  res = pow(gamma, 1. / alpha) * stable + gamma * delta;
  res = (!FFFF(stable)) ? res :
                          TEST;
  return (res);
}

/*****************************************************************************/
/*!
 **  Generate random numbers according to a stable distribution (alpha=1)
 **
 ** \return  Stable value with unit parameters

 ** \param[in]  beta   value of the beta parameter
 ** \param[in]  gamma  value of the gamma parameter
 ** \param[in]  delta  value of the delta parameter
 *****************************************************************************/
double law_stable_a1(double beta, double gamma, double delta)
{
  double stable, res;
  stable = law_stable_standard_a1gd(beta);
  res = gamma * (stable + delta + 2. / GV_PI * beta * log(gamma));
  res = (!FFFF(stable)) ? res :
                          TEST;
  return (res);
}

/*****************************************************************************/
/*!
 **  Generate random numbers according to a stable distribution
 **
 ** \return  Stable value with unit parameters
 **
 ** \param[in]  alpha  value of the alpha parameter
 ** \param[in]  beta   value of the beta parameter
 ** \param[in]  gamma  value of the gamma parameter
 ** \param[in]  delta  value of the delta parameter
 **
 *****************************************************************************/
double law_stable(double alpha, double beta, double gamma, double delta)
{
  double res;
  if (alpha == 1)
    res = law_stable_a1(beta, gamma, delta);
  else
    res = law_stable_a(alpha, beta, gamma, delta);
  return (res);
}

/*****************************************************************************/
/*!
 **  Generate random numbers according to a beta distribution (first kind)
 **
 ** \return  Beta random value (first kind)
 **
 ** \param[in]  parameter1 first parameter of the beta distribution
 ** \param[in]  parameter2 first parameter of the beta distribution
 **
 *****************************************************************************/
double law_beta1(double parameter1, double parameter2)
{
  double a, b, res;

  a = law_gamma(parameter1,1.);
  b = law_gamma(parameter2,1.);

  res = (!FFFF(a) && !FFFF(b)) ? a / (a + b) : TEST;
  return (res);
}

/*****************************************************************************/
/*!
 **  Generate random numbers according to a beta distribution (second kind)
 **
 ** \return  Beta random value (second kind)
 **
 ** \param[in]  parameter1 first parameter of the beta distribution
 ** \param[in]  parameter2 first parameter of the beta distribution
 **
 *****************************************************************************/
double law_beta2(double parameter1, double parameter2)
{
  double a, b, res;

  a = law_gamma(parameter1,1.);
  b = law_gamma(parameter2,1.);

  res = (!FFFF(a) && !FFFF(b)) ? a / b : TEST;
  return (res);
}

/*****************************************************************************/
/*!
 **  Density function of a gaussian distribution
 **
 ** \return  Gaussian density function
 **
 ** \param[in]  value raw value
 **
 *****************************************************************************/
double law_df_gaussian(double value)

{
  double val;

  val = (ABS(value) > 10) ? 0 : exp(-value * value / 2.);
  return (val / sqrt(2. * GV_PI));
}

/*****************************************************************************/
/*!
 **  Density function of a (non-normalized) gaussian distribution
 **
 ** \return  Gaussian density function
 **
 ** \param[in]  value Raw value
 ** \param[in]  mean  Mean value
 ** \param[in]  std   Standard deviation
 **
 *****************************************************************************/
double law_dnorm(double value, double mean, double std)
{
  double val, center;

  center = (value - mean) / std;

  if (ABS(center) > 10)
    val = 0.;
  else
    val = exp(-center * center / 2.);

  return (val / sqrt(2. * GV_PI) / std);
}

/*****************************************************************************/
/*!
 **  Cumulated density function of a gaussian distribution
 **
 ** \return  Gaussian cumulated density function
 **
 ** \param[in]  value raw value
 **
 ** \remark  Handbook P932 (26.2.17)  precision <7.5 E-08
 **
 *****************************************************************************/
double law_cdf_gaussian(double value)

{
  static double b[] = { 0.319381530,
                        -0.356563782,
                        1.781477937,
                        -1.821255978,
                        1.330274429 };
  static double p = 0.2316419;
  double u, v, x, t;
  int i;

  x = fabs(value);
  t = 1. / (1. + p * x);
  u = 0.;
  if (x < 10)
  {
    v = 1.;
    for (i = 0; i < 5; i++)
    {
      v *= t;
      u += b[i] * v;
    }
    u = u * exp(-x * x / 2.) / sqrt(2. * GV_PI);
  }

  if (value >= 0) u = 1. - u;

  return (u);
}

/*****************************************************************************/
/*!
 **  Inverse cumulated density function of a gaussian distribution
 **
 ** \return  Inverse of gaussian cumulated density function
 **
 ** \param[in]  value cumulative density
 **
 *****************************************************************************/
double law_invcdf_gaussian(double value)

{
  static double b[] = { 0.319381530,
                        -0.356563782,
                        1.781477937,
                        -1.821255978,
                        1.330274429 };
  static double c[] = { 2.515517, 0.802853, 0.010328 };
  static double d[] = { 1.432788, 0.189269, 0.001308 };
  static double p = 0.2316419;
  double xmin, xmax, v, t, freq, x = 0.;

  if (value <= 0.) return (-10.);
  if (value >= 1.) return (10.);
  v = (value < 0.5) ? 1 - value :
                      value;
  t = sqrt(-2 * log(1 - v));
  xmin = t
      - (c[0] + t * (c[1] + t * c[2])) / (1.
          + t * (d[0] + t * (d[1] + t * d[2])));
  xmin = xmin - 0.001;
  xmax = xmin + 0.002;

  while ((xmax - xmin) > 0.0000001)
  {
    x = (xmin + xmax) / 2.;
    t = 1. / (1. + p * x);
    freq = t * (b[0] + t * (b[1] + t * (b[2] + t * (b[3] + t * b[4]))));
    if (freq * exp(-x * x / 2) / sqrt(2. * GV_PI) > 1. - v)
      xmin = x;
    else
      xmax = x;
  }

  if (value < 0.5)
    return (-x);
  else
    return (x);
}

/*****************************************************************************/
/*!
 **  Generates a gaussian value which lies in an interval
 **
 ** \return  The gaussian value
 **
 ** \param[in]  binf lower bound of the interval
 ** \param[in]  bsup upper bound of the interval
 **
 *****************************************************************************/
double law_gaussian_between_bounds(double binf, double bsup)
{
  double atab[4], btab[4], ptab[4];
  double a, b, aa, bb, total, a2, b2, c2, u, x;
  int k, n, isim, type, ok, rank, itab[4];
  static double seuil = 2.;
  static double large = 20.;
  static double sqe = 1.6487212707;

  x = 0.;

  /* Loop on the intervals */

  n = 0;
  a = (FFFF(binf)) ? -large : binf;
  b = (FFFF(bsup)) ? large : bsup;
  aa = a;
  bb = b;

  if (aa < -seuil)
  {
    atab[n] = aa;
    btab[n] = MIN(bb, -seuil);
    itab[n++] = 1;
    aa = -seuil;
    if (aa >= bb) goto label_norme;
  }
  if (aa < 0)
  {
    atab[n] = aa;
    btab[n] = MIN(bb, 0.);
    itab[n++] = 2;
    aa = 0.;
    if (aa >= bb) goto label_norme;
  }
  if (aa < seuil)
  {
    atab[n] = aa;
    btab[n] = MIN(bb, seuil);
    itab[n++] = 3;
    aa = seuil;
    if (aa >= bb) goto label_norme;
  }
  atab[n] = aa;
  btab[n] = bb;
  itab[n++] = 4;

  label_norme:
  total = 0.;
  double wgt;
  for (k = 0; k < n; k++)
  {
    aa = atab[k];
    bb = btab[k];
    type = itab[k];
    wgt = 0.;
    if (ABS(aa - bb) > 0.) switch (type)
    {
      case 1:
        wgt = (exp(-aa * aa / 2.) - exp(-bb * bb / 2.)) / bb;
        break;
      case 2:
        wgt = sqe * (exp(bb) - exp(aa));
        break;
      case 3:
        wgt = sqe * (exp(-aa) - exp(-bb));
        break;
      case 4:
        wgt = (exp(-aa * aa / 2.) - exp(-bb * bb / 2.)) / aa;
        break;
    }
    total += wgt;
    ptab[k] = total;
  }

  /* Normalization */

  if (total <= 0.)
  {
    rank = (int) ((double) n * law_uniform(0., 1.));
    x = atab[rank];
    return (x);
  }

  for (k = 0; k < n; k++)
    ptab[k] /= total;

  /* Simulation by acceptance and reject method */

  ok = 1;
  while (ok)
  {
    isim = 0;
    u = law_uniform(0., 1.);
    while (ptab[isim] < u)
      isim++;

    /* Simulate a value in the interval */

    type = itab[isim];
    aa = atab[isim];
    bb = btab[isim];
    a2 = aa * aa;
    b2 = bb * bb;

    u = law_uniform(0., 1.);
    switch (type)
    {
      case 1:
        c2 = 1. - exp((b2 - a2) / 2.);
        x = -sqrt(b2 - 2. * log(1. - u * c2));
        break;

      case 2:
        x = log(exp(aa) * (1. - u) + exp(bb) * u);
        break;

      case 3:
        x = -log(exp(-aa) * (1. - u) + exp(-bb) * u);
        break;

      case 4:
        c2 = 1. - exp((a2 - b2) / 2.);
        x = sqrt(a2 - 2. * log(1. - u * c2));
        break;
    }

    /* Accept or reject the value */

    u = law_uniform(0., 1.);
    switch (type)
    {
      case 1:
        if (x * u >= bb) ok = 0;
        break;

      case 2:
        if (log(u) <= -(x + 1.) * (x + 1.) / 2.) ok = 0;
        break;

      case 3:
        if (log(u) <= -(x - 1.) * (x - 1.) / 2.) ok = 0;
        break;

      case 4:
        if (x * u <= aa) ok = 0;
        break;
    }
  }
  return (x);
}

/*****************************************************************************/
/*!
 **  Density function of a bigaussian distribution
 **
 ** \return  Gaussian density function
 **
 ** \param[in]  vect   Array of values (Dimension = 2)
 ** \param[in]  mean   Array of means (Dimension = 2)
 ** \param[in]  correl Correlation matrix (Dimension: 2*2)
 **
 *****************************************************************************/
double law_df_bigaussian(VectorDouble &vect,
                         VectorDouble &mean,
                         MatrixSquareSymmetric& correl)
{
  VectorDouble xc(2);
  xc[0] = vect[0] - mean[0];
  xc[1] = vect[1] - mean[1];
  double detv = correl.getValue(0,0) * correl.getValue(1,1) - correl.getValue(0,1) * correl.getValue(1,0);
  double det2 = (correl.getValue(1,1) * xc[0] * xc[0] - 2. * correl.getValue(0,1) * xc[0] * xc[1] + correl.getValue(0,0) * xc[1] * xc[1]);
  double logres = 2. * log(2. * GV_PI) + log(detv) + det2 / detv;
  double density = exp(-logres / 2.);
  return (density);
}

/*****************************************************************************/
/*!
 **  Density function of a quadrigaussian distribution
 **
 ** \return  Gaussian density function
 **
 ** \param[in]  vect   Array of values (Dimension = nvar)
 ** \param[in]  correl Correlation matrix (Dimension: nvar*nvar)
 **
 *****************************************************************************/
double law_df_quadgaussian(VectorDouble& vect, MatrixSquareSymmetric& correl)
{
  int nvar = (int) vect.size();
  double density = -2. * log(2 * GV_PI);

  if (correl.computeEigen()) return TEST;

  VectorDouble eigval = correl.getEigenValues();
  for (int ivar = 0; ivar < nvar; ivar++)
    density -= 0.5 * log(eigval[ivar]);

  MatrixSquareSymmetric invcor = correl;
  if (invcor.invert()) return TEST;

  density -= 0.5 * invcor.normVec(vect);
  density = exp(density);
  return density;
}

/*****************************************************************************/
/*!
 **  Density function of a multigaussian distribution
 **
 ** \return  Gaussian density function
 **
 ** \param[in]  vect   Array of values (Dimension = nvar)
 ** \param[in]  correl Correlation matrix (Dimension: nvar*nvar)
 **
 *****************************************************************************/
double law_df_multigaussian(VectorDouble& vect, MatrixSquareSymmetric& correl)

{
  int nvar = (int) vect.size();
  double density = -0.5 * nvar * log(2 * GV_PI);

  if (correl.computeEigen()) return TEST;
  VectorDouble eigval = correl.getEigenValues();

  for (int i = 0; i < nvar; i++)
    density -= 0.5 * log(eigval[i]);

  MatrixSquareSymmetric invcor(correl);
  if (invcor.invert()) return TEST;

  density -= invcor.normVec(vect);
  density = exp(density);
  return (density);
}

VectorDouble law_df_poisson_vec(VectorInt is, double parameter)
{
  int size = (int) is.size();
  VectorDouble res(size);
  for (int ii = 0; ii < size; ii++)
    res[ii] = law_df_poisson(is[ii], parameter);
  return res;
}

double law_df_poisson(int i, double parameter)
{
  return (exp(-parameter) * pow(parameter, i) / ut_factorial(i));
}

/****************************************************************************/
/*!
 ** Generate random number according to a poisson distribution
 **
 ** \return Poisson random value
 **
 ** \param[in] parameter parameter of the Poisson distribution
 **
 ** \remarks  Method Ahrens-Dieter (1973)
 **
 *****************************************************************************/
int law_poisson(double parameter)
{
  if (Random_Old_Style)
  {
    double x, p, q;
    int n, ok;

    int k = 0;
    double t = parameter;

    while (t >= 16)
    {
      n = (int) floor(0.875 * t);
      x = law_gamma((double) n, 1.);
      if (FFFF(x)) return (ITEST);

      if (x > t)
      {
        p = t / x;
        ok = 1;
        while (ok)
        {
          if (law_uniform(0., 1.) <= p) k++;
          if (n-- <= 1) ok = 0;
        }
        return (k);
      }
      else
      {
        k += n;
        t -= x;
      }
    }

    p = 1;
    q = exp(-t);

    ok = 1;
    while (ok)
    {
      p *= law_uniform(0., 1.);
      if (p < q) ok = 0;
      k++;
    }
    return (k - 1);
  }
  else
  {
    std::poisson_distribution<int> d(parameter);
    return d(Random_gen);
  }
}

/****************************************************************************/
/*!
 ** Define a random path
 **
 ** \param[in] nech       : Number of samples
 **
 *****************************************************************************/
VectorInt law_random_path(int nech)
{
  VectorDouble order(nech);
  VectorInt path(nech);

  for (int i = 0; i < nech; i++)
  {
    path[i] = i;
    order[i] = law_uniform(0., 1.);
  }
  VH::arrangeInPlace(0, path, order, true, nech);
  return path;
}

/*****************************************************************************/
/*!
 **  Generates a binomial value
 **
 ** \return  The binomial value
 **
 ** \param[in]  n    Number of trials
 ** \param[in]  p    Event probability
 **
 *****************************************************************************/
int law_binomial(int n, double p)
{
  const double q = 1 - p;
  if (n * p < 30.0) /* Algorithm BINV */
  {
    const double s = p / q;
    const double a = (n + 1) * s;
    double r = exp(n * log(q)); /* pow() causes a crash on AIX */
    int x = 0;
    double u = law_uniform(0., 1.);
    while (1)
    {
      if (u < r) return x;
      u -= r;
      x++;
      r *= (a / x) - s;
    }
  }
  else /* Algorithm BTPE */
  {
    /* Step 0 */
    const double fm = n * p + p;
    const int m = (int) fm;
    const double p1 = floor(2.195 * sqrt(n * p * q) - 4.6 * q) + 0.5;
    const double xm = m + 0.5;
    const double xl = xm - p1;
    const double xr = xm + p1;
    const double c = 0.134 + 20.5 / (15.3 + m);
    const double a = (fm - xl) / (fm - xl * p);
    const double b = (xr - fm) / (xr * q);
    const double lambdal = a * (1.0 + 0.5 * a);
    const double lambdar = b * (1.0 + 0.5 * b);
    const double p2 = p1 * (1 + 2 * c);
    const double p3 = p2 + c / lambdal;
    const double p4 = p3 + c / lambdar;
    while (1)
    {
      /* Step 1 */
      int y;
      int k;
      double u = law_uniform(0., 1.);
      double v = law_uniform(0., 1.);
      u *= p4;
      if (u <= p1) return (int) (xm - p1 * v + u);
      /* Step 2 */
      if (u > p2)
      {
        /* Step 3 */
        if (u > p3)
        {
          /* Step 4 */
          y = (int) (xr - log(v) / lambdar);
          if (y > n) continue;
          /* Go to step 5 */
          v = v * (u - p3) * lambdar;
        }
        else
        {
          y = (int) (xl + log(v) / lambdal);
          if (y < 0) continue;
          /* Go to step 5 */
          v = v * (u - p2) * lambdal;
        }
      }
      else
      {
        const double x = xl + (u - p1) / c;
        v = v * c + 1.0 - fabs(m - x + 0.5) / p1;
        if (v > 1) continue;
        /* Go to step 5 */
        y = (int) x;
      }
      /* Step 5 */
      /* Step 5.0 */
      k = ABS(y - m);
      if (k > 20 && k < 0.5 * n * p * q - 1.0)
      {
        /* Step 5.2 */
        double rho = (k / (n * p * q))
            * ((k * (k / 3.0 + 0.625) + 0.1666666666666) / (n * p * q) + 0.5);
        double t = -k * k / (2 * n * p * q);
        double A = log(v);
        if (A < t - rho)
          return y;
        else if (A > t + rho)
          continue;
        else
        {
          /* Step 5.3 */
          double x1 = y + 1;
          double f1 = m + 1;
          double z = n + 1 - m;
          double w = n - y + 1;
          double x2 = x1 * x1;
          double f2 = f1 * f1;
          double z2 = z * z;
          double w2 = w * w;
          if (A > xm * log(f1 / x1)
              + (n - m + 0.5) * log(z / w)
              + (y - m) * log(w * p / (x1 * q))
              + (13860. - (462. - (132. - (99. - 140. / f2) / f2) / f2) / f2) / f1
                / 166320.
              + (13860. - (462. - (132. - (99. - 140. / z2) / z2) / z2) / z2) / z
                / 166320.
              + (13860. - (462. - (132. - (99. - 140. / x2) / x2) / x2) / x2) / x1
                / 166320.
              + (13860. - (462. - (132. - (99. - 140. / w2) / w2) / w2) / w2) / w
                / 166320.) continue;
          return y;
        }
      }
      else
      {
        /* Step 5.1 */
        int i;
        const double s = p / q;
        const double aa = s * (n + 1);
        double f = 1.0;
        for (i = m; i < y; f *= (aa / (++i) - s))
          ;
        for (i = y; i < m; f /= (aa / (++i) - s))
          ;
        if (v > f) continue;
        return y;
      }
    }
  }
  /* Never get here */
  return -1;
}

/****************************************************************************/
/*!
 **  Sample a multivariate empirical distribution
 **
 ** \return  Address to the newly created array
 **
 ** \param[in]  tabin   Input array
 ** \param[in]  mode    Describes the way 'tabin' and 'tabout' are filled
 **                     1: by column; 2: by row
 ** \param[in]  nvar    Number of variables (input and output)
 ** \param[in]  nechin  Number of samples in the input array
 ** \param[in]  nechout Number of created samples
 ** \param[in]  niter   Maximum number of iterations
 ** \param[in]  nconst  Number of constraints
 ** \param[in]  consts  Array of constraints (optional)
 **                     (Dimension: nconst * (nvar+1))
 **                     This array is entered by line
 ** \param[in]  seed    Value for the seed generator
 ** \param[in]  percent Dimension of the convolution kernel
 **                     expressed as a percentage of the dispersion st. dev.
 **                     Should be between 0 and 100
 **
 ** \remarks The input array must be isotopic (otherwise, an error is issued)
 **
 ** \remarks The resulting array if dimensionned to nvar * nechout
 ** \remarks The created (and returned) array must be freed by the calling
 ** \remarks function
 **
 ** \remarks Consider nvar1 = nvar + 1
 ** \remarks Sample temp[1:nvar1] is authorized for Constraint 'iconst' if:
 ** \remarks  Sum_ivar1^{1:nvar1) consts[iconst,ivar1) * temp[ivar1] > 0
 **
 *****************************************************************************/
double* law_exp_sample(double *tabin,
                       int mode,
                       int nvar,
                       int nechin,
                       int nechout,
                       int niter,
                       int nconst,
                       double *consts,
                       int seed,
                       double percent)
{
  double *tabout, *mean, *stdv, *mini, *maxi, *temp, value, rab, total;
  int error, iechin, selec, nvarin, nvarout, nvar1, flag_cont, flag_ok;

  /* Initializations */

  error = 1;
  tabout = mean = stdv = temp = mini = maxi = nullptr;
  law_set_random_seed(seed);
  nvarin = nvarout = nvar;
  nvar1 = nvar + 1;

  /* Internal core allocation */

  temp = (double*) mem_alloc(sizeof(double) * nvar1, 1);
  mean = (double*) mem_alloc(sizeof(double) * nvarin, 1);
  stdv = (double*) mem_alloc(sizeof(double) * nvarin, 1);
  mini = (double*) mem_alloc(sizeof(double) * nvarin, 1);
  maxi = (double*) mem_alloc(sizeof(double) * nvarin, 1);
  for (int ivar = 0; ivar < nvarin; ivar++)
  {
    mean[ivar] = stdv[ivar] = 0.;
    mini[ivar] = +1.e30;
    maxi[ivar] = -1.e30;
  }

  /* Count the number of active isotopic samples */

  for (iechin = 0; iechin < nechin; iechin++)
  {
    for (int ivar = 0; ivar < nvarin; ivar++)
    {
      value = (mode == 1) ? TABIN_BY_COL(iechin, ivar) :
                            TABIN_BY_LINE(iechin, ivar);
      if (FFFF(value))
      {
        messerr("The sample %d of the Training Data Base", iechin + 1);
        messerr("is not isotopic. This is an error");
        goto label_end;
      }
      if (value < mini[ivar]) mini[ivar] = value;
      if (value > maxi[ivar]) maxi[ivar] = value;
      mean[ivar] += value;
      stdv[ivar] += value * value;
    }
  }

  /* Normalize the statistics */

  for (int ivar = 0; ivar < nvarin; ivar++)
  {
    mean[ivar] /= (double) nechin;
    stdv[ivar] = stdv[ivar] / (double) nechin - mean[ivar] * mean[ivar];
    stdv[ivar] = (stdv[ivar] > 0) ? sqrt(stdv[ivar]) :
                                    0.;
    stdv[ivar] *= percent / 100.;
  }

  /* Core allocation */

  tabout = (double*) mem_alloc(sizeof(double) * nechout * nvarout, 0);
  if (tabout == nullptr) goto label_end;

  /* Generate the samples */

  for (int iechout = 0; iechout < nechout; iechout++)
  {

    /* Loop on the trials */

    flag_cont = 1;
    for (int iter = 0; iter < niter && flag_cont; iter++)
    {

      /* Get the closest experimental sample (the reference) */

      selec = (int) law_uniform(1., (double) nechin);
      iechin = (int) (selec + 0.5);
      if (iechin < 0) iechin = 0;
      if (iechin >= nechin) iechin = nechin - 1;

      /* Generate values around the reference */

      for (int ivar = 0; ivar < nvarout; ivar++)
      {
        rab = stdv[ivar] * law_gaussian();
        if (mode == 1)
          temp[ivar] = TABIN_BY_COL(iechin,ivar) + rab;
        else
          temp[ivar] = TABIN_BY_LINE(iechin,ivar) + rab;
      }
      temp[nvar] = 1;

      /* Check if the generated vector is authorized or not */

      flag_ok = 1;
      if (nconst > 0 && consts != nullptr)
      {
        for (int iconst = 0; iconst < nconst && flag_ok; iconst++)
        {
          total = 0.;
          for (int ivar1 = 0; ivar1 < nvar1; ivar1++)
            total += CONSTS(iconst,ivar1) * temp[ivar1];
          if (total < 0) flag_ok = 0;
        }
      }

      /* Store the results */

      if (flag_ok)
      {
        for (int ivar = 0; ivar < nvar; ivar++)
        {
          if (mode == 1)
            TABOUT_BY_COL(iechout,ivar) = temp[ivar];
          else
            TABOUT_BY_LINE(iechout,ivar) = temp[ivar];
        }
        flag_cont = 0;
      }
    }
    if (flag_cont)
    {
      messerr("No valid sample has been found after %d iterations", niter);
      messerr("Please check the consistency between your Training Data Base");
      messerr("and the set of constraints");
      messerr("Statistics on the Training Date Base");
      for (int ivar = 0; ivar < nvar; ivar++)
      {
        messerr(" %d - Mini=%lf - Maxi=%lf - Mean=%lf - Stdv=%lf", ivar + 1,
                mini[ivar], maxi[ivar], mean[ivar], stdv[ivar]);
      }
      print_matrix("Constraints", 0, 0, nvar1, nconst, NULL, consts);
      goto label_end;
    }
  }

  /* Return code */

  error = 0;

  label_end:
  mem_free((char* ) temp);
  mem_free((char* ) mean);
  mem_free((char* ) stdv);
  mem_free((char* ) mini);
  mem_free((char* ) maxi);
  if (error) tabout = (double*) mem_free((char* ) tabout);

  return (tabout);
}

