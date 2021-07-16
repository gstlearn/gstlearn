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
#include "geoslib_e.h"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"

/*! \cond */
#define INTRESX(ic,i)        (ctables->CT[ic]->res[(i)])
#define COVAL(ctables,iconf) (ctables->cmin + iconf * ctables->dc)
#define NELEM(ctables)       ((ctables->flag_cumul) ?               \
                              ctables->ndisc + 1 : ctables->ndisc)

static double c_b11 = 1.;
/*! \endcond */

/*****************************************************************************/
/*!
**  Normal distribution probabilities accurate to 1.e-15.
**     Z = no. of standard deviations from the mean.
**   Based upon algorithm 5666 for the error function, from:
**    Hart, J.F. et al, 'Computer Approximations', Wiley 1968
**
*****************************************************************************/
static double st_mvnphi(double *z)

{
  /* System generated locals */
  double ret_val, d__1;

  /* Local variables */
  static double zabs, p, expntl;

  zabs = ABS(*z);

  /* |Z| > 37 */

  if (zabs > 37.) {
    p = 0.;
  } else {

    /* |Z| <= 37 */

    /* Computing 2nd power */
    d__1 = zabs;
    expntl = exp(-(d__1 * d__1) / 2);

    /* |Z| < CUTOFF = 10/SQRT(2) */

    if (zabs < 7.071067811865475) {
      p = expntl * ((((((zabs * .03526249659989109 + .7003830644436881)
                        * zabs + 6.37396220353165) * zabs + 33.912866078383) *
                      zabs + 112.0792914978709) * zabs + 221.2135961699311) *
                    zabs + 220.2068679123761) /
        (((((((zabs * .08838834764831844 + 1.755667163182642) * zabs +
              16.06417757920695) * zabs + 86.78073220294608) * zabs +
            296.5642487796737) * zabs + 637.3336333788311) * zabs +
          793.8265125199484) * zabs + 440.4137358247522);

      /* |Z| >= CUTOFF. */

    } else {
      p = expntl / (zabs + 1 / (zabs + 2 / (zabs + 3 / (zabs + 4 / (
                                                          zabs + .65))))) / 2.506628274631001;
    }
  }
  if (*z > 0.) {
    p = 1 - p;
  }
  ret_val = p;
  return ret_val;
}

/*****************************************************************************/
/*!
**  Multivariate Normal Probability (local function)
**
*****************************************************************************/
static void st_mvnlms(double *a,
                      double *b,
                      int    *infin,
                      double *lower,
                      double *upper)
{
  *lower = 0.;
  *upper = 1.;
  if (*infin >= 0) {
    if (*infin != 0) {
      *lower = st_mvnphi(a);
    }
    if (*infin != 1) {
      *upper = st_mvnphi(b);
    }
  }
}

/*****************************************************************************/
/*!
** st_dkswap
**
*****************************************************************************/
static void st_dkswap(double *x,
                      double *y)
{
  static double t;

  t = *x;
  *x = *y;
  *y = t;
}

/*****************************************************************************/
/*!
** Swaps rows and columns P and Q in situ, with P <= Q
**
*****************************************************************************/
static void st_rcswp(int *p,
                     int *q,
                     double *a,
                     double *b,
                     int *infin,
                     int *n,
                     double *c)
{
  /* System generated locals */
  int i__1;

  /* Local variables */
  static int i, j, ii, jj;

  /* Parameter adjustments */
  --c;
  --infin;
  --b;
  --a;

  /* Function Body */
  st_dkswap(&a[*p], &a[*q]);
  st_dkswap(&b[*p], &b[*q]);
  j = infin[*p];
  infin[*p] = infin[*q];
  infin[*q] = j;
  jj = *p * (*p - 1) / 2;
  ii = *q * (*q - 1) / 2;
  st_dkswap(&c[jj + *p], &c[ii + *q]);
  i__1 = *p - 1;
  for (j = 1; j <= i__1; ++j) {
    st_dkswap(&c[jj + j], &c[ii + j]);
  }
  jj += *p;
  i__1 = *q - 1;
  for (i = *p + 1; i <= i__1; ++i) {
    st_dkswap(&c[jj + *p], &c[ii + i]);
    jj += i;
  }
  ii += *q;
  i__1 = *n;
  for (i = *q + 1; i <= i__1; ++i) {
    st_dkswap(&c[ii + *p], &c[ii + *q]);
    ii += i;
  }
}

/*****************************************************************************/
/*!
** Subroutine to sort integration limits and determine Cholesky factor
**
*****************************************************************************/
static void st_covsrt(int *n,
                      double *lower,
                      double *upper,
                      double *correl,
                      int *infin,
                      double *y,
                      int *infis,
                      double *a,
                      double *b,
                      double *cov,
                      int *infi)
{
  /* System generated locals */
  int i__1, i__2, i__3, i__4;
  double d__1;

  /* Local variables */
  static double amin, bmin, dmin_, emin;
  static int jmin;
  static double d, e;
  static int i, j, k, l, m;
  static double sumsq, aj, bj;
  static int ii, ij, il;
  static double cvdiag, yl, yu;
  static double sum;

  /* Parameter adjustments */
  --infi;
  --cov;
  --b;
  --a;
  --y;
  --infin;
  --correl;
  --upper;
  --lower;
  amin = bmin = 0.;

  /* Function Body */
  ij = 0;
  ii = 0;
  *infis = 0;
  i__1 = *n;
  for (i = 1; i <= i__1; ++i) {
    a[i] = 0.;
    b[i] = 0.;
    infi[i] = infin[i];
    if (infi[i] < 0) {
      ++(*infis);
    } else {
      if (infi[i] != 0) {
        a[i] = lower[i];
      }
      if (infi[i] != 1) {
        b[i] = upper[i];
      }
    }
    i__2 = i - 1;
    for (j = 1; j <= i__2; ++j) {
      ++ij;
      ++ii;
      cov[ij] = correl[ii];
    }
    ++ij;
    cov[ij] = 1.;
  }

  /* First move any doubly infinite limits to innermost positions. */

  if (*infis < *n) {
    i__1 = *n - *infis + 1;
    for (i = *n; i >= i__1; --i) {
      if (infi[i] >= 0) {
        i__2 = i - 1;
        for (j = 1; j <= i__2; ++j) {
          if (infi[j] < 0) {
            st_rcswp(&j, &i, &a[1], &b[1], &infi[1], n, &cov[1]);
            goto L10;
          }
        }
      }
    L10:
      ;
    }

    /* Sort remaining limits and determine Cholesky factor. */

    ii = 0;
    i__1 = *n - *infis;
    for (i = 1; i <= i__1; ++i) {

      /* Determine the integration limits for variable with minimum */
      /* expected probability and interchange that variable with Ith. */

      dmin_ = 0.;
      emin = 1.;
      jmin = i;
      cvdiag = 0.;
      ij = ii;
      i__2 = *n - *infis;
      for (j = i; j <= i__2; ++j) {
        if (cov[ij + j] > 1e-10) {
          sumsq = sqrt(cov[ij + j]);
          sum = 0.;
          i__3 = i - 1;
          for (k = 1; k <= i__3; ++k) {
            sum += cov[ij + k] * y[k];
          }
          aj = (a[j] - sum) / sumsq;
          bj = (b[j] - sum) / sumsq;
          st_mvnlms(&aj, &bj, &infi[j], &d, &e);
          if (emin + d >= e + dmin_) {
            jmin = j;
            amin = aj;
            bmin = bj;
            dmin_ = d;
            emin = e;
            cvdiag = sumsq;
          }
        }
        ij += j;
      }
      if (jmin > i) {
        st_rcswp(&i, &jmin, &a[1], &b[1], &infi[1], n, &cov[1]);
      }
      cov[ii + i] = cvdiag;

      /* Compute Ith column of Cholesky factor. */
      /* Compute expected value for Ith integration variable and */
      /* scale Ith covariance matrix row and limits. */

      if (cvdiag > 0.) {
        il = ii + i;
        i__2 = *n - *infis;
        for (l = i + 1; l <= i__2; ++l) {
          cov[il + i] /= cvdiag;
          ij = ii + i;
          i__3 = l;
          for (j = i + 1; j <= i__3; ++j) {
            cov[il + j] -= cov[il + i] * cov[ij + i];
            ij += j;
          }
          il += l;
        }
        if (emin > dmin_ + 1e-10) {
          yl = 0.;
          yu = 0.;
          if (infi[i] != 0) {
            /* Computing 2nd power */
            d__1 = amin;
            yl = -exp(-(d__1 * d__1) / 2) / 2.506628274631001;
          }
          if (infi[i] != 1) {
            /* Computing 2nd power */
            d__1 = bmin;
            yu = -exp(-(d__1 * d__1) / 2) / 2.506628274631001;
          }
          y[i] = (yu - yl) / (emin - dmin_);
        } else {
          if (infi[i] == 0) {
            y[i] = bmin;
          }
          if (infi[i] == 1) {
            y[i] = amin;
          }
          if (infi[i] == 2) {
            y[i] = (amin + bmin) / 2;
          }
        }
        i__2 = i;
        for (j = 1; j <= i__2; ++j) {
          ++ii;
          cov[ii] /= cvdiag;
        }
        a[i] /= cvdiag;
        b[i] /= cvdiag;
      } else {
        il = ii + i;
        i__2 = *n - *infis;
        for (l = i + 1; l <= i__2; ++l) {
          cov[il + i] = 0.;
          il += l;
        }

        /* If the covariance matrix diagonal entry is zero, */
        /* permute limits and/or rows, if necessary. */

        for (j = i - 1; j >= 1; --j) {
          if ((d__1 = cov[ii + j], ABS(d__1)) > 1e-10) {
            a[i] /= cov[ii + j];
            b[i] /= cov[ii + j];
            if (cov[ii + j] < 0.) {
              st_dkswap(&a[i], &b[i]);
              if (infi[i] != 2) {
                infi[i] = 1 - infi[i];
              }
            }
            i__2 = j;
            for (l = 1; l <= i__2; ++l) {
              cov[ii + l] /= cov[ii + j];
            }
            i__2 = i - 1;
            for (l = j + 1; l <= i__2; ++l) {
              if (cov[(l - 1) * l / 2 + j + 1] > 0.) {
                ij = ii;
                i__3 = l;
                for (k = i - 1; k >= i__3; --k) {
                  i__4 = k;
                  for (m = 1; m <= i__4; ++m) {
                    st_dkswap(&cov[ij - k + m], &cov[ij + m]
                      );
                  }
                  st_dkswap(&a[k], &a[k + 1]);
                  st_dkswap(&b[k], &b[k + 1]);
                  m = infi[k];
                  infi[k] = infi[k + 1];
                  infi[k + 1] = m;
                  ij -= k;
                }
                goto L20;
              }
            }
            goto L20;
          }
          cov[ii + j] = 0.;
        }
      L20:
        ii += i;
        y[i] = 0.;
      }
    }
  }
}

/*****************************************************************************/
/*!
** Produces the normal deviate Z corresponding to a given lower
** tail area of P.
** The hash sums below are the sums of the mantissas of the
** coefficients.   They are included for use in checking
** transcription.
**
*****************************************************************************/
static double st_phinvs(double *p)

{
  double ret_val, d__1, d__2;
  static double q, r;

  q = (*p * 2 - 1) / 2;
  if (ABS(q) <= .425) {
    r = .180625 - q * q;
    ret_val = q * (((((((r * 2509.0809287301226727 +
                         33430.575583588128105) * r + 67265.770927008700853) * r +
                       45921.953931549871457) * r + 13731.693765509461125) * r +
                     1971.5909503065514427) * r + 133.14166789178437745) * r +
                   3.387132872796366608) /
      (((((((r * 5226.495278852854561 +
             28729.085735721942674) * r + 39307.89580009271061) * r +
           21213.794301586595867) * r + 5394.1960214247511077) * r +
         687.1870074920579083) * r + 42.313330701600911252) * r + 1);

  } else {

    /* Computing MIN */

    d__1 = *p, d__2 = 1 - *p;
    r = MIN(d__1,d__2);
    if (r > 0.) {
      r = sqrt(-log(r));
      if (r <= 5.) {
        r += -1.6;
        ret_val = (((((((r * 7.7454501427834140764e-4 +
                         .0227238449892691845833) * r + .24178072517745061177)
                       * r + 1.27045825245236838258) * r +
                      3.64784832476320460504) * r + 5.7694972214606914055) *
                    r + 4.6303378461565452959) * r +
                   1.42343711074968357734) /
          (((((((r * 1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                r + .0151986665636164571966) * r +
               .14810397642748007459) * r + .68976733498510000455) *
             r + 1.6763848301838038494) * r +
            2.05319162663775882187) * r + 1);
      } else {
        r += -5.;
        ret_val = (((((((r * 2.01033439929228813265e-7 +
                         2.71155556874348757815e-5) * r +
                        .0012426609473880784386) * r + .026532189526576123093)
                      * r + .29656057182850489123) * r +
                     1.7848265399172913358) * r + 5.4637849111641143699) *
                   r + 6.6579046435011037772) /
          (((((((r * 2.04426310338993978564e-15 + 1.4215117583164458887e-7)
                * r + 1.8463183175100546818e-5) * r +
               7.868691311456132591e-4) * r +
              .0148753612908506148525) * r + .13692988092273580531)
            * r + .59983220655588793769) * r + 1);
      }
    } else {
      ret_val = 9.;
    }
    if (q < 0.) {
      ret_val = -ret_val;
    }
  }
  return ret_val;
}

/*****************************************************************************/
/*!
** A function for computing bivariate normal probabilities.
**  Yihong Ge
**  Department of Computer Science and Electrical Engineering
**  Washington State University
**  Pullman, WA 99164-2752
**   and
**  Alan Genz
**  Department of Mathematics
**  Washington State University
**  Pullman, WA 99164-3113
**  Email : alangenz@wsu.edu
**
** \return  Calculate the probability that X is larger than SH and Y is
** \return  larger than SK.
**
** \param[in] sh integration limit
** \param[in] sk integration limit
** \param[in] r  REAL, correlation coefficient
**
*****************************************************************************/
static double st_bvu(double *sh,
                     double *sk,
                     double *r)
{
  /* Initialized data */

  static struct {
    double e_1[3];
    double fill_2[7];
    double e_3[6];
    double fill_4[4];
    double e_5[10];
  } equiv_99 = {
    { .1713244923791705, .3607615730481384, .4679139345726904} ,
    {0},
    {.04717533638651177, .1069393259953183, .1600783285433464,
     .2031674267230659, .2334925365383547, .2491470458134029 },
    {0},
    {.01761400713915212, .04060142980038694, .06267204833410906,
     .08327674157670475, .1019301198172404,  .1181945319615184,
     .1316886384491766,  .1420961093183821,  .1491729864726037,
     .1527533871307259 }
  };

#define w ((double *)&equiv_99)

  static struct {
    double e_1[3];
    double fill_2[7];
    double e_3[6];
    double fill_4[4];
    double e_5[10];
  } equiv_100 = {
    { -.9324695142031522, -.6612093864662647, -.238619186083197 },
    {0},
    {-.9815606342467191, -.904117256370475,  -.769902674194305,
     -.5873179542866171, -.3678314989981802, -.1252334085114692 },
    {0},
    { -.9931285991850949, -.9639719272779138, -.9122344282513259,
      -.8391169718222188, -.7463319064601508, -.636053680726515,
      -.5108670019508271, -.3737060887154196, -.2277858511416451,
      -.07652652113349733 }
  };

#define x ((double *)&equiv_100)

  /* System generated locals */
  int i__1;
  double ret_val, d__1, d__2, d__3, d__4;

  /* Local variables */
  static double a, b, c, d, h;
  static int i;
  static double k;
  static int lg;
  static double as;
  static int ng;
  static double bs, hk, hs, sn, rs, xs;
  static double bvn, asr;

  if (ABS(*r) < .3) {
    ng = 1;
    lg = 3;
  } else if (ABS(*r) < .75) {
    ng = 2;
    lg = 6;
  } else {
    ng = 3;
    lg = 10;
  }
  h = *sh;
  k = *sk;
  hk = h * k;
  bvn = 0.;
  if (ABS(*r) < .925) {
    hs = (h * h + k * k) / 2;
    asr = asin(*r);
    i__1 = lg;
    for (i = 1; i <= i__1; ++i) {
      sn = sin(asr * (x[i + ng * 10 - 11] + 1) / 2);
      bvn += w[i + ng * 10 - 11] * exp((sn * hk - hs) / (1 - sn * sn));
      sn = sin(asr * (-x[i + ng * 10 - 11] + 1) / 2);
      bvn += w[i + ng * 10 - 11] * exp((sn * hk - hs) / (1 - sn * sn));
    }
    d__1 = -h;
    d__2 = -k;
    bvn = bvn * asr / 12.566370614359172 + st_mvnphi(&d__1) * st_mvnphi(&d__2)
      ;
  } else {
    if (*r < 0.) {
      k = -k;
      hk = -hk;
    }
    if (ABS(*r) < 1.) {
      as = (1 - *r) * (*r + 1);
      a = sqrt(as);
      /* Computing 2nd power */
      d__1 = h - k;
      bs = d__1 * d__1;
      c = (4 - hk) / 8;
      d = (12 - hk) / 16;
      bvn = a * exp(-(bs / as + hk) / 2) * (1 - c * (bs - as) * (1 - d * bs / 5) / 3 + c * d * as * as / 5);
      if (hk > -160.) {
        b = sqrt(bs);
        d__1 = -b / a;
        bvn -= exp(-hk / 2) * sqrt(6.283185307179586) * st_mvnphi(&d__1)
          * b * (1 - c * bs * (1 - d * bs / 5) / 3);
      }
      a /= 2;
      i__1 = lg;
      for (i = 1; i <= i__1; ++i) {
        /* Computing 2nd power */
        d__1 = a * (x[i + ng * 10 - 11] + 1);
        xs = d__1 * d__1;
        rs = sqrt(1 - xs);
        bvn += a * w[i + ng * 10 - 11] * (exp(-bs / (xs * 2) - hk / (
                                                rs + 1)) / rs - exp(-(bs / xs + hk) / 2) * (c * xs * (
                                                                                              d * xs + 1) + 1));
        /* Computing 2nd power */
        d__1 = -x[i + ng * 10 - 11] + 1;
        xs = as * (d__1 * d__1) / 4;
        rs = sqrt(1 - xs);
        bvn += a * w[i + ng * 10 - 11] * exp(-(bs / xs + hk) / 2) * (
          exp(-hk * (1 - rs) / ((rs + 1) * 2)) / rs - (c * xs *
                                                       (d * xs + 1) + 1));
      }
      bvn = -bvn / 6.283185307179586;
    }
    if (*r > 0.) {
      d__1 = -MAX(h,k);
      bvn += st_mvnphi(&d__1);
    }
    if (*r < 0.) {
      /* Computing MAX */
      d__3 = -h;
      d__4 = -k;
      d__1 = 0., d__2 = st_mvnphi(&d__3) - st_mvnphi(&d__4);
      bvn = -bvn + MAX(d__1,d__2);
    }
  }
  ret_val = bvn;
  return ret_val;
}
#undef x
#undef w

/*****************************************************************************/
/*!
** A function for computing bivariate normal probabilities
**
** \return bivariate normal probability
**
** \param[in] lower array of lower integration limits
** \param[in] upper array of upper integration limits
** \param[in] infin array of integration limits flag (0,1,2)
** \param[in] correl correlation coefficient
**
*****************************************************************************/
static double st_bvnmvn(double *lower,
                        double *upper,
                        int *infin,
                        double *correl)
{
  double ret_val=0., d__1, d__2, d__3, d__4;

  /* Parameter adjustments */
  --infin;
  --upper;
  --lower;

  /* Function Body */
  if (infin[1] == 2 && infin[2] == 2) {
    ret_val = (st_bvu(&lower[1], &lower[2], correl) -
               st_bvu(&upper[1], &lower[2], correl) -
               st_bvu(&lower[1], &upper[2], correl) +
               st_bvu(&upper[1], &upper[2], correl));
  } else if (infin[1] == 2 && infin[2] == 1) {
    ret_val = (st_bvu(&lower[1], &lower[2], correl) -
               st_bvu(&upper[1], &lower[2], correl));
  } else if (infin[1] == 1 && infin[2] == 2) {
    ret_val = (st_bvu(&lower[1], &lower[2], correl) -
               st_bvu(&lower[1], &upper[2], correl));
  } else if (infin[1] == 2 && infin[2] == 0) {
    d__1 = -upper[1];
    d__2 = -upper[2];
    d__3 = -lower[1];
    d__4 = -upper[2];
    ret_val = st_bvu(&d__1, &d__2, correl) - st_bvu(&d__3, &d__4, correl);
  } else if (infin[1] == 0 && infin[2] == 2) {
    d__1 = -upper[1];
    d__2 = -upper[2];
    d__3 = -upper[1];
    d__4 = -lower[2];
    ret_val = st_bvu(&d__1, &d__2, correl) - st_bvu(&d__3, &d__4, correl);
  } else if (infin[1] == 1 && infin[2] == 0) {
    d__1 = -upper[2];
    d__2 = -(*correl);
    ret_val = st_bvu(&lower[1], &d__1, &d__2);
  } else if (infin[1] == 0 && infin[2] == 1) {
    d__1 = -upper[1];
    d__2 = -(*correl);
    ret_val = st_bvu(&d__1, &lower[2], &d__2);
  } else if (infin[1] == 1 && infin[2] == 1) {
    ret_val = st_bvu(&lower[1], &lower[2], correl);
  } else if (infin[1] == 0 && infin[2] == 0) {
    d__1 = -upper[1];
    d__2 = -upper[2];
    ret_val = st_bvu(&d__1, &d__2, correl);
  }
  return ret_val;
}

/*****************************************************************************/
/*!
** Integrand subroutine
**
*****************************************************************************/
static double st_mvndfn_0(int n__,
                          int *n,
                          double *w,
                          double *correl,
                          double *lower,
                          double *upper,
                          int *infin,
                          int *infis,
                          double *d,
                          double *e)
{
  int i__1, i__2;
  double ret_val, d__1, d__2;


  /* Local variables */
  static int infa, infb, infi[100];
  static double a[100], b[100];
  static int i, j;
  static double y[100], ai, bi, di, ei;
  static int ij, ik;
  static double cov[5050], sum;


  /* Parameter adjustments */
  if (w) {
    --w;
  }
  if (correl) {
    --correl;
  }
  if (lower) {
    --lower;
  }
  if (upper) {
    --upper;
  }
  if (infin) {
    --infin;
  }

  /* Function Body */
  switch(n__) {
    case 1: goto L_mvndnt;
  }

  ret_val = 1.;
  infa = 0;
  infb = 0;
  ik = 1;
  ij = 0;
  i__1 = *n + 1;
  for (i = 1; i <= i__1; ++i) {
    sum = 0.;
    i__2 = i - 1;
    for (j = 1; j <= i__2; ++j) {
      ++ij;
      if (j < ik) {
        sum += cov[ij - 1] * y[j - 1];
      }
    }
    if (infi[i - 1] != 0) {
      if (infa == 1) {
        /* Computing MAX */
        d__1 = ai, d__2 = a[i - 1] - sum;
        ai = MAX(d__1,d__2);
      } else {
        ai = a[i - 1] - sum;
        infa = 1;
      }
    }
    if (infi[i - 1] != 1) {
      if (infb == 1) {
        /* Computing MIN */
        d__1 = bi, d__2 = b[i - 1] - sum;
        bi = MIN(d__1,d__2);
      } else {
        bi = b[i - 1] - sum;
        infb = 1;
      }
    }
    ++ij;
    if (i == *n + 1 || cov[ij + ik] > 0.) {
      i__2 = (infa << 1) + infb - 1;
      st_mvnlms(&ai, &bi, &i__2, &di, &ei);
      if (di >= ei) {
        ret_val = 0.;
        return ret_val;
      } else {
        ret_val *= ei - di;
        if (i <= *n) {
          d__1 = di + w[ik] * (ei - di);
          y[ik - 1] = st_phinvs(&d__1);
        }
        ++ik;
        infa = 0;
        infb = 0;
      }
    }
  }
  return ret_val;

  /*     Entry point for intialization. */


L_mvndnt:
  ret_val = 0.;

  /*     Initialization and computation of covariance Cholesky factor. */

  st_covsrt(n, &lower[1], &upper[1], &correl[1], &infin[1], y, infis, a, b,
            cov, infi);
  if (*n - *infis == 1) {
    st_mvnlms(a, b, infi, d, e);
  } else if (*n - *infis == 2) {
    /* Computing 2nd power */
    d__1 = cov[1];
    *d = sqrt(d__1 * d__1 + 1);
    if (infi[1] != 0) {
      a[1] /= *d;
    }
    if (infi[1] != 1) {
      b[1] /= *d;
    }
    d__1 = cov[1] / *d;
    *e = st_bvnmvn(a, b, infi, &d__1);
    *d = 0.;
    ++(*infis);
  }
  return ret_val;
}

/*****************************************************************************/
/*!
** st_mvndnt
**
*****************************************************************************/
static double st_mvndnt(int *n,
                        double *correl,
                        double *lower,
                        double *upper,
                        int *infin,
                        int *infis,
                        double *d,
                        double *e)
{
  return st_mvndfn_0(1, n, (double *)0, correl, lower, upper, infin,
                     infis, d, e);
}

/*****************************************************************************/
/*!
** This subroutine generates a new quasi-random Richtmeyer vector.
**  A reference is
**  "Methods of Numerical Integration", P.J. Davis and P. Rabinowitz,
**  Academic Press, 1984, pp. 482-483.
**
** \param[in] s the number of dimensions
**
** \param[out] quasi a new quasi-random S-vector
**
*****************************************************************************/
static void st_dkrcht(int *s,
                      double *quasi)
{
  /* Initialized data */

  static int olds = 0;
  static int prime[80] = {
    2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,
    59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,
    149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,
    233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,
    331,337,347,349,353,359,367,373,379,383,389,397,401,409
  };

  int i__1;
  double d__1;

  /* Local variables */
  static double psqt[80];
  static int i, n[49], hisum;
  static double rn;

  /* Parameter adjustments */
  --quasi;

  /* Function Body */
  if (*s != olds || *s < 1) {
    olds = *s;
    n[0] = 0;
    hisum = 0;
    i__1 = *s;
    for (i = 1; i <= i__1; ++i) {
      rn = (double) prime[i - 1];
      psqt[i - 1] = sqrt(rn);
    }
  }
  i__1 = hisum;
  for (i = 0; i <= i__1; ++i) {
    ++n[i];
    if (n[i] < 2) {
      goto L10;
    }
    n[i] = 0;
  }
  ++hisum;
  if (hisum > 48) {
    hisum = 0;
  }
  n[hisum] = 1;
L10:
  rn = 0.;
  for (i = hisum; i >= 0; --i) {
    rn = n[i] + rn * 2;
  }
  i__1 = *s;
  for (i = 1; i <= i__1; ++i) {
    d__1 = rn * psqt[i - 1];
    quasi[i] = fmod(d__1, c_b11);
  }
}

/*****************************************************************************/
/*!
** st_dksmrc
**
*****************************************************************************/
static void st_dksmrc(int *ndim,
                      int *klim,
                      double *sumkro,
                      int *prime,
                      double *vk,
                      double (*functn)(int*,double *),
                      double *x)
{
  /* System generated locals */
  int i__1, i__2;
  double d__1;

  /* Local variables */
  static int j, k, nk, jp;
  static double xt;

  /* Parameter adjustments */
  --x;
  --vk;

  /* Function Body */
  *sumkro = 0.;
  nk = MIN(*ndim,*klim);
  i__1 = nk - 1;
  for (j = 1; j <= i__1; ++j) {
    jp = (int) (j + law_uniform(0.,1.) * (nk + 1 - j));
    xt = vk[j];
    vk[j] = vk[jp];
    vk[jp] = xt;
  }
  i__1 = *ndim;
  for (j = 1; j <= i__1; ++j) {
    x[*ndim + j] = law_uniform(0.,1.);
  }
  i__1 = *prime;
  for (k = 1; k <= i__1; ++k) {
    i__2 = nk;
    for (j = 1; j <= i__2; ++j) {
      d__1 = k * vk[j];
      x[j] = fmod(d__1, c_b11);
    }
    if (*ndim > *klim) {
      i__2 = *ndim - *klim;
      st_dkrcht(&i__2, &x[*klim + 1]);
    }
    i__2 = *ndim;
    for (j = 1; j <= i__2; ++j) {
      xt = x[j] + x[*ndim + j];
      if (xt > 1.) {
        xt += -1;
      }
      x[j] = (d__1 = xt * 2 - 1, ABS(d__1));
    }
    *sumkro += ((*functn)(ndim, &x[1]) - *sumkro) / ((k << 1) - 1);
    i__2 = *ndim;
    for (j = 1; j <= i__2; ++j) {
      x[j] = 1 - x[j];
    }
    *sumkro += ((*functn)(ndim, &x[1]) - *sumkro) / (k << 1);
  }
}

/*****************************************************************************/
/*!
** Automatic Multidimensional Integration Subroutine
**        AUTHOR: Alan Genz
**                Department of Mathematics
**                Washington State University
**                Pulman, WA 99164-3113
**                Email: AlanGenz@wsu.edu
**        Last Change: 5/15/98
** KRBVRC computes an approximation to the integral
**     1  1     1
**     I  I ... I       F(X)  dx(NDIM)...dx(2)dx(1)
**     0  0     0
** DKBVRC uses randomized Korobov rules for the first 20 variables.
** The primary references are
**  "Randomization of Number Theoretic Methods for Multiple Integration"
** *    R. Cranley and T.N.L. Patterson, SIAM J Numer Anal, 13, pp. 904-14,
** *  and
**  "Optimal Parameters for Multidimensional Integration",
**   P. Keast, SIAM J Numer Anal, 10, pp.831-838.
** If there are more than 20 variables, the remaining variables are
** integrated using Richtmeyer rules. A reference is
**  "Methods of Numerical Integration", P.J. Davis and P. Rabinowitz,
**   Academic Press, 1984, pp. 482-483.
**
** \param[in] ndim    Number of variables, must exceed 1, but not exceed 40
**                    Integer minimum number of function evaluations allowed.
** \param[in] minvls  must not exceed MAXVLS.  If MINVLS < 0 then the
**                    routine assumes a previous call has been made with
**                   the same integrand and continues that calculation.
** \param[in] maxvls  Integer maximum number of function evaluations allowed.
** \param[in] functn  EXTERNALLY declared user defined function to be integrated.
**                    It must have parameters (NDIM,Z), where Z is a real array
**                    of dimension NDIM.
** \param[in] abseps  Required absolute accuracy.
** \param[in] releps  Required relative accuracy.
**
** \param[out] minvls  Actual number of function evaluations used.
** \param[out] abserr  Estimated absolute accuracy of FINEST.
** \param[out] finest  Estimated value of integral.
** \param[out] inform = 0 for normal status, when
**                        ABSERR <= MAX(ABSEPS, RELEPS*ABS(FINEST))
**                     and
**                        INTVLS <= MAXCLS.
**                     = 1 If MAXVLS was too small to obtain the required
**                     accuracy. In this case a value FINEST is returned with
**                     estimated absolute accuracy ABSERR.
**
*****************************************************************************/
static void st_dkbvrc(int *ndim,
                      int *minvls,
                      int *maxvls,
                      double (*functn)(int*, double*),
                      double *abseps,
                      double *releps,
                      double *abserr,
                      double *finest,
                      int *inform)
{
  /* Initialized data */

  static int p[25] = {
    31,47,73,113,173,263,397,593,907,1361,2053,3079,
    4621,6947,10427,15641,23473,35221,52837,79259,118891,178349,
    267523,401287,601942 };
  static int c[475] = {
    12,13,27,35,64,111,163,
    246,347,505,794,1189,1763,2872,4309,6610,9861,10327,19540,34566,
    31929,40701,103650,165843,130365,9,11,28,27,66,42,154,189,402,220,
    325,888,1018,3233,3758,6977,3647,7582,19926,9579,49367,69087,
    125480,90647,236711,9,17,10,27,28,54,83,242,322,601,960,259,1500,
    1534,4034,1686,4073,7124,11582,12654,10982,77576,59978,59925,
    110235,13,10,11,36,28,118,43,102,418,644,528,1082,432,2941,1963,
    3819,2535,8214,11113,26856,3527,64590,46875,189541,125699,12,15,
    11,22,44,20,82,250,215,612,247,725,1332,2910,730,2314,3430,9600,
    24585,37873,27066,39397,77172,67647,56483,12,15,20,29,44,31,92,
    250,220,160,247,811,2203,393,642,5647,9865,10271,8726,38806,13226,
    33179,83021,74795,93735,12,15,11,29,55,31,150,102,339,206,338,636,
    126,1796,1502,3953,2830,10193,17218,29501,56010,10858,126904,
    68365,234469,12,15,11,20,67,72,59,250,339,206,366,965,2240,919,
    2246,3614,9328,10800,419,17271,18911,38935,14541,167485,60549,12,
    15,28,45,10,17,76,280,339,206,847,497,1719,446,3834,5115,4320,
    9086,4918,3663,40574,43129,56299,143918,1291,12,15,13,5,10,94,76,
    118,337,422,753,497,1284,919,1511,423,5913,2365,4918,10763,20767,
    35468,43636,74912,93937,12,22,13,5,10,14,47,196,218,134,753,1490,
    878,919,1102,423,10365,4409,4918,18955,20767,35468,11655,167289,
    245291,12,15,28,5,10,14,11,118,315,518,236,1490,1983,1117,1102,
    5408,8272,13812,15701,1298,9686,2196,52680,75517,196061,3,15,13,
    21,10,11,11,191,315,134,334,392,266,103,1522,7426,3706,5661,17710,
    26560,47603,61518,88549,8148,258647,3,6,13,21,10,14,100,215,315,
    134,334,1291,266,103,1522,423,6186,9344,4037,17132,47603,61518,
    29804,172106,162489,3,6,13,21,38,14,131,121,315,518,461,508,266,
    103,3427,423,7806,9344,4037,17132,11736,27945,101894,126159,
    176631,12,6,14,21,38,14,116,121,167,652,711,508,266,103,3427,487,
    7806,10362,15808,4753,11736,70975,113675,35867,204895,7,15,14,21,
    10,94,116,49,167,382,652,1291,747,103,3928,6227,7806,9344,11401,
    4753,41601,70975,48040,35867,73353,7,15,14,21,10,10,116,49,167,
    206,381,1291,747,103,915,2660,8610,9344,19398,8713,12888,86478,
    113675,35867,172319,12,9,14,21,10,10,116,49,167,158,381,508,127,
    103,915,6227,2563,8585,25950,18624,32948,86478,34987,121694,28881
  };

  int i__1, i__2;
  double d__1, d__2;

  /* Local variables */
  static int i;
  static double x[200];
  static int klimi;
  static double value;
  static int np=0;
  static double vk[20], difint, finval;
  static double varprd;
  static int sampls=0;
  static double varest=0, varsqr;
  static int intvls;

  *inform = 1;
  intvls = 0;
  klimi = 20;
  if (*minvls >= 0) {
    *finest = 0.;
    varest = 0.;
    sampls = 8;
    for (i = 1; i <= 25; ++i) {
      np = i;
      if (*minvls < (sampls << 1) * p[i - 1]) {
        goto L10;
      }
    }
    /* Computing MAX */
    i__1 = 8, i__2 = *minvls / (p[np - 1] << 1);
    sampls = MAX(i__1,i__2);
  }
L10:
  vk[0] = 1. / p[np - 1];
  i__1 = MIN(*ndim,20);
  for (i = 2; i <= i__1; ++i) {
    /* Computing MIN */
    i__2 = *ndim - 1;
    d__1 = c[np + MIN(i__2,19) * 25 - 26] * vk[i - 2];
    vk[i - 1] = fmod(d__1, c_b11);
  }
  finval = 0.;
  varsqr = 0.;
  i__1 = sampls;
  for (i = 1; i <= i__1; ++i) {
    st_dksmrc(ndim, &klimi, &value, &p[np - 1], vk, functn, x);
    difint = (value - finval) / i;
    finval += difint;
    /* Computing 2nd power */
    d__1 = difint;
    varsqr = (i - 2) * varsqr / i + d__1 * d__1;
  }
  intvls += (sampls << 1) * p[np - 1];
  varprd = varest * varsqr;
  *finest += (finval - *finest) / (varprd + 1);
  if (varsqr > 0.) {
    varest = (varprd + 1) / varsqr;
  }
  *abserr = sqrt(varsqr / (varprd + 1)) * 3;
  /* Computing MAX */
  d__1 = *abseps, d__2 = ABS(*finest) * *releps;
  if (*abserr > MAX(d__1,d__2)) {
    if (np < 25) {
      ++np;
    } else {

      /* Computing MIN */

      i__1 = sampls * 3 / 2, i__2 = (*maxvls - intvls) / (p[np - 1] <<
                                                          1);
      sampls = MIN(i__1,i__2);
      sampls = MAX(8,sampls);
    }
    if (intvls + (sampls << 1) * p[np - 1] <= *maxvls) {
      goto L10;
    }
  } else {
    *inform = 0;
  }
  *minvls = intvls;
}

/*****************************************************************************/
/*!
** st_mvndfn
**
*****************************************************************************/
static double st_mvndfn(int *n,
                        double *w)
{
  return st_mvndfn_0(0, n, w, (double *)0, (double *)0, (double *)
                     0, (int *)0, (int *)0, (double *)0, (double *)0);
}

/****************************************************************************/
/*!
**  Multivariate Normal Probability
**
** \param[in]  n      Number of variables
** \param[in]  lower  Array of lower integration limits
** \param[in]  upper  Array of upper integration limits
** \param[in]  infin  Array of integration limit flags
** \li                 <0 for ]-Infinity; +Infinity[
** \li                 =0 for ]-Infinity; upper]
** \li                 =1 for [lower; +Infinity[
** \li                 =2 for [lower; upper]
** \param[in]  correl Array of correlation coefficients
** \param[in]  maxpts Maximum number of function values allowed
** \param[in]  abseps Absolute error tolerance
** \param[in]  releps Relative error tolerance
**
** \param[out]  error   Estimated absolute error with 90% confidence level
** \param[out]  value   Estimated value for the integral
** \param[out]  inform  Returned code
**
** \remark The array correl must be entered as the non-diagonal upper part
** \remark of the correlation matrix, entered by line
**
** \remark This subroutine uses an algorithm given in the paper:
** \remark "Numerical Computation of Multivariate Normal Probabilities", in
** \remark  J. of Computational and Graphical Stat., 1(1992), pp. 141-149, by
** \remark         Alan Genz
** \remark      Department of Mathematics
** \remark      Washington State University
** \remark      Pullman, WA 99164-3113
** \remark      Email : AlanGenz@wsu.edu
**
*****************************************************************************/
GEOSLIB_API void mvndst(int     n,
                        double *lower,
                        double *upper,
                        int    *infin,
                        double *correl,
                        int     maxpts,
                        double  abseps,
                        double  releps,
                        double *error,
                        double *value,
                        int    *inform)
{
  int seed_memo;
  int i__1;

  /* Local variables */
  static double d, e;
  static int infis;
  int ivls;

  /* Parameter adjustments */
  --correl;
  --infin;
  --upper;
  --lower;

  /* Set the seed for random number generator */
  seed_memo = law_get_random_seed();
  law_set_random_seed(4323151);

  /* Function Body */
  if (n > 100 || n < 1) {
    *inform = 2;
    *value = 0.;
    *error = 1.;
  } else {
    *inform = (int) st_mvndnt(&n, &correl[1], &lower[1], &upper[1], &
                              infin[1], &infis, &d, &e);
    if (n - infis == 0) {
      *value = 1.;
      *error = 0.;
    } else if (n - infis == 1) {
      *value = e - d;
      *error = (float)2e-16;
    } else {

      /*        Call the lattice rule integration subroutine */

      ivls = 0;
      i__1 = n - infis - 1;
      st_dkbvrc(&i__1, &ivls, &maxpts, st_mvndfn, &abseps, &releps,
                error, value, inform);
    }
  }
  law_set_random_seed(seed_memo);
}

/****************************************************************************/
/*!
**  Calculate the quadri-variable gaussian integral
**
** \return  Integration value
**
** \param[in]  lower        Array of lower bounds
** \param[in]  upper        Array of upper bounds
** \param[in]  correl       Correlation matrix (Dimension: 4*4)
** \param[in]  maxpts       Maximum number of function values allowed
** \param[in]  abseps       Absolute error tolerance
** \param[in]  releps       Relative error tolerance
**
** \param[out]  error       Estimated absolute error with 90% confidence level
** \param[out]  value       Estimated value for the integral
** \param[out]  inform      Returned code
**
*****************************************************************************/
GEOSLIB_API void mvndst4(double *lower,
                         double *upper,
                         double *correl,
                         int     maxpts,
                         double  abseps,
                         double  releps,
                         double *error,
                         double *value,
                         int    *inform)
{
  int    i,j,ecr,infin[4];
  double corloc[6];

  /* Initializations */

  for (i=ecr=0; i<4; i++)
  {
    infin[i] = mvndst_infin(lower[i],upper[i]);
    for (j=0; j<i; j++,ecr++) corloc[ecr] = M_R(correl,4,i,j);
  }

  mvndst(4,lower,upper,infin,corloc,maxpts,abseps,releps,
         error,value,inform);

  return;
}

/****************************************************************************/
/*!
**  Calculate the multigaussian integral (non-normalized)
**
** \return  Integral value
**
** \param[in]  lower        Array of lower bounds (Dimension: nvar)
** \param[in]  upper        Array of upper bounds (Dimension: nvar)
** \param[in]  means        Array of means (Dimension: 2)
** \param[in]  correl       Correlation matrix (Dimension: 2*2)
** \param[in]  maxpts       Maximum number of function values allowed
** \param[in]  abseps       Absolute error tolerance
** \param[in]  releps       Relative error tolerance
**
** \param[out]  error       Estimated absolute error with 90% confidence level
** \param[out]  value       Estimated value for the integral
** \param[out]  inform      Returned code
**
*****************************************************************************/
GEOSLIB_API void mvndst2n(double *lower,
                          double *upper,
                          double *means,
                          double *correl,
                          int     maxpts,
                          double  abseps,
                          double  releps,
                          double *error,
                          double *value,
                          int    *inform)
{
  int    i,infin[2];
  double scale,covar,low[2],upp[2];

  /* Initializations */

  for (i=0; i<2; i++)
  {
    low[i] = lower[i];
    upp[i] = upper[i];
    scale  = sqrt(M_R(correl,2,i,i));
    if (IS_GAUSS_DEF(low[i])) low[i] = (low[i] - means[i]) / scale;
    if (IS_GAUSS_DEF(upp[i])) upp[i] = (upp[i] - means[i]) / scale;
    infin[i] = mvndst_infin(low[i],upp[i]);
  }
  covar = correl[1] / sqrt(M_R(correl,2,0,0) * M_R(correl,2,1,1));

  mvndst(2,low,upp,infin,&covar,maxpts,abseps,releps,error,value,inform);
  return;
}

/****************************************************************************/
/*!
**  Set the flags for the bound of numerical integration
**
** \return  Flag for the bound
**
** \param[in]  low Lower integration bound
** \param[in]  sup Upper integration bound
**
*****************************************************************************/
GEOSLIB_API int mvndst_infin(double low,
                             double sup)
{
  if (low == THRESH_INF && sup == THRESH_SUP) return(-1);
  if (low == THRESH_INF) return(0);
  if (sup == THRESH_SUP) return(1);
  return(2);
}

/*****************************************************************************/
/*!
**  Calculation of the logarithm of the gamma function
**
** \return  Logarithm of the gamma function
**
** \param[in]  parameter raw value 
**
*****************************************************************************/
GEOSLIB_API double loggamma(double parameter)

{
  static double cval[2][8] =
    {
      { 4.120843185847770,85.68982062831317,243.175243524421,
        -261.7218583856145,-922.2613728801522,-517.6383498023218,
        -77.41064071332953,-2.20884399721618 },
      { 1.,45.64677187585908,377.8372484823942,951.323597679706,
        846.0755362020782,262.3083470269460,24.43519662506312,
        0.40977929210926}
    };
  double sval[2],x,xe,p,dalgam;
  int m,k;

  x  = parameter;
  xe = floor(x);
  if (x-xe>0.5) xe += 1.;
  m = (int)(ceil(xe)-1);

  xe = x;
  if (m == -1) xe = x+1.;
  if (m >   0) xe = x-m;

  sval[0] = sval[1] = 0.;
  for (k=0; k<8; k++)
  {
    sval[0] = xe*sval[0] + cval[0][k];
    sval[1] = xe*sval[1] + cval[1][k];
  }

  dalgam = (xe-1.)*sval[0]/sval[1];
  if (m <= -1) return(dalgam-log(x));
  if (m ==  0) return(dalgam);
  if (m ==  1) return(dalgam+log(xe));

  if (m < 33)
  {
    p=1;
    for (k=0; k<m; k++) p *= (xe+k);
    return(dalgam+log(p));
  }
  else
  {
    for (k=0; k<m; k++) dalgam += log(xe+k);
    return(dalgam);
  }
}

/****************************************************************************
 **
 ** FUNCTION: st_tableone_manage
 **
 ** PURPOSE: Management of one item of the CTables structure
 **
 ** RETURNS: For Initialization: 1 if actually performed; 0 otherwise
 ** RETURNS: For Freeing: 1 if actually performed; 0 otherwise
 **
 ** IN_ARGS:  ctables : Pointer to the CTables structure
 ** IN_ARGS:  mode    : 1 for allocation; -1 for deallocation
 ** IN_ARGS:  rank    : rank of the CTable structure
 **
 ** OUT_ARGS: nb_used : Number of defined items for this discretization level
 **                     0 if the discretization level has not been used
 ** OUT_ARGS: nb_max  : Number of pixels in the Discretized covariance array
 **
 *****************************************************************************/
static void st_tableone_manage(CTables *ctables,
                               int      mode,
                               int      rank,
                               int     *nb_used,
                               int     *nb_max)
{
  int nelem,size,number;

  // Initializations

  nelem    = NELEM(ctables);
  size     = (ctables->ndim == 2) ? nelem * nelem : nelem * nelem * nelem;
  *nb_used = 0;
  *nb_max  = size;
    
  // Dispatch
  
  if (mode > 0)
  {
    
    // Allocation

    if (ctables->CT[rank] == (CTable *) NULL)
    {
      if (ctables->CT[rank] == NULL)
      {
        ctables->CT[rank] = (CTable *) mem_alloc(sizeof(CTable),1);
        ctables->CT[rank]->res = (double *)
          mem_alloc(sizeof(double) * size,1);
        for (int i=0; i<size; i++) INTRESX(rank,i) = TEST;
        return;
      }
    }
  }
  else
  {
    if (ctables->CT[rank] != (CTable *) NULL)
    {
      number = 0;
      for (int i=0; i<size; i++)
        if (FFFF(INTRESX(rank,i))) number++;
      ctables->CT[rank]->res = (double *)
        mem_free((char *) ctables->CT[rank]->res);
      *nb_used = number;
      ctables->CT[rank] = (CTable *) mem_alloc(sizeof(CTable),1);
      return;
    }
  }
  return;
}

/****************************************************************************
 **
 ** FUNCTION: ct_tableone_covrank
 **
 ** PURPOSE: Returns the rank of the taget Ctable
 **
 ** RETURNS: Rank of the Target Ctable Rank
 **
 ** IN_ARGS:  ctables : Pointer to the CTables structure
 ** IN_ARGS:  cova    : Covariance value
 **
 ** OUT_ARGS: cround  : Round discretized covariance value
 **
 *****************************************************************************/
GEOSLIB_API int ct_tableone_covrank(CTables *ctables,
                                    double   cova,
                                    double  *cround)
{
  double dc,ecart;
  int    iconf,nconf;

  nconf = ctables->nconf;
  ecart = (cova - ctables->cmin);
  dc    = ctables->dc;
  iconf = (int) (0.5 + ecart / dc);

  if (iconf < 0)      iconf = 0;
  if (iconf >= nconf) iconf = nconf-1;

  *cround = COVAL(ctables,iconf);

  return(iconf);
}

/****************************************************************************
 **
 ** FUNCTION: ct_INTRES2
 **
 ** PURPOSE: Get the value of the Discretized Table for given Configuration
 ** PURPOSE: If not available, the integral is calculated on the fly
 ** PURPOSE: Limitated to the 2-D case
 **
 ** IN_ARGS:  ctables : Pointer to the CTables structure
 ** IN_ARGS:  iconf0  : Rank of the Target Ctable Rank
 ** IN_ARGS:  idisc0  : Index along first dimension
 ** IN_ARGS:  jdisc0  : Index along second dimension
 **
 *****************************************************************************/
GEOSLIB_API double ct_INTRES2(CTables *ctables,
                              int      iconf0,
                              int      idisc0,
                              int      jdisc0)
{
  double  lower[3],upper[3],error,value,cova;
  int     infin[3],inform,nelem,iad,nb_used,nb_max;
  static  double abseps = 1.e-8;
  static  double releps = 0.;
  static  int    maxpts = 25000;

  // Check if integral has already been defined 

  if (ctables->CT[iconf0] == (CTable *) NULL) 
    st_tableone_manage(ctables,1,iconf0,&nb_used,&nb_max);
  
  // Dispatch
  
  nelem = NELEM(ctables);
  iad   = idisc0 + nelem * jdisc0;
  cova  = COVAL(ctables,iconf0);

  // Check if the value has already been calculated

  value = INTRESX(iconf0,iad);
  if (! FFFF(value)) return(value);

  // First call to this value, calculate it

  if (! ctables->flag_cumul)
  {

    // Pixelated case

    lower[0]  = ctables->v[idisc0];
    upper[0]  = ctables->v[idisc0+1];
    infin[0]  = 2;
    if (lower[0] == THRESH_INF) infin[0] = 0;
    if (upper[0] == THRESH_SUP) infin[0] = 1;
    
    lower[1]  = ctables->v[jdisc0];
    upper[1]  = ctables->v[jdisc0+1];
    infin[1]  = 2;
    if (lower[1] == THRESH_INF) infin[1] = 0;
    if (upper[1] == THRESH_SUP) infin[1] = 1;
    
    mvndst(2,lower,upper,infin,&cova,maxpts,abseps,releps,
           &error,&value,&inform);
    if (inform) messageAbort("Error in function 'mvndst'");
  }
  else
  {

    // Cumulative case
      
    lower[0] = THRESH_INF;
    lower[1] = THRESH_INF;
    lower[2] = THRESH_INF;
      
    upper[0] = ctables->v[idisc0];
    infin[0] = 0;
    if (upper[0] == THRESH_SUP) infin[0] = 1;
    
    upper[1]  = ctables->v[jdisc0];
    infin[1]  = 0;
    if (upper[1] == THRESH_SUP) infin[1] = 1;
    
    mvndst(2,lower,upper,infin,&cova,maxpts,abseps,releps,
           &error,&value,&inform);
    if (inform) messageAbort("Error in function 'mvndst'");
  }
    
  // Store it int the table

  INTRESX(iconf0,iad) = value;
  return(value);
}
  
/****************************************************************************
 **
 ** FUNCTION: ct_INTRES3
 **
 ** PURPOSE: Get the value of the Discretized Table for given Configuration
 ** PURPOSE: If not available, the integral is calculated on the fly
 ** PURPOSE: Limitated to the 3-D case
 **
 ** IN_ARGS:  ctables : Pointer to the CTables structure
 ** IN_ARGS:  iconf0  : Rank of the Target Ctable Rank
 ** IN_ARGS:  idisc0  : Index along first dimension
 ** IN_ARGS:  jdisc0  : Index along second dimension
 ** IN_ARGS:  kdisc0  : Index along third dimension
 **
 *****************************************************************************/
GEOSLIB_API double ct_INTRES3(CTables *ctables,
                              int      iconf0,
                              int      idisc0,
                              int      jdisc0,
                              int      kdisc0)
{
  double  lower[3],upper[3],error,value,cova;
  int     infin[3],inform,nelem,iad,nb_used,nb_max;
  static  double abseps = 1.e-8;
  static  double releps = 0.;
  static  int    maxpts = 25000;

  // Check if integral has already been defined 

  if (ctables->CT[iconf0] == (CTable *) NULL) 
    st_tableone_manage(ctables,1,iconf0,&nb_used,&nb_max);
  
  // Dispatch
  
  nelem = NELEM(ctables);
  iad   = idisc0 + nelem * (jdisc0 + nelem * kdisc0);
  cova  = COVAL(ctables,iconf0);

  // Check if the value has already been calculated

  value = INTRESX(iconf0,iad);
  if (! FFFF(value)) return(value);

  // First call to this value, calculate it

  if (! ctables->flag_cumul)
  {
    
    // Pixelated case 
    
    lower[0]  = ctables->v[idisc0];
    upper[0]  = ctables->v[idisc0+1];
    infin[0]  = 2;
    if (lower[0] == THRESH_INF) infin[0] = 0;
    if (upper[0] == THRESH_SUP) infin[0] = 1;
      
    lower[1]  = ctables->v[jdisc0];
    upper[1]  = ctables->v[jdisc0+1];
    infin[1]  = 2;
    if (lower[1] == THRESH_INF) infin[1] = 0;
    if (upper[1] == THRESH_SUP) infin[1] = 1;
    
    lower[2]  = ctables->v[kdisc0];
    upper[2]  = ctables->v[kdisc0+1];
    infin[2]  = 2;
    if (lower[2] == THRESH_INF) infin[2] = 0;
    if (upper[2] == THRESH_SUP) infin[2] = 1;
    
    mvndst(3,lower,upper,infin,&cova,maxpts,abseps,releps,
           &error,&value,&inform);
    if (inform) messageAbort("Error in function 'mvndst'");
  }
  else
  {

    // Cumulative case

    upper[0]  = ctables->v[idisc0];
    infin[0]  = 0;
    if (upper[0] == THRESH_SUP) infin[0] = 1;
    
    upper[1]  = ctables->v[jdisc0];
    infin[1]  = 0;
    if (upper[1] == THRESH_SUP) infin[1] = 1;
    
    upper[2]  = ctables->v[kdisc0];
    infin[2]  = 0;
    if (upper[2] == THRESH_SUP) infin[2] = 1;
    
    mvndst(3,lower,upper,infin,&cova,maxpts,abseps,releps,
           &error,&value,&inform);
    if (inform) messageAbort("Error in function 'mvndst'");
  }
  // Store it int the table

  INTRESX(iconf0,iad) = value;
  return(value);
}

/****************************************************************************
 **
 ** FUNCTION:  ct_tables_print
 **
 ** PURPOSE: Print the CTables structure
 **
 ** IN_ARGS:  ctables    : Address to CTables
 ** IN_ARGS:  flag_print : Verbose option (0, 1 or 2)
 **
 *****************************************************************************/
GEOSLIB_API void ct_tables_print(CTables *ctables,
                                 int      flag_print)
{
  int ndisc,ndim,nconf,nelem;

  ndisc = ctables->ndisc;
  ndim  = ctables->ndim;
  nconf = ctables->nconf;
  nelem = NELEM(ctables);

  mestitle(0,"Precalculation of Gaussian integral");
  message("Point-statistics number                     = %d\n",ndim);
  message("Number of Covariance Discretizations steps  = %d\n",nconf);
  message("Lower Bound of Covariance Discretization    = %lf\n",ctables->cmin);
  message("Upper Bound of Covariance Discretization    = %lf\n",ctables->cmax);
  message("Covariance Discretization Interval          = %lf\n",ctables->dc);
  if (ctables->flag_cumul)
    message("Storing the integral from -Inf to the threshold per class\n");
  else
    message("Storing the integral per discretized class\n");
  message("Covariance is discretized between %lf and %lf\n",
          ctables->cmin,ctables->cmax);
  
  message("\n");
  message("Number of Probability Discretizations       = %d\n",ndisc);
  if (ctables->v != (double *) NULL)
    print_matrix("List of Gaussian Thresholds",
                 0,1,1,ctables->ndisc+1,NULL,ctables->v);

  if (flag_print > 0)
  {
    mestitle(2,"List of the configurations already calculated");

    for (int iconf=0; iconf<ctables->nconf; iconf++)
    {
      if (ctables->CT[iconf] == NULL) continue;

      if (flag_print > 0)
        message("- Configuration %d/%d (Cov=%lf)\n",
                iconf+1,ctables->nconf,COVAL(ctables,iconf));

      if (flag_print == 2)
      {
        if (ndim == 2)
          print_matrix(NULL,0,1,nelem,nelem,NULL,&INTRESX(iconf,0));
        else
          print_matrix(NULL,0,1,nelem,nelem*nelem,NULL,&INTRESX(iconf,0));
      }
    }
    message("\n");
  }
  return;
}

/****************************************************************************
 **
 ** FUNCTION:  ct_tables_manage
 **
 ** PURPOSE: Management of the CTables structure
 **
 ** RETURNS: Pointer to the newly allocated CTables structure
 **
 ** IN_ARGS:  mode        : 1 for allocation; -1 for deallocation
 ** IN_ARGS:  verbose     : Verbose flag
 ** IN_ARGS:  flag_cumul  : 1 for storing gauss integral from -inf
 ** IN_ARGS:                0 for storing gauss integral per pixel
 ** IN_ARGS:  ndim        : Space dimension for discretization
 ** IN_ARGS:  nconf       : Number of configurations
 ** IN_ARGS:  ndisc       : Number of Discretization steps
 ** IN_ARGS:  cmin        : Minimum value allowed for the correlation
 ** IN_ARGS:  cmax        : Maximum value allowed for the correlation
 ** IN_ARGS:  ctables_old : Address to CTables to be freed
 **
 *****************************************************************************/
GEOSLIB_API CTables *ct_tables_manage(int      mode,
                                      int      verbose,
                                      int      flag_cumul,
                                      int      ndim,
                                      int      nconf,
                                      int      ndisc,
                                      double   cmin,
                                      double   cmax,
                                      CTables *ctables_old)
{
  CTables *ctables;
  double  *v;
  int      n_used,nb_used,nb_max;

  /* Dispatch */

  if (mode > 0)
  {
    // Allocation

    if (verbose)
      message("Allocating CTables (%dx%d) for %d possible configurations\n",
              ndisc,ndisc,nconf);
    ctables = (CTables *) mem_alloc(sizeof(CTables),1);
    ctables->flag_cumul = flag_cumul;
    ctables->ndim  = ndim;
    ctables->nconf = nconf;
    ctables->ndisc = ndisc;
    ctables->cmin  = cmin;
    ctables->cmax  = cmax;
    ctables->dc    = (ctables->cmax - ctables->cmin) / (double) (nconf-1);
    ctables->dp    = 1. / (double) ndisc;

    ctables->CT = (CTable **) mem_alloc(sizeof(CTable *) * ctables->nconf,1);
    for (int iconf=0; iconf<ctables->nconf; iconf++)
      ctables->CT[iconf] = (CTable *) NULL;

    // Define the array of thresholds

    v = ctables->v = (double *) mem_alloc(sizeof(double) * (ndisc+1),1);
    v[0]     = THRESH_INF;
    v[ndisc] = THRESH_SUP;
    for (int idisc=0; idisc<ndisc; idisc++)
      v[idisc] = law_invcdf_gaussian((double) idisc * ctables->dp);
  }
  else
  {

    // Desallocation

    ctables = ctables_old;
    if (ctables == (CTables *) NULL) return(ctables);
    ctables->v = (double *) mem_free((char *) ctables->v);
    if (verbose)
      message("Freeing CTables from %d configuration(s)\n",ctables->nconf);

    n_used = 0;
    for (int iconf=0; iconf<ctables->nconf; iconf++)
    {
      st_tableone_manage(ctables,-1,iconf,&nb_used,&nb_max);
      if (nb_used > 0) 
      {
        if (verbose)
          message("Configuration %3d - Number of items used: %d/%d\n",
                  iconf+1,nb_used,nb_max);
        n_used++;
      }
    }
    if (verbose)
      message("Total of configurations actually used: %d\n",n_used);
    ctables = (CTables *) mem_free((char *) ctables);
  }
  return(ctables);
}

/****************************************************************************
 **
 ** FUNCTION: st_tableone_getrank
 **
 ** PURPOSE: Get the starting and ending ranks of the integration
 ** PURPOSE: in the array of thresholds
 **
 ** IN_ARGS:  ctables : Pointer to the CTables structure
 ** IN_ARGS:  low     : Lower bound
 ** IN_ARGS:  up      : Upper bound
 **
 ** OUT_ARGS: indmin  : index of the first discretized inside the interval
 ** OUT_ARGS: indmax  : index of the first discretized outside the interval
 **
 ** REMARKS:  The arguments (imin,imax) are returned so that we can easily
 ** REMARKS:  use them in a loop
 **
 *****************************************************************************/
static void st_tableone_getrank(CTables *ctables,
                                double   low,
                                double   up,
                                int     *indmin,
                                int     *indmax)
{
  double v1;
  int    nelem;

  // Initializations 

  nelem = NELEM(ctables);

  *indmin = *indmax = -1;
  
  // Dispatch

  if (! ctables->flag_cumul)
  {

    // Pixelated case

    for (int idisc=0; idisc<nelem-1; idisc++)
    {
      v1 = (ctables->v[idisc] + ctables->v[idisc+1]) / 2.;
      if (v1 < low) continue;
      if (*indmin < 0) *indmin = idisc;
      if (v1 > up)  
      {
        *indmax = idisc;
        break;
      }
    }
    if (*indmax < 0) *indmax = nelem - 1;
  }
  else
  {

    // Cumulative case

    for (int idisc=0; idisc<nelem; idisc++)
    {
      v1 = ctables->v[idisc];
      if (v1 < low) continue;
      if (*indmin < 0) *indmin = MAX(0,idisc - 1);
      if (v1 >= up)  
      {
        *indmax = MIN(nelem-1,idisc);
        break;
      }
    }
    if (*indmax < 0) *indmax = nelem;
  }

  return;
}

/****************************************************************************
 **
 ** FUNCTION: ct_tableone_calculate
 **
 ** PURPOSE: Evaluate the integral over a rectangle from the integral
 ** PURPOSE: values calculated at pixel size
 **
 ** RETURNS: Integral value
 **
 ** IN_ARGS:  ctables : Pointer to the CTables structure
 ** IN_ARGS:  iconf0  : Rank of the relevant internal structure
 ** IN_ARGS:  lows    : Array of lower values
 ** IN_ARGS:  ups     : Array of upper values
 **
 ** REMARKS: The code has been designed for boundaries of facies areas 
 ** REMAKRS: parallel to the main axes
 **
 *****************************************************************************/
GEOSLIB_API double ct_tableone_calculate(CTables *ctables,
                                         int      iconf0,
                                         double  *lows,
                                         double  *ups)
{
  double result;
  int    ndim,i1min,i1max,i2min,i2max,i3min,i3max;

  // Initializations 

  ndim   = ctables->ndim;
  result = 0;

  // Dispatch

  if (! ctables->flag_cumul)
  {

    // Pixelated case

    if (ndim == 2)
    {

      // 2-D case
      
      st_tableone_getrank(ctables,lows[0],ups[0],&i1min,&i1max);
      st_tableone_getrank(ctables,lows[1],ups[1],&i2min,&i2max);
      for (int idisc=i1min; idisc<i1max; idisc++)
        for (int jdisc=i2min; jdisc<i2max; jdisc++)
          result += ct_INTRES2(ctables,iconf0,idisc,jdisc);
    }
    else
    {
      st_tableone_getrank(ctables,lows[0],ups[0],&i1min,&i1max);
      st_tableone_getrank(ctables,lows[1],ups[1],&i2min,&i2max);
      st_tableone_getrank(ctables,lows[2],ups[2],&i3min,&i3max);

      // 3-D Case
      
      for (int idisc=i1min; idisc<i1max; idisc++)
        for (int jdisc=i2min; jdisc<i2max; jdisc++)
          for (int kdisc=i3min; kdisc<i3max; kdisc++)
            result += ct_INTRES3(ctables,iconf0,idisc,jdisc,kdisc);
    }
  }
  else
  {
    
    // Cumulative case

    if (ctables->ndim == 2)
    {

      // 2-D Case

      st_tableone_getrank(ctables,lows[0],ups[0],&i1min,&i1max);
      st_tableone_getrank(ctables,lows[1],ups[1],&i2min,&i2max);
      result = (ct_INTRES2(ctables,iconf0,i1max,i2max) -
                ct_INTRES2(ctables,iconf0,i1min,i2max) -
                ct_INTRES2(ctables,iconf0,i1max,i2min) + 
                ct_INTRES2(ctables,iconf0,i1min,i2min));
    }
    else
    {

      // 3-D Case

      st_tableone_getrank(ctables,lows[0],ups[0],&i1min,&i1max);
      st_tableone_getrank(ctables,lows[1],ups[1],&i2min,&i2max);
      st_tableone_getrank(ctables,lows[2],ups[2],&i3min,&i3max);
      message("Not calculated yet\n");
    }
  }
  return(result);
}

/****************************************************************************
 **
 ** FUNCTION: ct_tableone_calculate_by_rank
 **
 ** PURPOSE: Evaluate the integral over a rectangle from the integral
 ** PURPOSE: values calculated at pixel size
 **
 ** RETURNS: Integral value
 **
 ** IN_ARGS:  ctables : Pointer to the CTables structure
 ** IN_ARGS:  iconf0  : Rank of the relevant internal structure
 ** IN_ARGS:  rklows  : Array of indices of lower values
 ** IN_ARGS:  rkups   : Array of indices upper values
 **
 ** REMARKS: The arguments 'rklows' and ;rkups' are double although
 ** REMARKS: they designate ranks.
 **
 *****************************************************************************/
GEOSLIB_API double ct_tableone_calculate_by_rank(CTables *ctables,
                                                 int      iconf0,
                                                 double  *rklows,
                                                 double  *rkups)
{
  double result;
  int    ndim;

  // Initializations 

  ndim   = ctables->ndim;
  result = 0;

  // Dispatch

  if (! ctables->flag_cumul)
  {

    // Pixelated case

    if (ndim == 2)
    {

      // 2-D case
      
      for (int idisc=(int) rklows[0]; idisc<(int) rkups[0]; idisc++)
        for (int jdisc=(int) rklows[1]; jdisc<(int) rkups[1]; jdisc++)
          result += ct_INTRES2(ctables,iconf0,idisc,jdisc);
    }
    else
    {

      // 3-D Case
      
      for (int idisc=(int) rklows[0]; idisc<(int) rkups[0]; idisc++)
        for (int jdisc=(int) rklows[1]; jdisc<(int) rkups[1]; jdisc++)
          for (int kdisc=(int) rklows[2]; kdisc<(int) rkups[3]; kdisc++)
            result += ct_INTRES3(ctables,iconf0,idisc,jdisc,kdisc);
    }
  }
  else
  {
    
    // Cumulative case

    if (ctables->ndim == 2)
    {

      // 2-D Case

      result = (ct_INTRES2(ctables,iconf0,(int) rkups[0] ,(int) rkups[1]) -
                ct_INTRES2(ctables,iconf0,(int) rklows[0],(int) rkups[1]) -
                ct_INTRES2(ctables,iconf0,(int) rkups[0] ,(int) rklows[1]) + 
                ct_INTRES2(ctables,iconf0,(int) rklows[0],(int) rklows[1]));
    }
    else
    {

      // 3-D Case

      message("Not calculated yet\n");
    }
  }
  return(result);
}

/****************************************************************************
 **
 ** FUNCTION: ct_tableone_getrank_from_proba
 **
 ** PURPOSE: Get the starting and ending ranks of the integration
 ** PURPOSE: in the array of thresholds (starting from proportions)
 **
 ** RETURNS: Index for the bound in the discretization table
 **
 ** IN_ARGS:  ctables  : Pointer to the CTables structure
 ** IN_ARGS:  gaussian : Gaussian bound
 **
 *****************************************************************************/
GEOSLIB_API int ct_tableone_getrank_from_proba(CTables *ctables,
                                               double   gaussian)
{
  double dp,proba,vmin,vmax;
  int    iad,nelem;

  // Initializations 

  nelem = NELEM(ctables);
  dp    = ctables->dp;

  proba = law_cdf_gaussian(gaussian);

  iad   = proba / dp;
  vmin  = dp * iad;
  vmax  = dp * (iad+1);
  if (vmax - proba < proba - vmin) iad++;

  iad = MIN(nelem,iad);
  
  return(iad);
}
