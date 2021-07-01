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

/*****************************************************************************/
/*!
**   This routine calculates Bessel functions J SUB(NB+ALPHA) (X) for
**   non-negative argument X, and non-negative order NB+ALPHA. 
**
** \return  Error return code : NCALC 
** \return  NCALC < -1:  An argument is out of range.
** \return  1 < NCALC < NB: Not all requested function values could be
** \return  calculated accurately. BY(I) contains correct function values 
**
** \param[in]  x      Working precision non-negative real argument for which 
**                    J's are to be calculated. 
** \param[in]  alpha  Working precision fractional part of order for which J's 
**                    are to be calculated.  0 <= ALPHA < 1.0.
** \param[in]  nb     Integer number of functions to be calculated, NB > 0
**                    The first function calculated is of order ALPHA, 
**                    and the last is of order (NB - 1 + ALPHA). 
**
** \param[out] b      Working precision output vector of length NB. If the 
**                    routine terminates normally (NCALC=NB), the vector by[] 
**                    contains the functions Y(ALPHA,X), ... ,Y(NB-1+ALPHA,X), 
**
** \remark  This program is based on a program written by David J. Sookne (2)
** \remark  that computes values of the Bessel functions J or I of real
** \remark  argument and integer order.  Modifications include the
** \remark  restriction of the computation to the J Bessel function of
** \remark  non-negative real argument, the extension of the computation to
** \remark  arbitrary positive order, and the elimination of most underflow.
** \remark  References: "A Note on Backward Recurrence Algorithms," Olver, F.
** \remark  W. J., and Sookne, D. J., Math. Comp. 26, 1972, pp 941-947.
** \remark  "Bessel Functions of Real Argument and Integer Order," Sookne, D.
** \remark  J., NBS Jour. of Res. B. 77B, 1973, pp 125-132. 
**
*****************************************************************************/
GEOSLIB_API int bessel_j(double  x,
                         double  alpha,
                         int     nb,
                         double *b)
{
  static double enten = 1e38;
  static double ensig = 1e17;
  static double rtnsig = 1e-4;
  static double enmten = 1.2e-37;
  static double xlarge = 1e4;
  static double fact[25] = {
    1.,1.,2.,6.,24.,120.,720.,5040.,40320.,
    362880.,3628800.,39916800.,479001600.,6227020800.,87178291200.,
    1.307674368e12,2.0922789888e13,3.55687428096e14,6.402373705728e15,
    1.21645100408832e17,2.43290200817664e18,5.109094217170944e19,
    1.12400072777760768e21,2.585201673888497664e22,
    6.2044840173323943936e23 };
  static double twopi1 = 6.28125;
  static double twopi2 = .001935307179586476925286767;

  static double capp, capq, pold, gnu, xin, sum, vcos, test, vsin;
  static double p, s, t, z, alpem, halfx, tempa, tempb, tempc, psave;
  static double plast, tover, t1, alp2em, em, en, xc, xk, xm, psavel;
  static int i, j, k, l, m, n, nstart, nbmx, nend, magx;
  int ncalc;

  /* Start */

  --b;
  magx = (int) (x);
  if (nb > 0 && x >= 0. && x <= xlarge && alpha >= 0. && alpha < 1.)
  {

    /* Initialize result array to 0.. */

    ncalc = nb;
    for (i = 1; i <= nb; ++i) b[i] = 0.;

    /* Branch to use 2-term ascending series for small X and asymptotic */
    /* form for large X when NB is not too large. */

    if (x < rtnsig)
    {

      /* Two-term ascending series for small X. */

      tempa = 1.;
      alpem = 1. + alpha;
      halfx = 0.;
      if (x > enmten) halfx = x / 2.;
      if (alpha != 0.)
        tempa = pow(halfx, alpha) /
          (alpha * exp(loggamma(alpha)));
      tempb = 0.;
      if (x > 0.) tempb = -halfx * halfx;
      b[1] = tempa + tempa * tempb / alpem;
      if (x != 0. && b[1] == 0.) ncalc = 0;
      if (nb != 1)
      {
        if (x <= 0.)
          for (n = 2; n <= nb; ++n) b[n] = 0.;
        else
        {

          /* Calculate higher order functions. */

          tempc = halfx;
          tover = (enmten + enmten) / x;
          if (tempb != 0.) tover = enmten / tempb;
          for (n = 2; n <= nb; ++n)
          {
            tempa /= alpem;
            alpem += 1.;
            tempa *= tempc;
            if (tempa <= tover * alpem) tempa = 0.;
            b[n] = tempa + tempa * tempb / alpem;
            if (b[n] == 0. && ncalc > n) ncalc = n - 1;
          }
        }
      }
    }
    else if (x > 25. && nb <= magx + 1)
    {

      /* Asymptotic series for X .GT. 21.0. */

      xc = sqrt(2. / GV_PI / x);
      xin = 1. / (64. * x * x);
      m = 11;
      if (x >=  35.) m = 8;
      if (x >= 130.) m = 4;
      xm = 4. * (double) m;

      /* Argument reduction for SIN and COS routines. */

      t = floor(x / (twopi1 + twopi2) + 0.5);
      z = x - t * twopi1 - t * twopi2 - (alpha + 0.5) * GV_PI / 2.;
      vsin = sin(z);
      vcos = cos(z);
      gnu = alpha + alpha;
      for (i = 1; i <= 2; ++i)
      {
        s = (xm - 1. - gnu) * (xm - 1. + gnu) * xin / 2.;
        t = (gnu - (xm - 3.)) * (gnu + (xm - 3.));
        capp = s * t / fact[m * 2];
        t1 = (gnu - (xm + 1.)) * (gnu + (xm + 1.));
        capq = s * t1 / fact[(m << 1) + 1];
        xk = xm;
        k = m + m;
        t1 = t;
        for (j = 2; j <= m; ++j)
        {
          xk -= 4.;
          s = (xk - 1. - gnu) * (xk - 1. + gnu);
          t = (gnu - (xk - 3.)) * (gnu + (xk - 3.));
          capp = (capp + 1. / fact[k - 2]) * s * t * xin;
          capq = (capq + 1. / fact[k - 1]) * s * t1 * xin;
          k += -2;
          t1 = t;
        }
        capp += 1.;
        capq = (capq + 1.) * (gnu * gnu - 1.) / (8. *x);
        b[i] = xc * (capp * vcos - capq * vsin);
        if (nb == 1) return(ncalc);
        t = vsin;
        vsin = -vcos;
        vcos = t;
        gnu += 2.;
      }

      /* If  NB .GT. 2, compute J(X,ORDER+I)  I = 2, NB-1 */

      if (nb > 2)
      {
        gnu = alpha + alpha + 2.;
        for (j = 3; j <= nb; ++j)
        {
          b[j] = gnu * b[j - 1] / x - b[j - 2];
          gnu += 2.;
        }
      }

      /* Use recurrence to generate results.  First initialize the */
      /* calculation of P*S. */

    }
    else
    {
      nbmx = nb - magx;
      n = magx + 1;
      en = 2. * n + (alpha + alpha);
      plast = 1.;
      p = en / x;

      /* Calculate general significance test. */

      test = ensig + ensig;
      if (nbmx >= 3)
      {

        /* Calculate P*S until N = NB-1.  Check for overflow. */

        tover = enten / ensig;
        nstart = magx + 2;
        nend = nb - 1;
        en = 2. * nstart - 2. + (alpha + alpha);
        for (k = nstart; k <= nend; ++k)
        {
          n = k;
          en += 2.;
          pold = plast;
          plast = p;
          p = en * plast / x - pold;
          if (p > tover)
          {

            /* To avoid overflow, divide P*S by TOVER. */
            /* Calculate P*S until ABS(P) .GT. 1. */

            tover = enten;
            p /= tover;
            plast /= tover;
            psave = p;
            psavel = plast;
            nstart = n + 1;
          L100:
            ++n;
            en += 2.;
            pold = plast;
            plast = p;
            p = en * plast / x - pold;
            if (p <= 1.) goto L100;
            tempb = en / x;

            /* Calculate backward test and find NCALC */
            /* the highest N such that the test is passed. */

            test = pold * plast * (0.5 - 0.5 / (tempb * tempb));
            test /= ensig;
            p = plast * tover;
            --n;
            en -= 2.;
            nend = MIN(nb,n);
            for (l = nstart; l <= nend; ++l)
            {
              pold = psavel;
              psavel = psave;
              psave = en * psavel / x - pold;
              if (psave * psavel > test)
              {
                ncalc = l - 1;
                goto L190;
              }
            }
            ncalc = nend;
            goto L190;
          }
        }
        n = nend;
        en = 2. * n + (alpha + alpha);

        /* Calculate special significance test for NBMX .GT. 2. */

        test = MAX(test,sqrt(plast * ensig) * sqrt(p + p));
      }

      /* Calculate P*S until significance test passes. */

    L140:
      ++n;
      en += 2.;
      pold = plast;
      plast = p;
      p = en * plast / x - pold;
      if (p < test) goto L140;

      /* Initialize the backward recursion and the normalization sum. */

    L190:
      ++n;
      en += 2.;
      tempb = 0.;
      tempa = 1. / p;
      m = (n << 1) - (n / 2 << 2);
      sum = 0.;
      em = floor(n / 2.);
      alpem = em - 1. + alpha;
      alp2em = em + em + alpha;
      if (m != 0) sum = tempa * alpem * alp2em / em;
      nend = n - nb;
      if (nend > 0)
      {

        /* Recur backward via difference equation, calculating */
        /* but not storing) B(N), until N = NB. */

        for (l = 1; l <= nend; ++l)
        {
          --n;
          en -= 2.;
          tempc = tempb;
          tempb = tempa;
          tempa = en * tempb / x - tempc;
          m = 2 - m;
          if (m != 0)
          {
            em -= 1.;
            alp2em = em + em + alpha;
            if (n == 1) goto L210;
            alpem = em - 1. + alpha;
            if (alpem == 0.) alpem = 1.;
            sum = (sum + tempa * alp2em) * alpem / em;
          }
        }
      }

      /* Store B(NB). */

    L210:
      b[n] = tempa;
      if (nend >= 0)
      {
        if (nb <= 1)
        {
          alp2em = alpha;
          if (alpha + 1. == 1.) alp2em = 1.;
          sum += b[1] * alp2em;
          goto L250;
        }
        else
        {

          /* Calculate and store B(NB-1). */

          --n;
          en -= 2.;
          b[n] = en * tempa / x - tempb;
          if (n == 1) goto L240;
          m = 2 - m;
          if (m != 0)
          {
            em -= 1.;
            alp2em = em + em + alpha;
            alpem = em - 1. + alpha;
            if (alpem == 0.) alpem = 1.;
            sum = (sum + b[n] * alp2em) * alpem / em;
          }
        }
      }
      nend = n - 2;
      if (nend != 0)
      {

        /* Calculate via difference equation and store B(N) */
        /* until N = 2. */

        for (l = 1; l <= nend; ++l)
        {
          --n;
          en -= 2.;
          b[n] = en * b[n + 1] / x - b[n + 2];
          m = 2 - m;
          if (m != 0)
          {
            em -= 1.;
            alp2em = em + em + alpha;
            alpem = em - 1. + alpha;
            if (alpem == 0.) alpem = 1.;
            sum = (sum + b[n] * alp2em) * alpem / em;
          }
        }
      }

      /* Calculate B(1). */

      b[1] = 2. * (alpha + 1.) * b[2] / x - b[3];
    L240:
      em -= 1.;
      alp2em = em + em + alpha;
      if (alp2em == 0.) alp2em = 1.;
      sum += b[1] * alp2em;

      /* Normalize.  Divide all B(N) by sum. */

    L250:
      if (alpha + 1. != 1.)
        sum *= exp(loggamma(alpha)) * pow(x/2., -alpha);
      tempa = enmten;
      if (sum > 1.) tempa *= sum;
      for (n = 1; n <= nb; ++n)
      {
        if (ABS(b[n]) < tempa) b[n] = 0.;
        b[n] /= sum;
      }
    }

    /* Error return -- X, NB, or ALPHA is out of range. */

  }
  else
  {
    b[1] = 0.;
    ncalc = MIN(nb,0) - 1;
  }

  return(ncalc);
}

/*****************************************************************************/
/*!
**   This routine calculates modified Bessel functions of the second
**   kind, K SUB(N+ALPHA) (X), for non-negative argument X and
**   non-negative order N+ALPHA 
**
** \return  Error return code : NCALC NCALC < -1:  An argument is out of
** \return  range. 0 < NCALC < NB: Not all requested function values could be
** \return  calculated accurately.  BK(I) contains correct function values
** \return  for I <= NCALC, and contains the ratios
** \return  K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array. 
**
** \param[in]  x      Working precision non-negative real argument for which 
**                    K's are to calculated. If K's are to be calculated, 
**                    X must not be greater than XMAX. 
** \param[in]  alpha  Working precision fractional part of order for which K's 
**                    are to be calculated.  0 <= ALPHA <1.0. 
** \param[in]  nb     Integer number of functions to be calculated, NB > 0. 
**                    The first function calculated is of order ALPHA, and the
**                    last is of order (NB - 1 + ALPHA). 
**
** \param[out] bk     Working precision output vector of length NB. 
**                    If the routine terminates normally (NCALC=NB), the vector
**                    BK contains the functions : 
**                    K(ALPHA,X), ... , K(NB-1+ALPHA,X), 
**
** \remark  This program is based on a program written by J. B. Campbell (2)
** \remark  that computes values of the Bessel functions K of real argument
** \remark  and real order. References: "On Temme's Algorithm for the
** \remark  Modified Bessel Functions of the Third Kind," Campbell, J. B.,
** \remark  TOMS 6(4), Dec. 1980, pp. 581-586. "A FORTRAN IV Subroutine for
** \remark  the Modified Bessel Functions of the Third Kind of Real Order and
** \remark  Real Argument," Campbell, J. B., Report NRC/ERB-925, National
** \remark  Research Council, Canada. 
**
*****************************************************************************/
GEOSLIB_API int bessel_k(double  x,
                         double  alpha,
                         int     nb,
                         double *bk)
{
  static double p[] = {
    .805629875690432845,20.4045500205365151,
    157.705605106676174,536.671116469207504,900.382759291288778,
    730.923886650660393,229.299301509425145,.822467033424113231 };
  static double q[] = {
    29.4601986247850434,277.577868510221208,
    1206.70325591027438,2762.91444159791519,3443.74050506564618,
    2210.63190113378647,572.267338359892221 };
  static double r[] = {
    -.48672575865218401848,13.079485869097804016,
    -101.96490580880537526,347.65409106507813131,
    3.495898124521934782e-4 };
  static double s[] = {
    -25.579105509976461286,212.57260432226544008,
    -610.69018684944109624,422.69668805777760407 };
  static double t[] = {
    1.6125990452916363814e-10,
    2.5051878502858255354e-8,2.7557319615147964774e-6,
    1.9841269840928373686e-4,.0083333333333334751799,
    .16666666666666666446 };
  static double estm[] = {
    52.0583,5.7607,2.7782,14.4303,185.3004,9.3715 };
  static double estf[] = {
    41.8341,7.1075,6.4306,42.511,1.35633,84.5096,20. };
  static double sqxmin = 1.49e-154;
  static double xinf   = 1.79e308;
  static double xmin   = 2.23e-308;
  static double xmax   = 705.342;
  static double tinyx  = 1e-10;
  static double a      = .11593151565841244881;
  static double d      = .797884560802865364;
  static double eps    = 2.22e-16;

  static double x2by4, twox, c, blpha, dm, ex, bk1, bk2, enu;
  static double ratio, wminf, d1, d2, d3, f0, f1, f2, p0, q0, t1, t2, twonu;
  static int i, j, k, m, iend, itemp, mplus1, ncalc;

  /* Parameter adjustments */

  --bk;
  ex = x;
  enu = alpha;
  ncalc = MIN(nb,0) - 2;
  if (nb > 0 && (enu >= 0. && enu < 1.) && (ex <= xmax) && ex > 0.)
  {
    k = 0;
    if (enu < sqxmin) enu = 0.;
    if (enu > 0.5)
    {
      k = 1;
      enu -= 1.;
    }
    twonu = enu + enu;
    iend = nb + k - 1;
    c = enu * enu;
    d3 = -c;
    if (ex <= 1.)
    {

      /*  Calculation of P0 = GAMMA(1+ALPHA) * (2/X)**ALPHA */
      /*                 Q0 = GAMMA(1-ALPHA) * (X/2)**ALPHA */

      d1 = 0.;
      d2 = p[0];
      t1 = 1.;
      t2 = q[0];
      for (i = 2; i <= 7; i += 2)
      {
        d1 = c * d1 + p[i - 1];
        d2 = c * d2 + p[i];
        t1 = c * t1 + q[i - 1];
        t2 = c * t2 + q[i];
      }
      d1 = enu * d1;
      t1 = enu * t1;
      f1 = log(ex);
      f0 = a + enu * (p[7] - enu * (d1 + d2) / (t1 + t2)) - f1;
      q0 = exp(-enu * (a - enu * (p[7] + enu * (d1 - d2) / (t1 - t2)) - f1));
      f1 = enu * f0;
      p0 = exp(f1);

      /*  Calculation of F0 = */

      d1 = r[4];
      t1 = 1.;
      for (i = 1; i <= 4; ++i)
      {
        d1 = c * d1 + r[i - 1];
        t1 = c * t1 + s[i - 1];
      }
      if (ABS(f1) <= 0.5)
      {
        f1 *= f1;
        d2 = 0.;
        for (i = 1; i <= 6; ++i) d2 = f1 * d2 + t[i - 1];
        d2 = f0 + f0 * f1 * d2;
      }
      else
        d2 = sinh(f1) / enu;
      f0 = d2 - enu * d1 / (t1 * p0);
      if (ex <= tinyx)
      {

        /*  X <= 1.0E-10 */
        /*  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X) */

        bk[1] = f0 + ex * f0;
        bk[1] -= ex * bk[1];
        ratio = p0 / f0;
        c = ex * xinf;
        if (k != 0)
        {

          /*  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X), */
          /*  ALPHA >= 1/2 */

          ncalc = -1;
          if (bk[1] >= c / ratio) return(ncalc);
          bk[1] = ratio * bk[1] / ex;
          twonu += 2.;
          ratio = twonu;
        }
        ncalc = 1;
        if (nb == 1) return(ncalc);

        /*  Calculate  K(ALPHA+L,X)/K(ALPHA+L-1,X),  L  =  1,2,..., NB-1 */

        ncalc = -1;
        for (i = 2; i <= nb; ++i)
        {
          if (ratio >= c) return(ncalc);
          bk[i] = ratio / ex;
          twonu += 2.;
          ratio = twonu;
        }
        ncalc = 1;
        goto L420;
      }
      else
      {

        /*  1.0E-10 < X <= 1.0 */

        c = 1.;
        x2by4 = ex * ex / 4.;
        p0 = 0.5 * p0;
        q0 = 0.5 * q0;
        d1 = -1.;
        d2 = bk1 = bk2 = 0.;
        f1 = f0;
        f2 = p0;

      L100:
        d1 += 2.;
        d2 += 1.;
        d3 = d1 + d3;
        c = x2by4 * c / d2;
        f0 = (d2 * f0 + p0 + q0) / d3;
        p0 /= d2 - enu;
        q0 /= d2 + enu;
        t1 = c * f0;
        t2 = c * (p0 - d2 * f0);
        bk1 += t1;
        bk2 += t2;
        if (ABS(t1 / (f1 + bk1)) > eps ||
            ABS(t2 / (f2 + bk2)) > eps) goto L100;
        bk1 = f1 + bk1;
        bk2 = 2. * (f2 + bk2) / ex;
        wminf = estf[0] * ex + estf[1];
      }
    }
    else if (eps * ex > 1.)
    {

      /*  X > 1./EPS */

      ncalc = nb;
      bk1 = 1. / (d * sqrt(ex));
      for (i = 1; i <= nb; ++i) bk[i] = bk1;
      return(ncalc);
    }
    else
    {

      /*  X > 1.0 */

      twox = ex + ex;
      blpha = ratio = 0.;
      if (ex <= 4.)
      {

        /*  Calculation of K(ALPHA+1,X)/K(ALPHA,X),  1.0 <= X <= 4.0 */

        d2 = floor(estm[0] / ex + estm[1]);
        m = (int) d2;
        d1 = d2 + d2;
        d2 -= 0.5;
        d2 *= d2;
        for (i = 2; i <= m; ++i)
        {
          d1 -= 2.;
          d2 -= d1;
          ratio = (d3 + d2) / (twox + d1 - ratio);
        }

        /*  Calculation of I(|ALPHA|,X) and I(|ALPHA|+1,X) by backward */
        /*    recurrence and K(ALPHA,X) from the wronskian */

        d2 = floor(estm[2] * ex + estm[3]);
        m = (int) d2;
        c = ABS(enu);
        d3 = c + c;
        d1 = d3 - 1.;
        f1 = xmin;
        f0 = (2. * (c + d2) / ex + 0.5 * ex / (c + d2 + 1.)) * xmin;
        for (i = 3; i <= m; ++i)
        {
          d2 -= 1.;
          f2 = (d3 + d2 + d2) * f0;
          blpha = (1. + d1 / d2) * (f2 + blpha);
          f2 = f2 / ex + f1;
          f1 = f0;
          f0 = f2;
        }
        f1 = (d3 + 2.) * f0 / ex + f1;
        d1 = 0.;
        t1 = 1.;
        for (i = 1; i <= 7; ++i)
        {
          d1 = c * d1 + p[i - 1];
          t1 = c * t1 + q[i - 1];
        }
        p0 = exp(c * (a + c * (p[7] - c * d1 / t1) - log(ex))) / ex;
        f2 = (c + 0.5 - ratio) * f1 / ex;
        bk1 = p0 + (d3 * f0 - f2 + f0 + blpha) / (f2 + f1 + f0) * p0;
        bk1 *= exp(-ex);
        wminf = estf[2] * ex + estf[3];
      }
      else
      {

        /*  Calculation of K(ALPHA,X) and K(ALPHA+1,X)/K(ALPHA,X) */
        /*  by backward recurrence, for  X > 4.0 */

        dm = floor(estm[4] / ex + estm[5]);
        m = (int) dm;
        d2 = dm - 0.5;
        d2 *= d2;
        d1 = dm + dm;
        for (i = 2; i <= m; ++i)
        {
          dm -= 1.;
          d1 -= 2.;
          d2 -= d1;
          ratio = (d3 + d2) / (twox + d1 - ratio);
          blpha = (ratio + ratio * blpha) / dm;
        }
        bk1 = 1. / ((d + d * blpha) * sqrt(ex));
        bk1 *= exp(-ex);
        wminf = estf[4]*(ex - ABS(ex - estf[6])) + estf[5];
      }

      /*  Calculation of K(ALPHA+1,X) from K(ALPHA,X) and */
      /*    K(ALPHA+1,X)/K(ALPHA,X) */

      bk2 = bk1 + bk1 * (enu + 0.5 - ratio) / ex;
    }

    /*  Calculation of 'NCALC', K(ALPHA+I,X), I  =  0, 1, ... , NCALC-1, */
    /*  K(ALPHA+I,X)/K(ALPHA+I-1,X), I  =  NCALC, NCALC+1, ... , NB-1 */

    ncalc = nb;
    bk[1] = bk1;
    if (iend == 0) return(ncalc);
    j = 2 - k;
    if (j > 0) bk[j] = bk2;
    if (iend == 1) return(ncalc);

    /* Computing MIN */
    m = MIN((int) (wminf - enu),iend);
    for (i = 2; i <= m; ++i)
    {
      t1 = bk1;
      bk1 = bk2;
      twonu += 2.;
      if (ex < 1.)
      {
        if (bk1 >= xinf / twonu * ex) goto L195;
      }
      else
      {
        if (bk1 / ex >= xinf / twonu) goto L195;
      }

      bk2 = twonu / ex * bk1 + t1;
      itemp = i;
      ++j;
      if (j > 0) bk[j] = bk2;
    }

  L195:
    m = itemp;
    if (m == iend) return(ncalc);
    ratio = bk2 / bk1;
    mplus1 = m + 1;
    ncalc = -1;
    for (i = mplus1; i <= iend; ++i)
    {
      twonu += 2.;
      ratio = twonu / ex + 1. / ratio;
      ++j;
      if (j > 1)
        bk[j] = ratio;
      else
      {
        if (bk2 >= xinf / ratio) return(ncalc);
        bk2 = ratio * bk2;
      }
    }

    /* Computing MAX */
    ncalc = MAX(mplus1 - k,1);
    if (ncalc == 1) bk[1] = bk2;
    if (nb == 1) return(ncalc);

  L420:
    j = ncalc + 1;
    for (i = j; i <= nb; ++i)
    {
      if (bk[ncalc] >= xinf / bk[i]) return(ncalc);
      bk[i] = bk[ncalc] * bk[i];
      ncalc = i;
    }
  }
  return(ncalc);
}
