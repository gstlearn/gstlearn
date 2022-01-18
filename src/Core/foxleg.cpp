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
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "Basic/Utilities.hpp"
#include "Basic/EJustify.hpp"
#include "Basic/File.hpp"
#include "Basic/String.hpp"

#include <math.h>

/*! \cond */
#define TOLVAL     1.e-09

#define IAD(ic,ipar)             ((ipar) + NPAR * (ic))
#define IADAC(ic,iparac)         ((iparac) + NPARAC * (ic))
#define JR(idat,ipar)            (Jr[(ipar) * NDAT + (idat)])
#define GAUSS(ipar,jpar)         (gauss[(ipar) * NPCT + (jpar)])
#define GAUSS_RED(iparac,jparac) (gauss_red[(iparac) * (NPARAC+NCONT) + (jparac)])
#define AI(ic,ipar,jpar)         (ai[(jpar) * NPAR2 + IAD(ic,ipar)])
#define AI_RED(ic,iparac,jparac) (ai_red[(jparac) * NPARAC2 + IADAC(ic,iparac)])
#define BORDS(ic,ipar)           (bords[IAD(ic,ipar)])
#define BORDS_RED(ic,iparac)     (bords_red[IADAC(ic,iparac)])
#define A(i,j)                   (a[(j) * size + (i)])
#define POSSIBLE(ipar)           (ind_util[ipar] > 0)
#define SIGNE(ic)                ((ic == 0) ? +1 : -1)
#define FLAG_ACTIVE(ic,iparac)   (flag_active[IADAC(ic,iparac)])
#define ACONT(ipar,icont)        (acont[(icont) * NPAR + (ipar)])
/*! \endcond */

static int VERBOSE_GQO = 0;

static int NPAR, NPAR2, NPARAC, NPARAC2, NDAT, NCONT, NPCT, NPCT2;
static int ITERATION, SOUSITER;
static void (*FUNC_EVALUATE)(int ndat,
                             int npar,
                             VectorDouble &param,
                             VectorDouble &work);

/****************************************************************************/
/*!
 **  Calculate the gradient
 **
 ** \param[in]  param     Current values of the parameters
 ** \param[in]  lower     Array of lower values
 ** \param[in]  upper     Array of upper values
 ** \param[in]  scale     Array of scaling values
 ** \param[in]  tabwgt    Array of weights
 **
 ** \param[out] Jr        Array of gradients
 ** \param[out] param1    Working array (Dimension: NPAR)
 ** \param[out] param2    Working array (Dimension: NPAR)
 ** \param[out] tabmod1   Working array (Dimension: NDAT)
 ** \param[out] tabmod2   Working array (Dimension: NDAT)
 **
 *****************************************************************************/
static void st_gradient(VectorDouble &param,
                        VectorDouble &lower,
                        VectorDouble &upper,
                        VectorDouble &scale,
                        VectorDouble &tabwgt,
                        VectorDouble &Jr,
                        VectorDouble &param1,
                        VectorDouble &param2,
                        VectorDouble &tabmod1,
                        VectorDouble &tabmod2)
{
  double ratio1, ratio2, weight, top, bot, epsloc, epsgrad;
  int idat, ipar, jpar;

  /* Calculate the gradients */

  epsgrad = EPSILON3;
  for (ipar = 0; ipar < NPAR; ipar++)
  {
    epsloc = ABS(epsgrad * scale[ipar]);
    epsloc = MAX(epsgrad, epsloc);
    for (jpar = 0; jpar < NPAR; jpar++)
      param1[jpar] = param2[jpar] = param[jpar];

    param1[ipar] = param[ipar] + epsloc;
    if (!FFFF(upper[ipar])) param1[ipar] = MIN(upper[ipar], param1[ipar]);
    param2[ipar] = param[ipar] - epsloc;
    if (!FFFF(lower[ipar])) param2[ipar] = MAX(lower[ipar], param2[ipar]);

    ratio1 = ratio2 = epsloc;
    if (!FFFF(upper[ipar])) ratio1 = MIN(epsloc, upper[ipar] - param[ipar]);
    if (!FFFF(lower[ipar])) ratio2 = MIN(epsloc, param[ipar] - lower[ipar]);

    FUNC_EVALUATE(NDAT, NPAR, param1, tabmod1);
    FUNC_EVALUATE(NDAT, NPAR, param2, tabmod2);

    bot = ratio1 + ratio2;
    for (idat = 0; idat < NDAT; idat++)
    {
      top = tabmod1[idat] - tabmod2[idat];
      weight = (!tabwgt.empty()) ? tabwgt[idat] : 1.;
      JR(idat,ipar) = (bot != 0.) ? weight * top / bot : 0.;
    }
  }
  return;
}

/****************************************************************************/
/*!
 **  Calculate the residuals between the model and the experimental
 **
 ** \return  The weighted mean squared error
 **
 ** \param[in]  param     Current values of the parameters
 ** \param[in]  tabexp    Array of values at control points
 ** \param[in]  tabwgt    Array of weights at control points
 **
 ** \param[out] tabmod    Working array (Dimension: NDAT)
 ** \param[out] residuals Array of residuals
 **
 *****************************************************************************/
static double st_residuals(VectorDouble &param,
                           VectorDouble &tabexp,
                           VectorDouble &tabwgt,
                           VectorDouble &tabmod,
                           VectorDouble &residuals)
{
  int idat;
  double value, msse, weight;

  /* Evaluate the Model at conditioning points */

  FUNC_EVALUATE(NDAT, NPAR, param, tabmod);

  /* Evaluate the residuals */

  msse = 0.;
  for (idat = 0; idat < NDAT; idat++)
  {
    weight = (!tabwgt.empty()) ? tabwgt[idat] :
                                 1.;
    value = weight * (tabmod[idat] - tabexp[idat]);
    msse += value * value;
    residuals[idat] = value;
  }

  return (msse / 2.);
}

/****************************************************************************/
/*!
 **  Calculate the Gauss matrix
 **
 ** \param[in]  Jr        Array of gradients
 **
 ** \param[out] gauss     Gaussian matrix
 **
 *****************************************************************************/
static void st_determine_gauss(VectorDouble &Jr, VectorDouble &gauss)
{
  int idat, ipar, jpar;
  double value;
  for (ipar = 0; ipar < NPAR; ipar++)
    for (jpar = 0; jpar < NPAR; jpar++)
    {
      value = 0.;
      for (idat = 0; idat < NDAT; idat++)
        value += JR(idat,ipar) * JR(idat, jpar);
      GAUSS(ipar,jpar) = value;
    }

  return;
}

/****************************************************************************/
/*!
 **  Calculate the Norm of the HGN vector
 **
 ** \return  The norm value
 **
 ** \param[in]  hgn          Working vector
 ** \param[in]  scale        Scaling values
 **
 *****************************************************************************/
static double st_norm_hgn(VectorDouble &hgn, VectorDouble &scale)
{
  double v1, norme;
  int ipar;

  norme = 0.;
  for (ipar = 0; ipar < NPAR; ipar++)
  {
    //    v1 = ABS(hgn[ipar]);
    v1 = ABS(hgn[ipar] / scale[ipar]);
    if (v1 > norme) norme = v1;
  }
  return (norme);
}

/****************************************************************************/
/*!
 **  Score of the Minimization under constraints
 **
 ** \param[in]  hgnadm     Admissible vector
 ** \param[in]  grad_red   Reduced Gradient matrix
 ** \param[in]  gauss_red  Reduced Gauss matrix
 **
 *****************************************************************************/
static double st_essai(VectorDouble &hgnadm,
                       VectorDouble &grad_red,
                       VectorDouble &gauss_red)
{
  double v1, v2, result;

  matrix_product(1, NPARAC, 1, hgnadm.data(), grad_red.data(), &v1);
  v2 = matrix_normA(hgnadm.data(), gauss_red.data(), NPARAC + NCONT, NPARAC);
  result = v1 + v2 / 2.;

  return (result);
}

/****************************************************************************/
/*!
 **  Solve the direct minimization problem
 **
 ** \return Error return code
 **
 ** \param[in] npar         Current number of parameters
 ** \param[in] grad         Gradient matrix
 ** \param[in] gauss        Gauss matrix
 ** \param[in] invhess      Inverse Hessian matrix
 ** \param[in] hgnc         Resulting hgnc array
 ** \param[in] flaginvsign  if 1, the result is multiplied by -1
 **
 *****************************************************************************/
static int st_solve_hgnc(int npar,
                         VectorDouble &grad,
                         VectorDouble &gauss,
                         VectorDouble &invhess,
                         VectorDouble &hgnc,
                         int flaginvsign)
{
  int i, j;
  double value, signe;
  VectorDouble TEMPMATRIX(NPCT2 * NPCT2);
  VectorDouble TEMPVECTOR(NPCT2);

  for (i = 0; i < npar; i++)
    for (j = 0; j < npar; j++)
      TEMPMATRIX[i * npar + j] = gauss[i * npar + j];

  for (i = 0; i < npar; i++)
  {
    value = gauss[i + i * npar];
    value = (value < EPSILON10) ? 1 :
                                  value;
    TEMPVECTOR[i] = grad[i] / sqrt(value);
    for (j = 0; j < npar; j++)
    {
      TEMPMATRIX[j * npar + i] /= sqrt(value);
      TEMPMATRIX[i * npar + j] /= sqrt(value);
    }
  }

  if (matrix_invgen(TEMPMATRIX.data(), npar, invhess.data(), NULL))
  {
    messerr("Error: Singularity in the Generalized Inverse");
    messerr("The Automatic Fitting Procedure failed");
    return (1);
  }

  matrix_product(npar, npar, 1, invhess.data(), TEMPVECTOR.data(), hgnc.data());

  signe = (flaginvsign) ? -1 :
                          1.;
  for (i = 0; i < npar; i++)
  {
    value = gauss[i + i * npar];
    value = (value < EPSILON10) ? 1 :
                                  value;
    hgnc[i] = signe * hgnc[i] / sqrt(value);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Add  ncont linear constraints (AX=B) to the linear system
 **  (of size npar) used for quadratic optimization
 **
 ** \param[in]  acont      Matrix A
 **
 ** \param[out] grad       left hand-side of the system
 ** \param[out] gauss      matrix of the system
 **
 *****************************************************************************/
static void st_fill_constraints(const VectorDouble &acont,
                                VectorDouble &grad,
                                VectorDouble &gauss)
{
  int icont, ipar;

  if (NCONT <= 0) return;
  for (icont = 0; icont < NCONT; icont++)
  {
    grad[NPAR + icont] = 0;
    for (ipar = 0; ipar < NPAR; ipar++)
    {
      GAUSS(ipar,NPAR+icont) = ACONT(ipar, icont);
      GAUSS(NPAR+icont,ipar) = ACONT(ipar, icont);
    }
  }
}

/****************************************************************************/
/*!
 **  Evaluate Gradient and Hermitian matrices
 **
 ** \param[in]  param     Current values of the parameters
 ** \param[in]  lower     Array of lower values
 ** \param[in]  upper     Array of upper values
 ** \param[in]  scale     Array of scaling values
 ** \param[in]  acont     Array of constraints
 ** \param[in]  tabwgt    Array of weights at control points
 **
 ** \param[out] residuals  Array of residuals
 ** \param[out] Jr         Array of gradients
 ** \param[out] grad       Gradient matrix
 ** \param[out] gauss      Gauss matrix
 ** \param[out] invhess    Inverse Hessian matrix
 ** \param[out] hgnc       Resulting hgnc array
 ** \param[out] param1     Working array (Dimension: NPAR)
 ** \param[out] param2     Working array (Dimension: NPAR)
 ** \param[out] tabmod1    Working array (Dimension: NDAT)
 ** \param[out] tabmod2    Working array (Dimension: NDAT)
 **
 *****************************************************************************/
static int st_calcul0(VectorDouble &param,
                      VectorDouble &lower,
                      VectorDouble &upper,
                      VectorDouble &scale,
                      const VectorDouble &acont,
                      VectorDouble &tabwgt,
                      VectorDouble &residuals,
                      VectorDouble &Jr,
                      VectorDouble &grad,
                      VectorDouble &gauss,
                      VectorDouble &invhess,
                      VectorDouble &hgnc,
                      VectorDouble &param1,
                      VectorDouble &param2,
                      VectorDouble &tabmod1,
                      VectorDouble &tabmod2)
{
  int error;

  /* Initializations */

  error = 0;

  st_gradient(param, lower, upper, scale, tabwgt, Jr, param1, param2, tabmod1,
              tabmod2);
  matrix_product(1, NDAT, NPAR, residuals.data(), Jr.data(), grad.data());
  st_determine_gauss(Jr, gauss);
  st_fill_constraints(acont, grad, gauss);
  error = st_solve_hgnc(NPAR + NCONT, grad, gauss, invhess, hgnc, 1);
  return (error);
}

/****************************************************************************/
/*!
 **  Calculate the number of new inactive constraints
 **
 ** \return  Number of highlited constraints
 **
 ** \param[in]  npar       Current number of parameters
 ** \param[in]  bords      Array containing the bounds
 ** \param[in]  ai         AI matrix
 ** \param[in]  hgnc       Resulting hgnc array
 **
 ** \param[out] flag       Working array
 ** \param[out] temp       Working array
 **
 *****************************************************************************/
static int st_possibilities(int npar,
                            VectorDouble &bords,
                            VectorDouble &ai,
                            VectorDouble &hgnc,
                            VectorInt &flag,
                            VectorDouble &temp)
{
  int ic, ipar, ipar2, flag_imposs, n_imposs;

  matrix_product(2 * npar, npar, 1, ai.data(), hgnc.data(), temp.data());

  n_imposs = 0;
  for (ic = ipar2 = 0; ic < 2; ic++)
    for (ipar = 0; ipar < npar; ipar++, ipar2++)
    {
      flag_imposs = ((ABS(bords[ipar2]) < TOLVAL)
          && (temp[ipar2] * SIGNE(ic) < 0));
      flag[ipar2] = (!flag_imposs);
      if (flag_imposs) n_imposs++;
    }
  return (n_imposs);
}

/****************************************************************************/
/*!
 **  Calculate the number of constraints
 **
 ** \param[in]  mode       Type of constraints
 ** \li                    -1 : the constraints which are non positive
 ** \li                     0 : the constraints must be zero
 ** \li                     1 : the constraints must be zero (cumul)
 ** \param[in]  bords_red  Reduced array containing the bounds
 ** \param[in]  ai_red     Reduced AI matrix
 ** \param[in]  hgnc       Resulting hgnc array
 **
 ** \param[out] consts     Array of constraints
 ** \param[out] flag       Array of indices with zero valid constraint
 ** \param[out] temp       Working array
 **
 *****************************************************************************/
static int st_define_constraints(int mode,
                                 VectorDouble &bords_red,
                                 VectorDouble &ai_red,
                                 VectorDouble &hgnc,
                                 VectorDouble &consts,
                                 VectorInt &flag,
                                 VectorDouble &temp)
{
  int ic, iparac, iparac2, number, flag_loc;

  /* Calculate the constraints */

  matrix_product(NPARAC2, NPARAC, 1, ai_red.data(), hgnc.data(), temp.data());

  for (ic = iparac2 = 0; ic < 2; ic++)
    for (iparac = 0; iparac < NPARAC; iparac++, iparac2++)
    {
      consts[iparac2] = (temp[iparac2] - bords_red[iparac2]) * SIGNE(ic);
      if (ABS(consts[iparac2]) < TOLVAL) consts[iparac2] = 0.;
    }

  /* Count the number of constraints */

  number = flag_loc = 0;
  for (ic = iparac2 = 0; ic < 2; ic++)
    for (iparac = 0; iparac < NPARAC; iparac++, iparac2++)
    {
      switch (mode)
      {
        case 0:
          flag_loc = (ABS(consts[iparac2]) < TOLVAL);
          break;

        case 1:
          flag_loc = flag[iparac2];
          if (!flag_loc) flag_loc = (ABS(consts[iparac2]) < TOLVAL);
          break;

        case -1:
          flag_loc = (consts[iparac2] < 0);
          break;
      }

      /* If a constraint is active for both lower and upper cases */
      /* It means that the lower and upper bounds are equal */
      /* Retain only one of the two constraints */
      if (flag_loc && ic == 1 && flag[iparac2 - NPARAC]) flag_loc = 0;
      flag[iparac2] = flag_loc;
      if (flag_loc) number++;
    }
  if (number > NPARAC) messageAbort("Fatal error in st_define_constraints");

  return (number);
}

/****************************************************************************/
/*!
**  Calculate the minimum of the criterion and update hgnadm
**
** \param[in]  flag         Array of indices with negative valid constraint
** \param[in]  bords_red    Reduced array containing the bounds
** \param[in]  top          Calculated top of the fraction
** \param[in]  bot          Calculated bottom of the fraction
** \param[in]  hgnc         Hgnc array
**
** \param[out] hgnadm       Admissible Hgn array
**
*****************************************************************************/
static void st_minimum(VectorInt& /*ind_util*/,
                       VectorInt& flag,
                       VectorDouble& bords_red,
                       const VectorDouble& top,
                       const VectorDouble& bot,
                       VectorDouble& hgnc,
                       VectorDouble& hgnadm)
{
  int ic, iparac, iparac2, jparac;
  double alpha, alpha_inf, bordval;

  jparac = -1;
  bordval = -1.e30;
  alpha_inf = 1.e30;

  for (ic = iparac2 = 0; ic < 2; ic++)
    for (iparac = 0; iparac < NPARAC; iparac++, iparac2++)
    {
      if (!flag[iparac2]) continue;
      alpha = BORDS_RED(ic, iparac);
      if (!top.empty()) alpha -= top[iparac2];

      /* Before dividing by bot, check that it not zero */
      if (bot[iparac2] != 0.) alpha /= bot[iparac2];
      alpha = ABS(alpha);
      if (alpha < alpha_inf)
      {
        alpha_inf = alpha;
        jparac = iparac;
        bordval = BORDS_RED(ic, iparac);
      }
    }
  if (jparac < 0) messageAbort("Fatal error in st_minimum");

  for (iparac = 0; iparac < NPARAC; iparac++)
    hgnadm[iparac] += alpha_inf * (hgnc[iparac] - hgnadm[iparac]);
  hgnadm[jparac] = bordval;

  return;
}

/****************************************************************************/
/*!
 **  Eliminate the useless constraints
 **
 ** \return  1 if the number of possibilities is reduced down to zero
 ** \return  or when the Generalized Inverse calculation failed
 **
 ** \param[in]  bords      Array containing the bounds
 ** \param[in]  ai         AI matrix
 ** \param[in]  grad       Gradient matrix
 ** \param[in]  gauss      Gaussian matrix
 ** \param[in]  invhess    Inverse Hessian matrix
 ** \param[in]  hgnc       hgnc array
 **
 ** \param[out]  ind_util   List of retained constraint indices
 ** \param[out]  bords_red  Reduced Bounds array
 ** \param[out]  ai_red     Reduced AI matrix
 ** \param[out]  grad_red   Reduced Gradient matrix
 ** \param[out]  gauss_red  Reduced Gauss matrix
 ** \param[out]  flag1      Working array
 ** \param[out]  flag2      Working array
 ** \param[out]  temp       Working array
 **
 *****************************************************************************/
static int st_suppress_unused_constraints(VectorDouble &bords,
                                          VectorDouble &ai,
                                          VectorDouble &grad,
                                          VectorDouble &gauss,
                                          VectorDouble &invhess,
                                          VectorDouble &hgnc,
                                          VectorInt &ind_util,
                                          VectorDouble &bords_red,
                                          VectorDouble &ai_red,
                                          VectorDouble &grad_red,
                                          VectorDouble &gauss_red,
                                          VectorInt &flag1,
                                          VectorInt &flag2,
                                          VectorDouble &temp)
{
  int n_imposs, ic, ipar, jpar, iparac, jparac, ipar2, iparac2;

  /* Get the set of constraints to be discarded */

  for (ipar2 = 0; ipar2 < NPAR2; ipar2++)
    flag1[ipar2] = 1;
  n_imposs = st_possibilities(NPAR, bords, ai, hgnc, flag2, temp);
  if (n_imposs >= NPAR) return (1);

  do
  {
    /* Patch the memory flag */
    for (ipar2 = iparac2 = 0; ipar2 < NPAR2; ipar2++)
    {
      if (flag1[ipar2] == 0) continue;
      flag1[ipar2] = flag2[iparac2];
      iparac2++;
    }

    /* Calculate the new ind array */
    for (ic = ipar2 = 0; ic < 2; ic++)
      for (ipar = 0; ipar < NPAR; ipar++, ipar2++)
        if (!flag1[ipar2]) ind_util[ipar] = 0;

    /* Update the value of variable NPARAC */
    NPARAC = 0;
    for (ipar = 0; ipar < NPAR; ipar++)
      if (POSSIBLE(ipar)) NPARAC++;
    NPARAC2 = 2 * NPARAC;
    /* This test has been added in order to avoid continuing */
    /* when no constraint remains */
    if (NPARAC <= 0) return (1);

    /* Reduce the arrays grad and gauss */
    for (ipar = iparac = 0; ipar < NPCT; ipar++)
    {
      if (!POSSIBLE(ipar)) continue;
      grad_red[iparac] = grad[ipar];
      for (jpar = jparac = 0; jpar < NPCT; jpar++)
      {
        if (!POSSIBLE(jpar)) continue;
        GAUSS_RED(iparac,jparac) = GAUSS(ipar, jpar);
        jparac++;
      }
      iparac++;
    }

    /* Reduce the arrays bords and ai */
    for (ic = 0; ic < 2; ic++)
      for (ipar = iparac = 0; ipar < NPAR; ipar++)
      {
        if (!POSSIBLE(ipar)) continue;
        BORDS_RED(ic,iparac) = BORDS(ic, ipar);
        for (jpar = jparac = 0; jpar < NPAR; jpar++)
        {
          if (!POSSIBLE(jpar)) continue;
          AI_RED(ic,iparac,jparac) = AI(ic, ipar, jpar);
          jparac++;
        }
        iparac++;
      }

    if (n_imposs > 0)
    {

      /* Update the Hessian and gradient matrices */

      if (st_solve_hgnc(NPARAC + NCONT, grad_red, gauss_red, invhess, hgnc, 1))
        return (1);

      /* Update the number of constraints */

      n_imposs = st_possibilities(NPARAC, bords_red, ai_red, hgnc, flag2, temp);
      if (n_imposs >= NPARAC) return (1);
    }
  }
  while (n_imposs);

  return (0);
}

/****************************************************************************/
/*!
 **  Update the constraints when no move has been performed
 **
 ** \param[in]  bords      Array containing the bounds
 ** \param[in]  ind_util   List of retained constraint indices
 **
 ** \param[out]  bords_red Reduced Bounds array
 **
 *****************************************************************************/
static void st_update_constraints(VectorDouble &bords,
                                  VectorInt &ind_util,
                                  VectorDouble &bords_red)
{
  int ipar, iparac, ic;

  /* Reduce the arrays bords */

  for (ic = 0; ic < 2; ic++)
    for (ipar = iparac = 0; ipar < NPAR; ipar++)
    {
      if (!POSSIBLE(ipar)) continue;
      BORDS_RED(ic,iparac) = BORDS(ic, ipar);
      iparac++;
    }
}

/****************************************************************************/
/*!
**  Minimization under constraints
**
** \return  Error returned code
**
** \param[in]  nactive      Number of active constraints
** \param[in]  flag_active  Array of indices with zero valid constraint
** \param[in]  bords_red    Reduced array containing the bounds
** \param[in]  ai_red       Reduced AI matrix
** \param[in]  grad_red     Reduced Gradient matrix
** \param[in]  gauss_red    Reduced Gauss matrix
**
** \param[out] lambda_neg   Index of the first negative lambda value
** \param[out] hgnc         Resulting Hgnc array
** \param[out] a            Minimization L.H.S. matrix
** \param[out] b            Minimization R.H.S. matrix
** \param[out] temp         Working array
**
*****************************************************************************/
static int st_establish_minimization(int     nactive,
                                     VectorInt& /*ind_util*/,
                                     VectorInt& flag_active,
                                     VectorDouble& bords_red,
                                     VectorDouble& ai_red,
                                     VectorDouble& grad_red,
                                     VectorDouble& gauss_red,
                                     int* lambda_neg,
                                     VectorDouble& hgnc,
                                     VectorDouble& a,
                                     VectorDouble& b,
                                     VectorDouble& temp)
{
  int size, i, ic, iparac, jparac, iparac2, iecr;

  /* Initialization */

  *lambda_neg = -1;
  size = NPARAC + NCONT + nactive;
  for (i = 0; i < size * size; i++)
    a[i] = 0.;
  for (i = 0; i < size; i++)
    b[i] = 0.;
  VectorDouble TEMPAUX(NPCT2 * NPCT2);

  /* Fill the L.H.S. and R.H.S. matrices */

  for (iparac = 0; iparac < NPARAC + NCONT; iparac++)
  {
    b[iparac] = -grad_red[iparac];
    for (jparac = 0; jparac < NPARAC + NCONT; jparac++)
      A(iparac,jparac) = GAUSS_RED(iparac, jparac);
  }

  for (ic = iecr = 0; ic < 2; ic++)
    for (iparac = 0; iparac < NPARAC; iparac++)
    {
      if (!FLAG_ACTIVE(ic, iparac)) continue;
      b[NPARAC + NCONT + iecr] = BORDS_RED(ic, iparac);
      for (jparac = 0; jparac < NPARAC; jparac++)
      {
        A(NPARAC + NCONT + iecr,jparac) = AI_RED(ic, iparac, jparac);
        A(jparac,NPARAC + NCONT + iecr) = AI_RED(ic, iparac, jparac);
      }
      iecr++;
    }

  /* Solve the system */

  if (st_solve_hgnc(size, b, a, TEMPAUX, temp, 0)) return (1);

  /* Store the final result */

  for (iparac = 0; iparac < NPARAC; iparac++)
    hgnc[iparac] = temp[iparac];
  for (ic = 0; ic < 2; ic++)
    for (iparac = 0; iparac < NPARAC; iparac++)
      if (FLAG_ACTIVE(ic, iparac)) hgnc[iparac] = BORDS_RED(ic, iparac);

  /* Check for the first obviously negative coefficient */

  for (ic = iecr = iparac2 = 0; ic < 2; ic++)
    for (iparac = 0; iparac < NPARAC; iparac++, iparac2++)
    {
      if (!FLAG_ACTIVE(ic, iparac)) continue;
      if (temp[NPARAC + iecr] < -TOLVAL)
      {
        *lambda_neg = iparac2;
        return (0);
      }
      iecr++;
    }
  return (0);
}

/****************************************************************************/
/*!
 **  Check that the algorithm is valid
 **
 ** \param[in]  ind_util   List of retained constraint indices
 ** \param[in]  hgnc       Working vector
 ** \param[in]  acont      Matrix of additional constraints
 **
 *****************************************************************************/
static void st_check(VectorInt &ind_util,
                     VectorDouble &hgnc,
                     const VectorDouble &acont)
{
  double temp;
  int ipar, icont, iparac;

  if (NCONT <= 0) return;

  for (icont = 0; icont < NCONT; icont++)
  {
    temp = 0;
    for (iparac = ipar = 0; ipar < NPAR; ipar++)
      if (POSSIBLE(ipar)) temp += ACONT(ipar,icont) * hgnc[iparac++];

    if (ABS(temp) > 1e-3)
      messageAbort(
          "The constraints are not fulfilled. This should never happen");
  }
  return;
}

/****************************************************************************/
/*!
 **  Minimization under constraints
 **
 ** \return  1 if the optimum has not been reached (various reasons)
 **
 ** \param[in]  ind_util   List of retained constraint indices
 ** \param[in]  bords_red  Reduced array containing the bounds
 ** \param[in]  ai_red     Reduced AI matrix
 ** \param[in]  grad_red   Reduced Gradient matrix
 ** \param[in]  gauss_red  Reduced Gauss matrix
 **
 ** \param[out]  consts       Array of constraints
 ** \param[out]  hgnc         Working vector
 ** \param[out]  hgnadm       Working array
 ** \param[out]  flag_active  Array of indices with zero valid constraint
 ** \param[out]  flag_actaux  Array of indices with negative valid constraint
 ** \param[out]  a            Minimization L.H.S. matrix
 ** \param[out]  b1           Minimization R.H.S. matrix
 ** \param[out]  b2           Minimization R.H.S. matrix
 ** \param[out]  b3           Minimization R.H.S. matrix
 ** \param[out]  temp         Working array
 ** \param[out]  acont        Constraint array
 **
 *****************************************************************************/
static int st_minimization_under_constraints(VectorInt &ind_util,
                                             VectorDouble &bords_red,
                                             VectorDouble &ai_red,
                                             VectorDouble &grad_red,
                                             VectorDouble &gauss_red,
                                             VectorDouble &consts,
                                             VectorDouble &hgnc,
                                             VectorDouble &hgnadm,
                                             VectorInt &flag_active,
                                             VectorInt &flag_actaux,
                                             VectorDouble &a,
                                             VectorDouble &b1,
                                             VectorDouble &b2,
                                             VectorDouble &b3,
                                             VectorDouble &temp,
                                             const VectorDouble &acont)
{
  int iparac, nactaux, sortie, nactive, lambda_neg;
  double min_adm_cur, min_adm_best;
  static int nitermax = 2000;

  /* Calculate the constraints vector */

  nactaux = st_define_constraints(-1, bords_red, ai_red, hgnc, consts,
                                  flag_actaux, temp);
  if (nactaux <= 0)
  {
    for (iparac = 0; iparac < NPARAC; iparac++)
      hgnadm[iparac] = hgnc[iparac];
    st_check(ind_util, hgnadm, acont);
    return (0);
  }

  /* Find an initial admissible point */

  for (iparac = 0; iparac < NPARAC; iparac++)
    hgnadm[iparac] = 0.;
  matrix_product(NPARAC2, NPARAC, 1, ai_red.data(), hgnc.data(), b1.data());
  st_minimum(ind_util, flag_actaux, bords_red, VectorDouble(), b1, hgnc,
             hgnadm);
  st_check(ind_util, hgnadm, acont);

  /* Calculate the constraints vector */

  nactive = st_define_constraints(0, bords_red, ai_red, hgnadm, consts,
                                  flag_active, temp);
  min_adm_best = st_essai(hgnadm, grad_red, gauss_red);
  if (VERBOSE_GQO && debug_query("converge"))
    message("GQO(  0) : Gain for initial solution  = %lg\n", -min_adm_best);

  sortie = SOUSITER = 0;
  while (!sortie && SOUSITER + 1 <= nitermax)
  {
    SOUSITER++;

    /* Load the minimization matrix */

    if (st_establish_minimization(nactive, ind_util, flag_active, bords_red,
                                  ai_red, grad_red, gauss_red, &lambda_neg,
                                  hgnc, a, b1, temp)) return (1);
    st_check(ind_util, hgnc, acont);

    /* Calculate the constraints vector */

    nactaux = st_define_constraints(-1, bords_red, ai_red, hgnc, consts,
                                    flag_actaux, temp);
    if (nactaux > 0)
    {
      for (iparac = 0; iparac < NPARAC; iparac++)
        b3[iparac] = hgnc[iparac] - hgnadm[iparac];
      matrix_product(NPARAC2, NPARAC, 1, ai_red.data(), hgnadm.data(),
                     b1.data());
      matrix_product(NPARAC2, NPARAC, 1, ai_red.data(), b3.data(), b2.data());
      st_minimum(ind_util, flag_actaux, bords_red, b1, b2, hgnc, hgnadm);
      st_check(ind_util, hgnadm, acont);

      nactive = st_define_constraints(1, bords_red, ai_red, hgnadm, consts,
                                      flag_active, temp);
      min_adm_cur = st_essai(hgnadm, grad_red, gauss_red);
      if (VERBOSE_GQO && debug_query("converge"))
        message("GQO(%3d) : Gain for infeasible case   = %lg\n", SOUSITER,
                -min_adm_cur);
      if (min_adm_cur >= min_adm_best) break;
      min_adm_best = min_adm_cur;
    }
    else
    {
      for (iparac = 0; iparac < NPARAC; iparac++)
        hgnadm[iparac] = hgnc[iparac];
      if (lambda_neg >= 0)
      {
        flag_active[lambda_neg] = 0;
        nactive--;
        if (VERBOSE_GQO && debug_query("converge"))
          message("GQO(%3d) : Gain for feasible case     = %lg\n", SOUSITER,
                  -st_essai(hgnadm, grad_red, gauss_red));
      }
      else
      {
        sortie = 1;
      }
    }
  }
  st_check(ind_util, hgnadm, acont);

  return (SOUSITER >= nitermax);
}

/****************************************************************************/
/*!
 **  Initialize the constraints
 **
 ** \param[in]  ind_util     List of retained constraint indices
 ** \param[in]  ai           AI matrix
 **
 *****************************************************************************/
static void st_constraints_init(VectorInt &ind_util, VectorDouble &ai)
{
  int ipar, jpar, ic;

  for (ipar = 0; ipar < NPCT; ipar++)
    ind_util[ipar] = 1;

  for (ic = 0; ic < 2; ic++)
    for (ipar = 0; ipar < NPAR; ipar++)
      for (jpar = 0; jpar < NPAR; jpar++)
        AI(ic,ipar,jpar) = (ipar == jpar);
}

/****************************************************************************/
/*!
 **  Evaluate the bounds and check if one constraint is on its boudn
 **
 ** \param[in]  param      Current values of the parameters
 ** \param[in]  lower      Array of lower values
 ** \param[in]  upper      Array of upper values
 ** \param[in]  scale      Array of scaling values
 ** \param[in]  delta      Radius of the trusting area
 **
 ** \param[out] bords      Value for the bounds
 **
 *****************************************************************************/
static void st_define_bounds(VectorDouble &param,
                             VectorDouble &lower,
                             VectorDouble &upper,
                             VectorDouble &scale,
                             double delta,
                             VectorDouble &bords)
{
  double dloc, diff;
  int ipar;

  for (ipar = 0; ipar < NPAR; ipar++)
  {
    dloc = delta * scale[ipar];

    /* Lower bound */

    if (FFFF(lower[ipar]))
      BORDS(0,ipar) = -dloc;
    else
    {
      diff = param[ipar] - lower[ipar];
      if (diff > dloc) diff = dloc;
      if (ABS(param[ipar] - diff) < dloc / 10.) diff /= 2.;
      BORDS(0,ipar) = -diff;
    }

    /* Upper bound */

    if (FFFF(upper[ipar]))
      BORDS(1,ipar) = dloc;
    else
    {
      diff = upper[ipar] - param[ipar];
      if (diff > dloc) diff = dloc;
      if (ABS(param[ipar] - diff) < dloc / 10.) diff /= 2.;
      BORDS(1,ipar) = diff;
    }
  }
}

/****************************************************************************/
/*!
 **  Display the title for FOXLEG trace
 **
 *****************************************************************************/
static void st_foxleg_debug_title(void)

{
  int ipar;
  static char string[10];

  if (!debug_query("converge")) return;
  mestitle(1, "Trajectory of parameters in Foxleg Algorithm");
  tab_prints(NULL, 1, EJustify::RIGHT, "Iteration");
  tab_prints(NULL, 1, EJustify::RIGHT, "Score");
  tab_prints(NULL, 1, EJustify::RIGHT, "Delta");
  for (ipar = 0; ipar < NPAR; ipar++)
  {
    (void) gslSPrintf(string, "Par-%d", ipar + 1);
    tab_prints(NULL, 1, EJustify::RIGHT, string);
  }
  message("\n");
}

/****************************************************************************/
/*!
 **  Display the current status for FOXLEG trace
 **
 *****************************************************************************/
static void st_foxleg_debug_current(double mscur,
                                    double delta,
                                    VectorDouble &param)
{
  int ipar;

  if (!debug_query("converge")) return;
  tab_printi(NULL, 1, EJustify::RIGHT, ITERATION);
  tab_printd(NULL, 1, EJustify::RIGHT, mscur);
  tab_printd(NULL, 1, EJustify::RIGHT, delta);
  for (ipar = 0; ipar < NPAR; ipar++)
    tab_printg(NULL, 1, EJustify::RIGHT, param[ipar]);
  message("\n");
}

/****************************************************************************/
/*!
 **  Display the FOXLEG score
 **
 *****************************************************************************/
static void st_foxleg_score(const Option_AutoFit &mauto,
                            double mscur,
                            double delta,
                            double arret)
{
  if (mauto.getVerbose() > 0)
  {
    mestitle(1, "Statistics for the Minimization Foxleg procedure");
    message("- Number of experimental values = %d\n", NDAT);
    message("- Number of parameters          = %d\n", NPAR);
    message("- Number of iterations          = %d/%d \n", ITERATION,
            mauto.getMaxiter());
    message("- Value of minimized function   = %lg\n", mscur);
    message("- Initial increment value       = %lg\n", mauto.getInitdelta());
    message("- Current increment value       = %lg\n", delta);
    message("- Increment Stopping Criterion  = %lg\n", mauto.getEpsdelta());
    message("- Stopping Value                = %lg\n", arret);
    message("- Stopping Criterion            = %lg\n", mauto.getTolstop());
    message("- Stopping Criterion (scaled)   = %lg\n", mauto.getTolred());
  }
}

/****************************************************************************/
/*!
 **  Interpolate linearly the vector of parameters
 **
 ** \param[in]  mscur     Current minimization value
 ** \param[in]  param     Current values of the parameters
 ** \param[in]  acont     Array of constraints
 ** \param[in]  tabexp    Array of values at control points
 ** \param[in]  tabwgt    Array of weights at control points
 ** \param[in]  bords     Array containing the bounds
 ** \param[in]  grad      Gradient matrix
 **
 ** \param[out] msaux      New minimization value
 ** \param[out] paramaux   New vector of parameters
 ** \param[out] residuals  Array of residuals
 ** \param[out] tabmod1    Working array (Dimension: NDAT)
 **
 *****************************************************************************/
static void st_linear_interpolate(double mscur,
                                  VectorDouble &param,
                                  const VectorDouble &acont,
                                  VectorDouble &tabexp,
                                  VectorDouble &tabwgt,
                                  VectorDouble &bords,
                                  VectorDouble &grad,
                                  double *msaux,
                                  VectorDouble &paramaux,
                                  VectorDouble &residuals,
                                  VectorDouble &tabmod1)
{
  double alpha, shift;
  int ipar, icont, flag_ok;

  alpha = 100.;
  while (1)
  {
    alpha /= 2.;
    flag_ok = 1;

    /* Modify the vector of parameters (avoiding edges) */

    for (ipar = 0; ipar < NPAR && flag_ok; ipar++)
    {
      shift = -alpha * grad[ipar];
      if ((shift < 0 && ABS(BORDS(0,ipar)) < TOLVAL) || (shift > 0
          && ABS(BORDS(1,ipar)) < TOLVAL)) shift = 0.;

      /* Discard any move if constrained minimization */
      for (icont = 0; icont < NCONT; icont++)
        if (ACONT(ipar,icont) != 0.) shift = 0.;

      paramaux[ipar] = param[ipar] + shift;
      if (shift < BORDS(0,ipar) - TOLVAL || shift > BORDS(1,ipar) + TOLVAL)
        flag_ok = 0;
    }
    if (!flag_ok) continue;

    /* Check that we are in the descent orientation */

    *msaux = st_residuals(paramaux, tabexp, tabwgt, tabmod1, residuals);
    if (*msaux <= mscur) return;
  }
}

/****************************************************************************/
/*!
 **  Check if the default, lower and upper bounds of the parameters
 **  are consistent or not
 **
 ** \return  Error returned code
 **
 ** \param[in]  param         Current values of the parameters
 ** \param[in]  lower         Array of lower values
 ** \param[in]  upper         Array of upper values
 **
 *****************************************************************************/
static int st_check_param(VectorDouble &param,
                          VectorDouble &lower,
                          VectorDouble &upper)
{
  int ipar;

  /* Check lower vs upper bounds */

  if (!lower.empty() && !upper.empty()) for (ipar = 0; ipar < NPAR; ipar++)
  {
    if (FFFF(lower[ipar]) || FFFF(upper[ipar])) continue;
    if (lower[ipar] <= upper[ipar]) continue;
    messerr("Error in parameter %d: lower bound (%lf) > upper bound (%lf)", ipar+1,lower[ipar],upper[ipar]);
    return (1);
  }

  /* Check lower bounds vs default values */

  if (!lower.empty() && !param.empty()) for (ipar = 0; ipar < NPAR; ipar++)
  {
    if (FFFF(lower[ipar]) || FFFF(param[ipar])) continue;
    if (param[ipar] < lower[ipar]) param[ipar] = lower[ipar];
  }

  /* Check upper bounds vs default values */

  if (!upper.empty() && !param.empty()) for (ipar = 0; ipar < NPAR; ipar++)
  {
    if (FFFF(upper[ipar]) || FFFF(param[ipar])) continue;
    if (param[ipar] > upper[ipar]) param[ipar] = upper[ipar];
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Foxleg algorithm
 **
 ** \return  Error returned code
 ** \return  0 : Correct termination
 ** \return -1 : The convergence has not been reached
 ** \return  1 : Core problem
 **
 ** \param[in]  ndat          Number of control points
 ** \param[in]  npar          Number of parameters to estimate
 ** \param[in]  ncont         Number of additional constraints
 ** \param[in]  acont         Matrix of additional constraints
 **                           (Dimension = ncont * npar)
 ** \param[in]  param         Current values of the parameters
 ** \param[in]  lower         Array of lower values
 ** \param[in]  upper         Array of upper values
 ** \param[in]  scale         Array of scale
 ** \param[in]  mauto         Option_AutoFit structure
 ** \param[in]  flag_title    Print the title after func_evaluate()
 ** \param[in]  func_evaluate Function for evaluating the model
 ** \param[in]  tabexp        Array of values at control points
 ** \param[in]  tabwgt        Array of weights at control points
 **
 ** \remark  The evaluation function func_evaluate() is called with the
 ** \remark  following arguments:
 ** \remark  ndat         : Number of control points
 ** \remark  npar         : Number of parameters
 ** \remark  param        : Vector of current values for the parameters
 ** \remark  work         : Output vector for the values at control points
 ** \remark                 evaluated using the current parameter values
 **
 ** \remark  Some additional constraints may be applied on the parameters
 ** \remark  When not used, we must set: ncont=0, acont=NULL
 **
 *****************************************************************************/
int foxleg_f(int ndat,
                             int npar,
                             int ncont,
                             const VectorDouble &acont,
                             VectorDouble &param,
                             VectorDouble &lower,
                             VectorDouble &upper,
                             VectorDouble &scale,
                             const Option_AutoFit &mauto,
                             int flag_title,
                             void (*func_evaluate)(int ndat,
                                                   int npar,
                                                   VectorDouble &param,
                                                   VectorDouble &work),
                             VectorDouble &tabexp,
                             VectorDouble &tabwgt)
{
  int error, ipar, iparac, flag_cont, ipct, jpct, flag_moved;
  double mscur, ms0, msaux, rho, arret, delta, denom;

  /* Preliminary checks */

  error = 1;
  arret = 0.;
  delta = mauto.getInitdelta();
  NDAT = ndat;
  NPAR = NPARAC = npar;
  NCONT = ncont;
  NPAR2 = NPARAC2 = npar * 2;
  NPCT = NPAR + NCONT;
  NPCT2 = NPAR2 + NCONT;
  FUNC_EVALUATE = func_evaluate;

  /* Core allocation */

  VectorInt ind_util(NPCT);
  VectorInt flag_active(NPAR2);
  VectorInt flag_actaux(NPAR2);
  VectorDouble a(NPCT2 * NPCT2);
  VectorDouble b1(NPCT2);
  VectorDouble b2(NPAR2);
  VectorDouble b3(NPAR2);
  VectorDouble temp(NPCT2);
  VectorDouble param1(NPAR);
  VectorDouble param2(NPAR);
  VectorDouble grad(NPCT);
  VectorDouble grad_red(NPCT);
  VectorDouble hgn(NPAR);
  VectorDouble hgnc(NPCT);
  VectorDouble hgnadm(NPCT);
  VectorDouble paramaux(NPAR);
  VectorDouble gauss(NPCT * NPCT);
  VectorDouble gauss_red(NPCT * NPCT);
  VectorDouble invhess(NPCT * NPCT);
  VectorDouble Jr(NPAR * NDAT);
  VectorDouble residuals(NDAT);
  VectorDouble tabmod1(NDAT);
  VectorDouble tabmod2(NDAT);
  VectorDouble ai(NPAR * NPAR2);
  VectorDouble ai_red(NPAR * NPAR2);
  VectorDouble consts(NPAR2);
  VectorDouble bords(NPAR2);
  VectorDouble bords_red(NPAR2);

  /* Preliminary check */

  if (st_check_param(param, lower, upper)) goto label_end;

  /* Initializations */

  for (ipar = 0; ipar < NPAR; ipar++)
  {
    paramaux[ipar] = param[ipar];
    hgn[ipar] = 0.;
  }
  for (ipar = 0; ipar < NPCT; ipar++)
    hgnc[ipar] = hgnadm[ipar] = 0.;
  for (ipct = 0; ipct < NPCT; ipct++)
    for (jpct = 0; jpct < NPCT; jpct++)
      GAUSS(ipct,jpct) = 0.;

  st_constraints_init(ind_util, ai);

  /* Calculate the gradient */
  ms0 = mscur = st_residuals(param, tabexp, tabwgt, tabmod1, residuals);
  if (st_calcul0(param, lower, upper, scale, acont, tabwgt, residuals, Jr, grad,
                 gauss, invhess, hgnc, param1, param2, tabmod1, tabmod2))
    goto label_end;

  st_foxleg_debug_title();

  /***********************/
  /* Iterative procedure */
  /***********************/

  ITERATION = 0;
  flag_moved = 1;
  flag_cont = (mauto.getMaxiter() <= 0) ? 0 :
                                          1;
  while (flag_cont)
  {
    ITERATION++;
    if (ITERATION > 1 && flag_title) st_foxleg_debug_title();
    st_foxleg_debug_current(mscur, delta, param);

    /* Calculate the bounds around the current parameter vector */

    st_define_bounds(param, lower, upper, scale, delta, bords);

    /* Flag out the impossible components of the parameter vector */

    if (flag_moved)
    {
      if (st_suppress_unused_constraints(bords, ai, grad, gauss, invhess, hgnc,
                                         ind_util, bords_red, ai_red, grad_red,
                                         gauss_red, flag_active, flag_actaux,
                                         temp)) break;
    }
    else
    {
      st_update_constraints(bords, ind_util, bords_red);
    }

    /* Perform the minimization under constraints */

    if (st_minimization_under_constraints(ind_util, bords_red, ai_red, grad_red,
                                          gauss_red, consts, hgnc, hgnadm,
                                          flag_active, flag_actaux, a, b1, b2,
                                          b3, temp, acont))
    {
      if (debug_query("converge"))
      {
        messerr("Convergence not reached in minimization under constraints");
        messerr("The process is resumed in the final minimization status");
      }
    }

    /* Check for an early stop */

    if (ITERATION >= mauto.getMaxiter() || delta < mauto.getEpsdelta())
      flag_cont = 0;
    arret = 0.;
    for (iparac = 0; iparac < NPARAC; iparac++)
      arret += ABS(hgnc[iparac]);
    arret /= (ms0 - mscur + mauto.getEpsdelta());
    if (arret < mauto.getTolstop()) flag_cont = 0;

    /* Update values for the next iteration */

    for (ipar = iparac = 0; ipar < NPAR; ipar++)
    {
      hgn[ipar] = 0.;
      paramaux[ipar] = param[ipar];
      if (!POSSIBLE(ipar)) continue;
      hgn[ipar] = hgnadm[iparac++];
      paramaux[ipar] += hgn[ipar];
    }

    denom = st_essai(hgnadm, grad_red, gauss_red);
    if (ABS(denom) < EPSILON10) goto label_ok;
    if (denom > 0)
      st_linear_interpolate(mscur, param, acont, tabexp, tabwgt, bords, grad,
                            &msaux, paramaux, residuals, tabmod1);
    else
      msaux = st_residuals(paramaux, tabexp, tabwgt, tabmod1, residuals);
    rho = (msaux - mscur) / denom;

    if (msaux < mscur)
    {
      SOUSITER++;
      mscur = msaux;
      for (ipar = 0; ipar < NPAR; ipar++)
        param[ipar] = paramaux[ipar];
      if (st_calcul0(param, lower, upper, scale, acont, tabwgt, residuals, Jr,
                     grad, gauss, invhess, hgnc, param1, param2, tabmod1,
                     tabmod2)) goto label_end;
      st_constraints_init(ind_util, ai);
      flag_moved = 1;
      if (denom < 0 && rho > 0.75)
        delta = MAX(delta, 3. * st_norm_hgn(hgn, scale));
    }
    else
      flag_moved = 0;
    if (rho < 0.25 || denom > 0) delta /= 2.;
  }

  /* Result printout */

  label_ok: if (ITERATION >= mauto.getMaxiter())
  {
    error = -1;
    if (debug_query("converge")) messerr("Convergence has not been reached");
  }
  else
  {
    error = 0;
  }

  set_keypair("Foxleg Value", 1, 1, 1, &mscur);
  st_foxleg_score(mauto, mscur, delta, arret);

  label_end: return (error);
}

/****************************************************************************/
/*!
 **  Add constraints to the Option_AutoFit structure
 **
 ** \return Error return code
 **
 ** \param[in]  mauto        Option_AutoFit structure
 ** \param[in]  constantSill Constant value for the Sill as a constraint
 **
 *****************************************************************************/
int opt_mauto_add_constraints(Option_AutoFit &mauto,
                                              double constantSill)
{
  mauto.setConstantSillValue(constantSill);

  return (0);
}

/****************************************************************************/
/*!
 **  Add constraints (all equal to 1) to the Option_AutoFit structure
 **
 ** \return Error return code
 **
 ** \param[in]  mauto       Option_AutoFit structure
 **
 *****************************************************************************/
int opt_mauto_add_unit_constraints(Option_AutoFit &mauto)
{
  mauto.setConstantSillValue(1.);
  return (0);
}

