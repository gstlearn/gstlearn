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

#include "Enum/EJustify.hpp"

#include "Basic/Utilities.hpp"
#include "Basic/File.hpp"
#include "Basic/String.hpp"
#include "Basic/OptDbg.hpp"

#include <math.h>

/*! \cond */
#define IAD(ic,ipar)             ((ipar) + NPAR * (ic))
#define IADAC(ic,iparac)         ((iparac) + NPARAC * (ic))
#define AI(ic,ipar,jpar)         (ai[(jpar) * NPAR2 + IAD(ic,ipar)])
#define AI_RED(ic,iparac,jparac) (ai_red[(jparac) * NPARAC2 + IADAC(ic,iparac)])
#define POSSIBLE(ipar)           (ind_util[ipar] > 0)
#define SIGNE(ic)                ((ic == 0) ? +1 : -1)
#define FLAG_ACTIVE(ic,iparac)   (flag_active[IADAC(ic,iparac)])
/*! \endcond */

static int VERBOSE_GQO = 0;

static int NPAR, NPAR2, NPARAC, NPARAC2, NDAT, NCONT, NPCT, NPCT2;
static int ITERATION, SOUSITER;
static int USE_EIGEN_LIBRARY = 0;
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
                        MatrixRectangular& Jr,
                        VectorDouble &param1,
                        VectorDouble &param2,
                        VectorDouble &tabmod1,
                        VectorDouble &tabmod2)
{

  /* Calculate the gradients */

  double epsgrad = EPSILON3;
  for (int ipar = 0; ipar < NPAR; ipar++)
  {
    double epsloc = ABS(epsgrad * scale[ipar]);
    epsloc = MAX(epsgrad, epsloc);
    for (int jpar = 0; jpar < NPAR; jpar++)
      param1[jpar] = param2[jpar] = param[jpar];

    param1[ipar] = param[ipar] + epsloc;
    if (!FFFF(upper[ipar])) param1[ipar] = MIN(upper[ipar], param1[ipar]);
    param2[ipar] = param[ipar] - epsloc;
    if (!FFFF(lower[ipar])) param2[ipar] = MAX(lower[ipar], param2[ipar]);

    double ratio1 = epsloc;
    double ratio2 = epsloc;
    if (!FFFF(upper[ipar])) ratio1 = MIN(epsloc, upper[ipar] - param[ipar]);
    if (!FFFF(lower[ipar])) ratio2 = MIN(epsloc, param[ipar] - lower[ipar]);

    FUNC_EVALUATE(NDAT, NPAR, param1, tabmod1);
    FUNC_EVALUATE(NDAT, NPAR, param2, tabmod2);

    double bot = ratio1 + ratio2;
    for (int idat = 0; idat < NDAT; idat++)
    {
      double top = tabmod1[idat] - tabmod2[idat];
      double weight = (!tabwgt.empty()) ? tabwgt[idat] : 1.;
      Jr.setValue(idat,ipar, (bot != 0.) ? weight * top / bot : 0.);
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
  /* Evaluate the Model at conditioning points */

  FUNC_EVALUATE(NDAT, NPAR, param, tabmod);

  /* Evaluate the residuals */

  double msse = 0.;
  for (int idat = 0; idat < NDAT; idat++)
  {
    double weight = (!tabwgt.empty()) ? tabwgt[idat] : 1.;
    double value = weight * (tabmod[idat] - tabexp[idat]);
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
static void st_determine_gauss(MatrixRectangular& Jr, MatrixSquareGeneral& gauss)
{
  for (int ipar = 0; ipar < NPAR; ipar++)
    for (int jpar = 0; jpar < NPAR; jpar++)
    {
      double value = 0.;
      for (int idat = 0; idat < NDAT; idat++)
        value += Jr.getValue(idat,ipar) * Jr.getValue(idat, jpar);
      gauss.setValue(ipar,jpar,value);
    }
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
  double norme = 0.;
  for (int ipar = 0; ipar < NPAR; ipar++)
  {
    double v1 = ABS(hgn[ipar] / scale[ipar]);
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
                       MatrixSquareGeneral& gauss_red)
{
  double v1 = VH::innerProduct(hgnadm, grad_red);
  double v2 = gauss_red.normVec(hgnadm);
  double result = v1 + v2 / 2.;

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
 ** \param[in] hgnc         Resulting hgnc array
 ** \param[in] flaginvsign  if 1, the result is multiplied by -1
 **
 *****************************************************************************/
static int st_solve_hgnc(int npar,
                         const VectorDouble &grad,
                         const MatrixSquareGeneral& gauss,
                         VectorDouble &hgnc,
                         int flaginvsign)
{
  VectorDouble tempMatVD(npar * npar,0.);
  VectorDouble tempVec(npar,0.);
  MatrixSquareSymmetric tempMat(npar, USE_EIGEN_LIBRARY);
  double eps = EPSILON10;

  double signe = (flaginvsign) ? -1 : 1.;

  for (int i = 0; i < npar; i++)
  {
    double vali = gauss.getValue(i, i);
    vali = (isZero(vali,eps)) ? 1 : sqrt(vali);
    tempVec[i] = grad[i] / vali;
    for (int j = 0; j < npar; j++)
    {
      double valj = gauss.getValue(j, j);
      valj = (isZero(valj,eps)) ? 1 : sqrt(valj);
      tempMat.setValue(i, j, gauss.getValue(i, j) / (vali * valj));
    }
  }

  if (tempMat.computeGeneralizedInverse(tempMat))
  {
    messerr("Error: Singularity in the Generalized Inverse");
    messerr("The Automatic Fitting Procedure failed");
    return (1);
  }

  matrix_product_safe(npar, npar, 1, tempMat.getValues().data(), tempVec.data(),
                      hgnc.data());

  for (int i = 0; i < npar; i++)
  {
    double value = gauss.getValue(i,i);
    value = (isZero(value,eps)) ? 1 : sqrt(value);
    hgnc[i] = signe * hgnc[i] / value;
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
static void st_fill_constraints(const MatrixRectangular& acont,
                                VectorDouble &grad,
                                MatrixSquareGeneral& gauss)
{
  if (NCONT <= 0) return;
  for (int icont = 0; icont < NCONT; icont++)
  {
    grad[NPAR + icont] = 0;
    for (int ipar = 0; ipar < NPAR; ipar++)
    {
      gauss.setValue(ipar,NPAR+icont,acont.getValue(ipar, icont));
      gauss.setValue(NPAR+icont,ipar,acont.getValue(ipar, icont));
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
                      const MatrixRectangular& acont,
                      VectorDouble &tabwgt,
                      VectorDouble &residuals,
                      MatrixRectangular& Jr,
                      VectorDouble &grad,
                      MatrixSquareGeneral& gauss,
                      VectorDouble &hgnc,
                      VectorDouble &param1,
                      VectorDouble &param2,
                      VectorDouble &tabmod1,
                      VectorDouble &tabmod2)
{
  st_gradient(param, lower, upper, scale, tabwgt, Jr, param1, param2, tabmod1, tabmod2);
  matrix_product_safe(1, NDAT, NPAR, residuals.data(), Jr.getValues().data(), grad.data());
  st_determine_gauss(Jr, gauss);
  st_fill_constraints(acont, grad, gauss);
  return st_solve_hgnc(NPAR + NCONT, grad, gauss, hgnc, 1);
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
                            MatrixRectangular &bords,
                            VectorDouble &ai,
                            VectorDouble &hgnc,
                            VectorInt &flag,
                            VectorDouble &temp)
{
  int flag_imposs;

  matrix_product_safe(2 * npar, npar, 1, ai.data(), hgnc.data(), temp.data());

  int n_imposs = 0;
  int ipar2 = 0;
  for (int ic = 0; ic < 2; ic++)
    for (int ipar = 0; ipar < npar; ipar++, ipar2++)
    {
      flag_imposs = ((ABS(bords.getValue(ic,ipar)) < EPSILON9) && (temp[ipar2] * SIGNE(ic) < 0));
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
                                 MatrixRectangular &bords_red,
                                 VectorDouble &ai_red,
                                 VectorDouble &hgnc,
                                 MatrixRectangular &consts,
                                 VectorInt &flag,
                                 VectorDouble &temp)
{
  int iparac2;

  /* Calculate the constraints */

  matrix_product_safe(NPARAC2, NPARAC, 1, ai_red.data(), hgnc.data(), temp.data());

  iparac2 = 0;
  for (int ic = 0; ic < 2; ic++)
    for (int iparac = 0; iparac < NPARAC; iparac++, iparac2++)
    {
      consts.setValue(ic,iparac,(temp[iparac2] - bords_red.getValue(ic,iparac)) * SIGNE(ic));
      if (ABS(consts.getValue(ic,iparac)) < EPSILON9) consts.setValue(ic,iparac,0.);
    }

  /* Count the number of constraints */

  int number = 0;
  int flag_loc = 0;
  iparac2 = 0;
  for (int ic = 0; ic < 2; ic++)
    for (int iparac = 0; iparac < NPARAC; iparac++, iparac2++)
    {
      switch (mode)
      {
        case 0:
          flag_loc = (ABS(consts.getValue(ic,iparac)) < EPSILON9);
          break;

        case 1:
          flag_loc = flag[iparac2];
          if (!flag_loc) flag_loc = (ABS(consts.getValue(ic,iparac)) < EPSILON9);
          break;

        case -1:
          flag_loc = (consts.getValue(ic,iparac) < 0);
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
                       MatrixRectangular& bords_red,
                       const VectorDouble& top,
                       const VectorDouble& bot,
                       VectorDouble& hgnc,
                       VectorDouble& hgnadm)
{
  int jparac = -1;
  double bordval = -1.e30;
  double alpha_inf = 1.e30;

  int iparac2 = 0;
  for (int ic = 0; ic < 2; ic++)
    for (int iparac = 0; iparac < NPARAC; iparac++, iparac2++)
    {
      if (!flag[iparac2]) continue;
      double alpha = bords_red.getValue(ic, iparac);
      if (!top.empty()) alpha -= top[iparac2];

      /* Before dividing by bot, check that it not zero */
      if (bot[iparac2] != 0.) alpha /= bot[iparac2];
      alpha = ABS(alpha);
      if (alpha < alpha_inf)
      {
        alpha_inf = alpha;
        jparac = iparac;
        bordval = bords_red.getValue(ic, iparac);
      }
    }
  if (jparac < 0) messageAbort("Fatal error in st_minimum");

  for (int iparac = 0; iparac < NPARAC; iparac++)
    hgnadm[iparac] += alpha_inf * (hgnc[iparac] - hgnadm[iparac]);
  hgnadm[jparac] = bordval;

  return;
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
static void st_update_bords(MatrixRectangular &bords,
                            VectorInt &ind_util,
                            MatrixRectangular &bords_red)
{
  for (int ic = 0; ic < 2; ic++)
  {
    int iparac = 0;
    for (int ipar = 0; ipar < NPAR; ipar++)
    {
      if (!POSSIBLE(ipar)) continue;
      bords_red.setValue(ic,iparac,bords.getValue(ic, ipar));
      iparac++;
    }
  }
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
static int st_suppress_unused_constraints(MatrixRectangular &bords,
                                          VectorDouble &ai,
                                          VectorDouble &grad,
                                          MatrixSquareGeneral& gauss,
                                          VectorDouble &hgnc,
                                          VectorInt &ind_util,
                                          MatrixRectangular &bords_red,
                                          VectorDouble &ai_red,
                                          VectorDouble &grad_red,
                                          MatrixSquareGeneral& gauss_red,
                                          VectorInt &flag1,
                                          VectorInt &flag2,
                                          VectorDouble &temp)
{
  int n_imposs, ic, ipar, jpar, iparac, jparac, ipar2, iparac2;

  // Blanking out the arrays

  grad_red.fill(0.);
  gauss_red.fill(0.);
  bords_red.fill(0.);
  ai_red.fill(0.);

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
        gauss_red.setValue(iparac,jparac,gauss.getValue(ipar, jpar));
        jparac++;
      }
      iparac++;
    }

    // Reduce the array 'bords'
    st_update_bords(bords, ind_util, bords_red);

    /* Reduce the arrays bords and ai */
    for (ic = 0; ic < 2; ic++)
      for (ipar = iparac = 0; ipar < NPAR; ipar++)
      {
        if (!POSSIBLE(ipar)) continue;
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

      if (st_solve_hgnc(NPARAC + NCONT, grad_red, gauss_red, hgnc, 1))
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
**  Minimization under constraints
**
** \return  Error returned code
**
** \param[in]  nactive      Number of active constraints
** \param[in]  ind_util   List of retained constraint indices
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
static int st_establish_minimization(int nactive,
                                     VectorInt& ind_util,
                                     VectorInt &flag_active,
                                     MatrixRectangular &bords_red,
                                     VectorDouble &ai_red,
                                     VectorDouble &grad_red,
                                     MatrixSquareGeneral &gauss_red,
                                     int *lambda_neg,
                                     VectorDouble &hgnc,
                                     MatrixSquareGeneral &a,
                                     VectorDouble &b,
                                     VectorDouble &temp)
{
  int size, ic, iparac, jparac, iparac2, iecr;
  DECLARE_UNUSED(ind_util);

  /* Initialization */

  *lambda_neg = -1;
  size = NPARAC + NCONT + nactive;
  a.fill(0.);
  b.fill(0.);

  /* Fill the L.H.S. and R.H.S. matrices */

  for (iparac = 0; iparac < NPARAC + NCONT; iparac++)
  {
    b[iparac] = -grad_red[iparac];
    for (jparac = 0; jparac < NPARAC + NCONT; jparac++)
      a.setValue(iparac,jparac,gauss_red.getValue(iparac, jparac));
  }

  for (ic = iecr = 0; ic < 2; ic++)
    for (iparac = 0; iparac < NPARAC; iparac++)
    {
      if (!FLAG_ACTIVE(ic, iparac)) continue;
      b[NPARAC + NCONT + iecr] = bords_red.getValue(ic, iparac);
      for (jparac = 0; jparac < NPARAC; jparac++)
      {
        a.setValue(NPARAC + NCONT + iecr,jparac,AI_RED(ic, iparac, jparac));
        a.setValue(jparac,NPARAC + NCONT + iecr,AI_RED(ic, iparac, jparac));
      }
      iecr++;
    }

  /* Solve the system */

  if (st_solve_hgnc(size, b, a, temp, 0)) return (1);

  /* Store the final result */

  for (iparac = 0; iparac < NPARAC; iparac++)
    hgnc[iparac] = temp[iparac];
  for (ic = 0; ic < 2; ic++)
    for (iparac = 0; iparac < NPARAC; iparac++)
      if (FLAG_ACTIVE(ic, iparac)) hgnc[iparac] = bords_red.getValue(ic, iparac);

  /* Check for the first obviously negative coefficient */

  for (ic = iecr = iparac2 = 0; ic < 2; ic++)
    for (iparac = 0; iparac < NPARAC; iparac++, iparac2++)
    {
      if (!FLAG_ACTIVE(ic, iparac)) continue;
      if (temp[NPARAC + iecr] < -EPSILON9)
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
                     const MatrixRectangular& acont)
{
  double temp;
  int ipar, icont, iparac;

  if (NCONT <= 0) return;

  for (icont = 0; icont < NCONT; icont++)
  {
    temp = 0;
    for (iparac = ipar = 0; ipar < NPAR; ipar++)
      if (POSSIBLE(ipar)) temp += acont.getValue(ipar,icont) * hgnc[iparac++];

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
                                             MatrixRectangular &bords_red,
                                             VectorDouble &ai_red,
                                             VectorDouble &grad_red,
                                             MatrixSquareGeneral& gauss_red,
                                             MatrixRectangular &consts,
                                             VectorDouble &hgnc,
                                             VectorDouble &hgnadm,
                                             VectorInt &flag_active,
                                             VectorInt &flag_actaux,
                                             MatrixSquareGeneral &a,
                                             VectorDouble &b1,
                                             VectorDouble &b2,
                                             VectorDouble &b3,
                                             VectorDouble &temp,
                                             const MatrixRectangular& acont)
{
  int iparac, nactaux, sortie, nactive, lambda_neg;
  double min_adm_cur, min_adm_best;
  static int nitermax = 2000;

  // Clean out arrays
  hgnadm.fill(0.);

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

  matrix_product_safe(NPARAC2, NPARAC, 1, ai_red.data(), hgnc.data(), b1.data());
  st_minimum(ind_util, flag_actaux, bords_red, VectorDouble(), b1, hgnc, hgnadm);
  st_check(ind_util, hgnadm, acont);

  /* Calculate the constraints vector */

  nactive = st_define_constraints(0, bords_red, ai_red, hgnadm, consts,
                                  flag_active, temp);
  min_adm_best = st_essai(hgnadm, grad_red, gauss_red);
  if (VERBOSE_GQO && OptDbg::query(EDbg::CONVERGE))
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
      matrix_product_safe(NPARAC2, NPARAC, 1, ai_red.data(), hgnadm.data(), b1.data());
      matrix_product_safe(NPARAC2, NPARAC, 1, ai_red.data(), b3.data(), b2.data());
      st_minimum(ind_util, flag_actaux, bords_red, b1, b2, hgnc, hgnadm);
      st_check(ind_util, hgnadm, acont);

      nactive = st_define_constraints(1, bords_red, ai_red, hgnadm, consts,
                                      flag_active, temp);
      min_adm_cur = st_essai(hgnadm, grad_red, gauss_red);
      if (VERBOSE_GQO && OptDbg::query(EDbg::CONVERGE))
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
        if (VERBOSE_GQO && OptDbg::query(EDbg::CONVERGE))
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
                             MatrixRectangular &bords)
{
  double dloc, diff;
  int ipar;

  for (ipar = 0; ipar < NPAR; ipar++)
  {
    dloc = delta * scale[ipar];

    /* Lower bound */

    if (FFFF(lower[ipar]))
      bords.setValue(0,ipar,-dloc);
    else
    {
      diff = param[ipar] - lower[ipar];
      if (diff > dloc) diff = dloc;
      if (ABS(param[ipar] - diff) < dloc / 10.) diff /= 2.;
      bords.setValue(0,ipar,-diff);
    }

    /* Upper bound */

    if (FFFF(upper[ipar]))
      bords.setValue(1,ipar,dloc);
    else
    {
      diff = upper[ipar] - param[ipar];
      if (diff > dloc) diff = dloc;
      if (ABS(param[ipar] - diff) < dloc / 10.) diff /= 2.;
      bords.setValue(1,ipar,diff);
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

  if (!OptDbg::query(EDbg::CONVERGE)) return;
  mestitle(1, "Trajectory of parameters in Foxleg Algorithm");
  tab_prints(NULL, "Iteration");
  tab_prints(NULL, "Score");
  tab_prints(NULL, "Delta");
  for (ipar = 0; ipar < NPAR; ipar++)
  {
    (void) gslSPrintf(string, "Par-%d", ipar + 1);
    tab_prints(NULL, string);
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

  if (!OptDbg::query(EDbg::CONVERGE)) return;
  tab_printi(NULL, ITERATION);
  tab_printd(NULL, mscur);
  tab_printd(NULL, delta);
  for (ipar = 0; ipar < NPAR; ipar++)
    tab_printg(NULL, param[ipar]);
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
  if (mauto.getVerbose())
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
                                  const MatrixRectangular& acont,
                                  VectorDouble &tabexp,
                                  VectorDouble &tabwgt,
                                  MatrixRectangular &bords,
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
      if ((shift < 0 && ABS(bords.getValue(0,ipar)) < EPSILON9) ||
          (shift > 0 && ABS(bords.getValue(1,ipar)) < EPSILON9)) shift = 0.;

      /* Discard any move if constrained minimization */
      for (icont = 0; icont < NCONT; icont++)
        if (acont.getValue(ipar,icont) != 0.) shift = 0.;

      paramaux[ipar] = param[ipar] + shift;
      if (shift < bords.getValue(0,ipar) - EPSILON9 ||
          shift > bords.getValue(1,ipar) + EPSILON9)
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
 ** \remark  When not used, we must set: ncont=0, acont=empty()
 **
 *****************************************************************************/
int foxleg_f(int ndat,
             int npar,
             int ncont,
             const MatrixRectangular& acont,
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
  int  iparac;
  double msaux;

  /* Preliminary checks */

  double arret = 0.;
  double delta = mauto.getInitdelta();
  NDAT = ndat;
  NPAR = NPARAC = npar;
  NCONT = ncont;
  NPAR2 = NPARAC2 = npar * 2;
  NPCT  = NPAR + NCONT;
  NPCT2 = NPAR2 + NCONT;
  FUNC_EVALUATE = func_evaluate;
  USE_EIGEN_LIBRARY = mauto.isUseEigenLibrary();

  /* Core allocation */

  VectorInt ind_util(NPCT, 0);
  VectorInt flag_active(NPAR2, 0);
  VectorInt flag_actaux(NPAR2, 0);

  VectorDouble b1(NPCT2, 0.);
  VectorDouble b2(NPAR2, 0.);
  VectorDouble b3(NPAR2, 0.);
  VectorDouble temp(NPCT2, 0.);
  VectorDouble param1(NPAR, 0.);
  VectorDouble param2(NPAR, 0.);
  VectorDouble grad(NPCT, 0.);
  VectorDouble grad_red(NPCT, 0.);
  VectorDouble hgn(NPAR, 0.);
  VectorDouble hgnc(NPCT, 0.);
  VectorDouble hgnadm(NPCT, 0.);
  VectorDouble paramaux(NPAR, 0.);
  VectorDouble residuals(NDAT, 0.);
  VectorDouble tabmod1(NDAT, 0.);
  VectorDouble tabmod2(NDAT, 0.);
  VectorDouble ai(NPAR * NPAR2, 0.);
  VectorDouble ai_red(NPAR * NPAR2, 0.);

  MatrixSquareGeneral a(NPCT2);
  MatrixSquareGeneral gauss(NPCT);
  MatrixSquareGeneral gauss_red(NPCT);
  MatrixRectangular Jr(NDAT, NPAR);
  MatrixRectangular consts(2,NPAR);
  MatrixRectangular bords(2,NPAR);
  MatrixRectangular bords_red(2,NPAR);

  /* Preliminary check */

  if (st_check_param(param, lower, upper)) return 1;

  /* Initializations */

  for (int ipar = 0; ipar < NPAR; ipar++)
    paramaux[ipar] = param[ipar];

  st_constraints_init(ind_util, ai);

  /* Calculate the gradient */
  double ms0 = st_residuals(param, tabexp, tabwgt, tabmod1, residuals);
  double mscur = ms0;
  if (st_calcul0(param, lower, upper, scale, acont, tabwgt, residuals, Jr, grad,
                 gauss, hgnc, param1, param2, tabmod1, tabmod2)) return 1;

  st_foxleg_debug_title();

  /***********************/
  /* Iterative procedure */
  /***********************/

  ITERATION = 0;
  bool flag_moved = true;
  bool flag_cont = (mauto.getMaxiter() <= 0) ? false : true;
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
      if (st_suppress_unused_constraints(bords, ai, grad, gauss, hgnc,
                                         ind_util, bords_red, ai_red, grad_red,
                                         gauss_red, flag_active, flag_actaux,
                                         temp)) break;
    }
    else
    {
      st_update_bords(bords, ind_util, bords_red);
    }

    /* Perform the minimization under constraints */

    if (st_minimization_under_constraints(ind_util, bords_red, ai_red, grad_red,
                                          gauss_red, consts, hgnc, hgnadm,
                                          flag_active, flag_actaux, a, b1, b2,
                                          b3, temp, acont))
    {
      if (OptDbg::query(EDbg::CONVERGE))
      {
        messerr("Convergence not reached in minimization under constraints");
        messerr("The process is resumed in the final minimization status");
      }
    }

    /* Check for an early stop */

    if (ITERATION >= mauto.getMaxiter() || delta < mauto.getEpsdelta())
      flag_cont = false;
    arret = 0.;
    for (iparac = 0; iparac < NPARAC; iparac++)
      arret += ABS(hgnc[iparac]);
    arret /= (ms0 - mscur + mauto.getEpsdelta());
    if (arret < mauto.getTolstop()) flag_cont = false;

    /* Update values for the next iteration */

    iparac = 0;
    for (int ipar = 0; ipar < NPAR; ipar++)
    {
      hgn[ipar] = 0.;
      paramaux[ipar] = param[ipar];
      if (!POSSIBLE(ipar)) continue;
      hgn[ipar] = hgnadm[iparac++];
      paramaux[ipar] += hgn[ipar];
    }

    double denom = st_essai(hgnadm, grad_red, gauss_red);
    if (isZero(denom)) goto label_ok;
    if (denom > 0)
      st_linear_interpolate(mscur, param, acont, tabexp, tabwgt, bords, grad,
                            &msaux, paramaux, residuals, tabmod1);
    else
      msaux = st_residuals(paramaux, tabexp, tabwgt, tabmod1, residuals);
    double rho = (msaux - mscur) / denom;

    if (msaux < mscur)
    {
      SOUSITER++;
      mscur = msaux;
      for (int ipar = 0; ipar < NPAR; ipar++)
        param[ipar] = paramaux[ipar];
      if (st_calcul0(param, lower, upper, scale, acont, tabwgt, residuals, Jr,
                     grad, gauss, hgnc, param1, param2, tabmod1,
                     tabmod2)) return 1;
      st_constraints_init(ind_util, ai);
      flag_moved = true;
      if (denom < 0 && rho > 0.75)
        delta = MAX(delta, 3. * st_norm_hgn(hgn, scale));
    }
    else
      flag_moved = false;
    if (rho < 0.25 || denom > 0) delta /= 2.;
  }

  /* Result printout */

  label_ok:
  if (ITERATION >= mauto.getMaxiter())
  {
    if (OptDbg::query(EDbg::CONVERGE)) messerr("Convergence has not been reached");
    return -1;
  }
  else
  {
    return 0;
  }

  set_keypair("Foxleg Value", 1, 1, 1, &mscur);
  st_foxleg_score(mauto, mscur, delta, arret);

  return 0;
}

/****************************************************************************/
/*!
 **  Add constraints to the Option_AutoFit structure
 **
 ** \return Error return code
 **
 ** \param[in]  constraints  Constraints structure
 ** \param[in]  constantSill Constant value for the Sill as a constraint
 **
 *****************************************************************************/
int add_sill_constraints(Constraints& constraints, double constantSill)
{
  constraints.setConstantSillValue(constantSill);

  return (0);
}

/****************************************************************************/
/*!
 **  Add constraints (all equal to 1) to the Option_AutoFit structure
 **
 ** \return Error return code
 **
 ** \param[in]  constraints   Constraints structure
 **
 *****************************************************************************/
int add_unit_sill_constraints(Constraints& constraints)
{
  constraints.setConstantSillValue(1.);
  return (0);
}
