/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_f.h"
#include "geoslib_enum.h"
#include "geoslib_old_f.h"
#include "Anamorphosis/AnamDiscreteDD.hpp"
#include "Anamorphosis/AnamDiscreteIR.hpp"
#include "Anamorphosis/AnamEmpirical.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Anamorphosis/AnamUser.hpp"
#include "Variogram/Vario.hpp"
#include "Polynomials/Hermite.hpp"
#include "Polynomials/MonteCarlo.hpp"
#include "Stats/Selectivity.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Db/Db.hpp"

#include <math.h>

/*! \cond */
#define EPS_TON    0.001

#define QT_EST    0
#define QT_STD    1

#define CT(i,j)         (ct[(i)*nclass+(j)])
#define CQ(i,j)         (cq[(i)*nclass+(j)])
#define CB(i,j)         (cb[(i)*nclass+(j)])
#define TAU(ip,jp)                (tau[(ip) * nbpoly + (jp)])
#define DD(ih,jh)                 (dd[(ih) * nbpoly + (jh)])
#define QT_VARS(i,j)              (qt_vars[(i) + 2 * (j)])
#define QT_FLAG(j)                (QT_VARS(QT_EST,j) > 0 || \
                                   QT_VARS(QT_STD,j) > 0)
/*! \endcond */

/*****************************************************************************/
/*!
 **  Interpolate the QT within an interval
 **
 ** \param[in]  zval     Cutoff value
 ** \param[in]  zi0      Lower cutoff of the interval
 ** \param[in]  zi1      Upper cutoff of the interval
 ** \param[in]  ti0      Lower tonnage of the interval
 ** \param[in]  ti1      Upper tonnage of the interval
 ** \param[in]  qi0      Lower metal quantity of the interval
 ** \param[in]  qi1      Upper metal quantity of the interval
 **
 ** \param[out] tval     Tonnage for the current cutoff
 ** \param[out] qval     Metal quantity for the current cutoff
 **
 *****************************************************************************/
static void st_interpolate_interval(double zval,
                                    double zi0,
                                    double zi1,
                                    double ti0,
                                    double ti1,
                                    double qi0,
                                    double qi1,
                                    double *tval,
                                    double *qval)
{
  double dzi, dti, u, aa0, zmoy;
  static double tol = 1.e-3;

  dzi = zi1 - zi0;
  dti = ti1 - ti0;
  zmoy = (qi1 - qi0) / (ti1 - ti0);
  aa0 = (zi1 - zmoy) / (zmoy - zi0);

  if (ABS(zval - zi0) < tol)
  {
    (*tval) = ti0;
    (*qval) = qi0;
    return;
  }

  if (ABS(zval - zi1) < tol)
  {
    (*tval) = ti1;
    (*qval) = qi1;
    return;
  }

  u = (zval - zi0) / dzi;
  (*tval) = (u <= 0.) ? ti0 : ti0 + dti * pow(u, 1. / aa0);
  (*qval) = (u <= 0.) ? qi0 :
      qi0 + zi0 * ((*tval) - ti0)
      + dzi * dti * pow(u, 1. + 1. / aa0) / (1. + aa0);
}

/*****************************************************************************/
/*!
 **  Interpolate the QT curves (Local estimation)
 **
 ** \param[in]  z_max    Maximum grade value (if defined)
 ** \param[in]  zcutmine Array of the requested cutoffs
 ** \param[in]  calest   Selectivity
 **
 ** \param[out] calcut   Interpolated Selectivity
 **
 *****************************************************************************/
static void st_interpolate_qt_local(double z_max,
                                    double *zcutmine,
                                    Selectivity& calest,
                                    Selectivity& calcut)
{
  double tval, qval;

  int nclass = calest.getNClass();
  int ncutmine = calcut.getNClass();
  VectorDouble zz(nclass + 2);
  VectorDouble TT(nclass + 2);
  VectorDouble QQ(nclass + 2);

  /* Load arrays */

  int ncleff = 1;
  TT[0] = QQ[0] = 0.;
  for (int iclass = 0; iclass < nclass; iclass++)
  {
    int jclass = nclass - iclass - 1;
    if (calest.getTest(jclass) <= TT[ncleff - 1]) continue;
    TT[ncleff] = calest.getTest(jclass);
    QQ[ncleff] = calest.getQest(jclass);
    ncleff++;
  }
  zz[0] = z_max;
  for (int iclass = 0; iclass < ncleff - 1; iclass++)
    zz[iclass + 1] = (QQ[iclass + 2] - QQ[iclass]) / (TT[iclass + 2] - TT[iclass]);
  zz[ncleff - 1] = 0.;
  if (FFFF(z_max)) zz[0] = 2 * zz[1];

  for (int icut = 0; icut < ncutmine; icut++)
  {
    double zval = zcutmine[icut];
    calcut.setZcut(icut, zval);

    /* Find interval [zz[iclass]; zz[iclass+1]] to which cutoffs belongs */

    int iclass = -1;
    for (int jclass = 0; jclass < ncleff && iclass < 0; jclass++)
      if ((zval - zz[jclass]) * (zval - zz[jclass + 1]) <= 0) iclass = jclass;

    /* Assuming that cutoffs belongs to the interval the class 'iclass' */

    double zi0 = zz[iclass];
    double zi1 = (iclass + 1 > ncleff - 1) ? 0. : zz[iclass + 1];
    double ti0 = TT[iclass];
    double ti1 = (iclass + 1 > ncleff - 1) ? 0. : TT[iclass + 1];
    double qi0 = QQ[iclass];
    double qi1 = QQ[iclass + 1];
    st_interpolate_interval(zval, zi0, zi1, ti0, ti1, qi0, qi1,
                            &tval, &qval);
    calcut.setTest(icut, tval);
    calcut.setQest(icut, qval);
  }
  return;
}

/*****************************************************************************/
/*!
 **  Interpolate the QT curves (Global estimation)
 **
 ** \param[in]  zcutmine Array of the requested cutoffs
 ** \param[in]  calest   Selectivity
 **
 ** \param[out] calcut   Interpolated Selectivity
 **
 *****************************************************************************/
static void st_interpolate_qt_global(double *zcutmine,
                                     Selectivity& calest,
                                     Selectivity& calcut)
{
  double tval, qval;

  int nclass = calest.getNClass();
  int ncutmine = calcut.getNClass();

  for (int icut = 0; icut < ncutmine; icut++)
  {
    double zval = zcutmine[icut];
    calcut.setZcut(icut, zval);

    /* Find interval [zmoy[iclass]; zmoy[iclass+1]] to which cutoffs belongs */

    int iclass = -1;
    for (int jclass = 0; jclass < nclass - 1 && iclass < 0; jclass++)
    {
      double valmin = MIN(calest.getZcut(jclass), calest.getZcut(jclass+1));
      double valmax = MAX(calest.getZcut(jclass), calest.getZcut(jclass+1));
      if (zval >= valmin && zval <= valmax) iclass = jclass;
    }

    if (iclass >= 0 && iclass < nclass)
    {

      /* Assuming that cutoffs belongs to the interval the class 'iclass' */

      double zi0 = calest.getZcut(iclass);
      double zi1 = (iclass + 1 > nclass - 1) ? 0. : calest.getZcut(iclass + 1);
      double ti0 = calest.getTest(iclass);
      double ti1 = (iclass + 1 > nclass - 1) ? 0. : calest.getTest(iclass + 1);
      double qi0 = calest.getQest(iclass);
      double qi1 = calest.getQest(iclass + 1);
      st_interpolate_interval(zval, zi0, zi1, ti0, ti1, qi0, qi1,
                              &tval, &qval);
      calcut.setTest(icut, tval);
      calcut.setQest(icut, qval);
    }
    else
    {
      calcut.setTest(icut, 0.);
      calcut.setQest(icut, 0.);
    }
  }
  return;
}

/****************************************************************************/
/*!
 **  Store the local results of the recovery
 **
 ** \param[in]  db          Db structure containing the factors (Z-locators)
 ** \param[in]  iech0       Rank of the target sample
 ** \param[in]  iptr        Rank for storing the results
 ** \param[in]  ncode       Number of stored results
 ** \param[in]  codes       Array of codes for stored results
 ** \param[in]  qt_vars     Array with the number of output variables
 ** \param[in]  zestim      Estimated grade
 ** \param[in]  zstdev      St. dev.
 ** \param[in]  calest      Selectivity
 **
 *****************************************************************************/
static void st_recovery_local(Db *db,
                              int iech0,
                              int iptr,
                              int ncode,
                              int *codes,
                              int *qt_vars,
                              double zestim,
                              double zstdev,
                              const Selectivity& calest)
{
  int jptr = iptr;
  int nclass = calest.getNClass();

  /* Store the recovered grade */

  if (codes[0] == 0)
  {
    if (QT_VARS(QT_EST,ANAM_QT_Z) > 0) db->setArray(iech0, jptr++, zestim);
    if (QT_VARS(QT_STD,ANAM_QT_Z) > 0) db->setArray(iech0, jptr++, zstdev);
  }

  /* Loop on the recovery functions */

  for (int icode = 0; icode < ncode; icode++)
  {

    /* Loop on the cutoff classes */

    for (int iclass = 0; iclass < nclass; iclass++)
    {
      double tval = calest.getTest(iclass);
      double qval = calest.getQest(iclass);
      double bval = calest.getBest(iclass);
      double mval = calest.getMest(iclass);
      double tstd = calest.getTstd(iclass);
      double qstd = calest.getQstd(iclass);

      switch (codes[icode])
      {
        case 1: /* Tonnage */
          if (QT_VARS(QT_EST,ANAM_QT_T) > 0) db->setArray(iech0, jptr++, tval);
          if (QT_VARS(QT_STD,ANAM_QT_T) > 0) db->setArray(iech0, jptr++, tstd);
          break;

        case 2: /* Metal Quantity */
          if (QT_VARS(QT_EST,ANAM_QT_Q) > 0) db->setArray(iech0, jptr++, qval);
          if (QT_VARS(QT_STD,ANAM_QT_Q) > 0) db->setArray(iech0, jptr++, qstd);
          break;

        case 3: /* Conventional Benefit */
          if (QT_VARS(QT_EST,ANAM_QT_B) > 0) db->setArray(iech0, jptr++, bval);
          break;

        case 4: /* Average recovered grade */
          if (QT_VARS(QT_EST,ANAM_QT_M) > 0) db->setArray(iech0, jptr++, mval);
          break;
      }
    }
  }
}

/****************************************************************************/
/*!
 **  Calculate the theoretical grade tonnage value
 **
 ** \return  Array f results (Dimension: 7 * nclass)
 **
 ** \param[in] anam         AAnam structure to be updated
 ** \param[in] nclass       Number of classes
 ** \param[in] zcut         Array of cutoffs
 ** \param[in] flag_correct 1 if Tonnage order relationship must be corrected
 ** \param[in] verbose      Verbose flag
 **
 ** \remark In the case of Discrete Anamorphosis, the number of classes
 ** \remark is defined by the number of cutoffs
 **
 *****************************************************************************/
Selectivity anam_selectivity(AAnam *anam,
                             int nclass,
                             VectorDouble zcut,
                             int flag_correct,
                             int verbose)
{
  Selectivity calest(nclass);

  /* Dispatch according to the anamorphosis */

  if (anam->getType() == EAnam::HERMITIAN)
  {
    AnamHermite *anam_hermite = dynamic_cast<AnamHermite*>(anam);
    calest = anam_hermite->calculateSelectivity(zcut);
  }
  else if (anam->getType() == EAnam::DISCRETE_DD)
  {
    AnamDiscreteDD *anam_discrete_DD = dynamic_cast<AnamDiscreteDD*>(anam);
    if (nclass != anam_discrete_DD->getNClass())
    {
      messerr(
          "Argument 'nclass' (%d) should be equal to the number of classes (%d)",
          nclass, anam_discrete_DD->getNClass());
      return calest;
    }
    calest = anam_discrete_DD->calculateSelectivity(flag_correct);
  }
  else if (anam->getType() == EAnam::DISCRETE_IR)
  {
    AnamDiscreteIR *anam_discrete_IR = dynamic_cast<AnamDiscreteIR*>(anam);
    if (nclass != anam_discrete_IR->getNClass())
    {
      messerr(
          "Argument 'nclass' (%d) should be equal to the number of classes (%d)",
          nclass, anam_discrete_IR->getNClass());
      return calest;
    }
    calest = anam_discrete_IR->calculateSelectivity(flag_correct);

  }
  else
  {
    messerr("This function is not programmed for this Anamorphosis");
    return calest;
  }

  if (verbose) calest.dumpGini();
  return calest;
}

/*****************************************************************************/
/*!
 **  Transform a point anamorphosis into a block anamorphosis
 **
 ** \return  Error return code
 **
 ** \param[in]  anam        Point anamorphosis -> Block anamorphosis [out]
 ** \param[in]  verbose     Verbose option
 ** \param[in]  cvv         Block variance
 ** \param[in]  coeff       Coefficient of change of support
 ** \param[in]  mu          Additional coefficient for Discrete case
 **
 ** \remark If 'coeff' is provided, it is used directly ('cvv' is ignored)
 ** \remark Otherwise, it is derived from 'cvv'
 **
 *****************************************************************************/
int anam_point_to_block(AAnam *anam,
                        int verbose,
                        double cvv,
                        double coeff,
                        double mu)
{
  if (anam == nullptr) return (1);
  double r_coef = 0.;
  AnamHermite *anam_hermite = dynamic_cast<AnamHermite*>(anam);
  AnamDiscreteDD *anam_discrete_DD = dynamic_cast<AnamDiscreteDD*>(anam);
  AnamDiscreteIR *anam_discrete_IR = dynamic_cast<AnamDiscreteIR*>(anam);

  /* Preliminary check */

  if (!FFFF(coeff) && (coeff < 0 || coeff > 1.))
  {
    messerr("Change of support coefficient (%lf) must lie between 0 and 1.",
            coeff);
    return 1;
  }

  /* Dispatch according to the anamorphosis type */

  switch (anam->getType().toEnum())
  {
    case EAnam::E_HERMITIAN:

      if (FFFF(coeff))
      {
        r_coef = anam->calculateR(cvv, 2.);
        if (verbose)
        {
          mestitle(1, "Calculation of the Change of Support Coefficient");
          message("Average Block covariance      = %lf\n", cvv);
          message("Change of support coefficient = %lf\n",
                  anam_hermite->getRCoef());
        }
      }
      else
        r_coef = coeff;
      break;

    case EAnam::E_DISCRETE_DD:
      anam_discrete_DD->setMu(mu);
      if (FFFF(coeff))
      {
        r_coef = anam->calculateR(cvv, 2.);
        if (verbose)
        {
          mestitle(1, "Calculation of the Change of Support Coefficient");
          message("Point Variance                = %lf\n",
                  anam_discrete_DD->getVariance());
          message("Average Block covariance      = %lf\n", cvv);
          message("Coefficient mu                = %lf\n",
                  anam_discrete_DD->getMu());
          message("Change of support coefficient = %lf\n",
                  anam_discrete_DD->getSCoef());
        }
      }
      else
        r_coef = coeff;
      break;

    case EAnam::E_DISCRETE_IR:
      if (FFFF(coeff))
      {
        r_coef = anam->calculateR(cvv, 2.);
        if (verbose)
        {
          mestitle(1, "Calculation of the Change of Support Coefficient");
          message("Average Block covariance      = %lf\n", cvv);
          message("Change of support coefficient = %lf\n",
                  anam_discrete_IR->getRCoef());
        }
      }
      else
        r_coef = coeff;
      break;

    default:
      messerr("The change of support is not defined for this Anamorphosis");
      return 1;
  }

  /* Update the Point Anamorphosis into Block Anamorphosis */

  anam->updatePointToBlock(r_coef);

  return 0;
}

/*****************************************************************************/
/*!
 **  Check if a sample must be considered or not
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  iech         Rank of the target sample
 ** \param[in]  nb_est       Number of columns for factor estimation
 ** \param[in]  cols_est     Array of columns for factor estimation
 ** \param[in]  nb_std       Number of columns for factor st. dev.
 ** \param[in]  cols_std     Array of columns for factor st. dev.
 **
 *****************************************************************************/
static int st_is_sample_skipped(Db *db,
                                int iech,
                                int nb_est,
                                int *cols_est,
                                int nb_std,
                                int *cols_std)
{
  double value;

  if (!db->isActive(iech)) return (1);

  if (cols_est != nullptr) for (int ivar = 0; ivar < nb_est; ivar++)
  {
    value = db->getArray(iech, cols_est[ivar]);
    if (FFFF(value)) return (1);
  }

  if (cols_std != nullptr) for (int ivar = 0; ivar < nb_std; ivar++)
  {
    value = db->getArray(iech, cols_std[ivar]);
    if (FFFF(value)) return (1);
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Print the contents of the qtvars structure
 **
 ** \param[in]  title        Title
 ** \param[in]  type         1 for estimation; 2 for stdev
 ** \param[in]  number       Number of cutoffs
 **
 *****************************************************************************/
static void st_print_qtvars(const char *title, int type, int number)
{
  message("- %s", title);
  if (type == 1)
    message(" (Estimation)");
  else
    message(" (St. Deviation)");
  message(": %d\n", number);
}

/*****************************************************************************/
/*!
 **  Check that 'ncut' is positive
 **
 ** \param[in]  ncut         Number of cutoffs
 **
 *****************************************************************************/
static bool isNcutValid(int ncut)
{
  if (ncut <= 0)
  {
    messerr("The computing option requires Cutoffs to be defiend");
    return false;
  }
  return true;
}

/*****************************************************************************/
/*!
 **  Check that 'proba' is defined
 **
 ** \param[in]  proba         Probability
 **
 *****************************************************************************/
static int st_check_proba(double proba)
{
  if (FFFF(proba))
  {
    messerr("The computing option requires Proba to be defined");
    return (1);
  }
  if (proba < 0 || proba > 1)
  {
    messerr("The computing option requires Proba to lie in [0,1]");
    return (1);
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Analyze the contents of the codes
 **
 ** \returns The number of different variables to be calculated
 **
 ** \param[in]  verbose      Verbose flag
 ** \param[in]  ncode        Number of operators
 ** \param[in]  codes        Array of codes for stored results
 ** \param[in]  nb_est       Number of columns for factor estimation
 ** \param[in]  nb_std       Number of columns for factor st. dev.
 ** \param[in]  ncut         Number of cutoffs
 ** \param[in]  proba        Probability value
 ** \param[in]  flag_inter   QT must be interpolated
 **
 ** \param[out] qt_vars      Array with the number of output variables
 **
 ** \remark When the number of cutoff is zero, the flag of T, Q, B and M
 ** \remark are set to 0
 ** \remark When QT are interpolated, no variance can be calculated
 **
 ** \remark When the resulting number of variables is zero, an error
 ** \remark message is issued
 **
 *****************************************************************************/
static int st_code_analyze(int verbose,
                           int ncode,
                           int *codes,
                           int nb_est,
                           int nb_std,
                           int ncut,
                           double proba,
                           int flag_inter,
                           int *qt_vars)
{
  int ntotal, flag_est, flag_std;

  /* Initializations */

  flag_est = nb_est > 0;
  flag_std = nb_std > 0 && !flag_inter;
  for (int i = 0; i < 2 * ANAM_N_QT; i++)
    qt_vars[i] = 0;
  ut_sort_int(0, ncode, NULL, codes);

  // Optional printout (title)

  if (verbose) mestitle(1, "List of options");

  /* Check for the presence of other codes */

  for (int icode = 0; icode < ncode; icode++)
  {
    switch (codes[icode])
    {
      case ANAM_QT_Z:
        if (flag_est)
        {
          QT_VARS(QT_EST,ANAM_QT_Z) = 1;
          if (verbose) st_print_qtvars("Average", 1, 1);
        }
        if (flag_std)
        {
          QT_VARS(QT_STD,ANAM_QT_Z) = 1;
          if (verbose) st_print_qtvars("Average", 2, 1);
        }
        break;

      case ANAM_QT_T:
        if (!isNcutValid(ncut)) return (0);
        if (flag_est)
        {
          QT_VARS(QT_EST,ANAM_QT_T) = ncut;
          if (verbose) st_print_qtvars("Tonnage", 1, ncut);
        }
        if (flag_std)
        {
          QT_VARS(QT_STD,ANAM_QT_T) = ncut;
          if (verbose) st_print_qtvars("Tonnage", 2, ncut);
        }
        break;

      case ANAM_QT_Q:
        if (!isNcutValid(ncut)) return (0);
        if (flag_est)
        {
          QT_VARS(QT_EST,ANAM_QT_Q) = ncut;
          if (verbose) st_print_qtvars("Metal Quantity", 1, ncut);
        }
        if (flag_std)
        {
          QT_VARS(QT_STD,ANAM_QT_Q) = ncut;
          if (verbose) st_print_qtvars("Metal Quantity", 2, ncut);
        }
        break;

      case ANAM_QT_B:
        if (!isNcutValid(ncut)) return (0);
        if (flag_est)
        {
          QT_VARS(QT_EST,ANAM_QT_B) = ncut;
          if (verbose) st_print_qtvars("Conventional Benefit", 1, ncut);
        }
        break;

      case ANAM_QT_M:
        if (!isNcutValid(ncut)) return (0);
        if (flag_est)
        {
          QT_VARS(QT_EST,ANAM_QT_M) = ncut;
          if (verbose) st_print_qtvars("Average Metal", 1, ncut);
        }
        break;

      case ANAM_QT_PROBA:
        if (!isNcutValid(ncut)) return (0);
        if (flag_est)
        {
          QT_VARS(QT_EST,ANAM_QT_PROBA) = ncut;
          if (verbose) st_print_qtvars("Probability", 1, ncut);
        }
        break;

      case ANAM_QT_QUANT:
        if (st_check_proba(proba)) return (0);
        if (flag_est)
        {
          QT_VARS(QT_EST,ANAM_QT_QUANT) = 1;
          if (verbose) st_print_qtvars("Quantile", 1, 1);
        }
        break;
    }
  }

  /* Count the total number of variables */

  ntotal = 0;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < ANAM_N_QT; j++)
      ntotal += QT_VARS(i, j);

  if (ntotal <= 0)
  {
    messerr("The number of variables calculated is zero");
    messerr("Check the recovery function (the number of cutoffs is %d)", ncut);
  }

  return (ntotal);
}

/*****************************************************************************/
/*!
 **  Calculate Experimental Grade-Tonnage curves from factors
 **  Case of Hermite Anamorphosis
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  anam         Point anamorphosis
 ** \param[in]  ncutmine     Number of required cutoffs
 ** \param[in]  cutmine      Array of the requested cutoffs
 ** \param[in]  nb_est       Number of columns for factor estimation
 ** \param[in]  cols_est     Array of columns for factor estimation
 ** \param[in]  nb_std       Number of columns for factor st. dev.
 ** \param[in]  cols_std     Array of columns for factor st. dev.
 ** \param[in]  iptr         Rank for storing the results
 ** \param[in]  ncode        Number of stored results
 ** \param[in]  codes        Array of codes for stored results
 ** \param[in]  qt_vars      Array of variables to be calculated
 **
 ** \param[out] calest       Selectivity structure
 **
 *****************************************************************************/
static int st_anam_factor2qt_hermitian(Db *db,
                                       AAnam *anam,
                                       int ncutmine,
                                       double *cutmine,
                                       int nb_est,
                                       int *cols_est,
                                       int nb_std,
                                       int *cols_std,
                                       int iptr,
                                       int ncode,
                                       int *codes,
                                       int *qt_vars,
                                       Selectivity& calest)
{
  AnamHermite *anam_hermite = dynamic_cast<AnamHermite*>(anam);
  anam_hermite->setFlagBound(1);
  int nbpoly = anam_hermite->getNbPoly();
  bool need_T = QT_FLAG(ANAM_QT_T) || QT_FLAG(ANAM_QT_B) ||
      QT_FLAG(ANAM_QT_M) || QT_FLAG(ANAM_QT_PROBA);
  bool need_Q = QT_FLAG(ANAM_QT_Q) || QT_FLAG(ANAM_QT_B) || QT_FLAG(ANAM_QT_M);

  /* Loop on the samples */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (st_is_sample_skipped(db, iech, nb_est, cols_est, nb_std, cols_std))
      continue;

    /* Z: Estimation */

    double zestim = 0.;
    if (QT_VARS(QT_EST,ANAM_QT_Z) > 0)
    {
      double total = anam_hermite->getPsiHn(0);
      for (int ivar = 0; ivar < nb_est; ivar++)
      {
        double value = db->getArray(iech, cols_est[ivar]);
        double coeff = anam_hermite->getPsiHn(ivar + 1);
        total += coeff * value;
      }
      zestim = total;
    }

    /* Z: Standard Deviation */

    double zstdev = 0.;
    if (QT_VARS(QT_STD,ANAM_QT_Z) > 0)
    {
      double total = 0.;
      for (int ivar = 0; ivar < nb_std; ivar++)
      {
        double value = db->getArray(iech, cols_std[ivar]);
        double coeff = anam_hermite->getPsiHn(ivar + 1);
        total += coeff * coeff * value;
      }
      zstdev = sqrt(total);
    }

    /* Loop on the cutoffs */

    for (int icut = 0; icut < ncutmine; icut++)
    {
      double yc = anam_hermite->RawToTransformValue(cutmine[icut]);
      if (need_T)
      {
        VectorDouble s_cc = hermiteCoefIndicator(yc, nbpoly);

        /* Tonnage estimation */

        if (QT_VARS(QT_EST,ANAM_QT_T) > 0)
        {
          double total = s_cc[0];
          for (int ivar = 0; ivar < nb_est; ivar++)
          {
            double value = db->getArray(iech, cols_est[ivar]);
            total += s_cc[ivar + 1] * value;
          }
          calest.setTest(icut, total);
        }

        /* Tonnage: Standard Deviation */

        if (QT_VARS(QT_STD,ANAM_QT_T) > 0)
        {
          double total = 0.;
          for (int ivar = 0; ivar < nb_std; ivar++)
          {
            double value = db->getArray(iech, cols_std[ivar]);
            total += s_cc[ivar + 1] * s_cc[ivar + 1] * value;
          }
          calest.setTstd(icut, sqrt(total));
        }
      }

      if (need_Q)
      {
        MatrixSquareGeneral TAU = hermiteIncompleteIntegral(yc, nbpoly);

        /* Metal Quantity: Estimation */

        if (QT_VARS(QT_EST,ANAM_QT_Q) > 0)
        {
          double total = 0.;
          for (int ivar = 0; ivar < nb_est; ivar++)
          {
            double value = db->getArray(iech, cols_est[ivar]);
            double fn = 0.;
            for (int jvar = 0; jvar < nbpoly; jvar++)
            {
              double coeff = anam_hermite->getPsiHn(jvar);
              fn += coeff * TAU.getValue(ivar, jvar);
            }
            total += fn * value;
          }
          calest.setQest(icut, total);
        }

        /* Metal Quantity: Standard Deviation */

        if (QT_VARS(QT_STD,ANAM_QT_Q) > 0)
        {
          double total = 0.;
          for (int ivar = 0; ivar < nb_std; ivar++)
          {
            double value = db->getArray(iech, cols_std[ivar]);
            double fn = 0.;
            for (int jvar = 0; jvar < nbpoly; jvar++)
            {
              double coeff = anam_hermite->getPsiHn(jvar);
              fn += coeff * TAU.getValue(ivar, jvar);
            }
            total += fn * fn * value;
          }
          calest.setQstd(icut, sqrt(total));
        }
      }
    }

    /* Storage */

    calest.calculateBenefitGrade();
    st_recovery_local(db, iech, iptr, ncode, codes, qt_vars, zestim, zstdev,
                      calest);
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Calculate Experimental Grade-Tonnage curves from factors
 **  Case of Discrete Diffusion
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  anam         Point anamorphosis
 ** \param[in]  ncutmine     Number of required cutoffs
 ** \param[in]  cutmine      Array of the requested cutoffs
 ** \param[in]  z_max        Maximum grade array (only for QT interpolation)
 ** \param[in]  flag_correct 1 if Tonnage order relationship must be corrected
 ** \param[in]  nb_est       Number of columns for factor estimation
 ** \param[in]  cols_est     Array of columns for factor estimation
 ** \param[in]  nb_std       Number of columns for factor st. dev.
 ** \param[in]  cols_std     Array of columns for factor st. dev.
 ** \param[in]  iptr         Rank for storing the results
 ** \param[in]  ncode        Number of stored results
 ** \param[in]  codes        Array of codes for stored results
 ** \param[in]  qt_vars      Array of variables to be calculated
 **
 ** \param[out] calest       Selectivity structure
 ** \param[out] calcut       Translated Selectivity
 **
 *****************************************************************************/
static int st_anam_factor2qt_discrete_DD(Db *db,
                                         AAnam *anam,
                                         int ncutmine,
                                         double *cutmine,
                                         double z_max,
                                         int flag_correct,
                                         int nb_est,
                                         int *cols_est,
                                         int nb_std,
                                         int *cols_std,
                                         int iptr,
                                         int ncode,
                                         int *codes,
                                         int *qt_vars,
                                         Selectivity& calest,
                                         Selectivity& calcut)
{
  AnamDiscreteDD *anam_discrete_DD = dynamic_cast<AnamDiscreteDD*>(anam);
  int nclass = anam_discrete_DD->getNClass();
  int nech = db->getSampleNumber();

  /* Modeling the diffusion process */

  VectorDouble chi = anam_discrete_DD->factors_mod();
  if (chi.empty()) return 1;
  VectorDouble ct = anam_discrete_DD->chi2I(chi, 1);
  VectorDouble cq = anam_discrete_DD->chi2I(chi, 2);
  VectorDouble cb = anam_discrete_DD->chi2I(chi, 3);
  for (int iclass = 0; iclass < nclass; iclass++)
    calest.setZcut(iclass, anam_discrete_DD->getDDStatZmoy(iclass));

  /* Calculate the Recovery Functions from the factors */

  for (int iech = 0; iech < nech; iech++)
  {
    if (st_is_sample_skipped(db, iech, nb_est, cols_est, nb_std, cols_std))
      continue;

    /* Tonnage: Estimation */

    for (int iclass = 0; iclass < nclass; iclass++)
    {
      double total = CT(iclass, 0);
      for (int ivar = 0; ivar < nb_est; ivar++)
      {
        double value = db->getArray(iech, cols_est[ivar]);
        total += value * CT(iclass, ivar + 1);
      }
      calest.setTest(iclass, total);
    }

    /* Correct Order relationship for Tonnages */

    if (flag_correct) calest.correctTonnageOrder();

    /* Tonnage: Standard Deviation */

    if (QT_VARS(QT_STD,ANAM_QT_T) > 0)
    {
      for (int iclass = 0; iclass < nclass; iclass++)
      {
        double total = 0.;
        for (int ivar = 0; ivar < nclass - 1; ivar++)
        {
          double value = (ivar < nb_std) ?
              db->getArray(iech, cols_std[ivar]) : 1.;
          double prod = value * CT(iclass, ivar + 1);
          total += prod * prod;
        }
        calest.setTstd(iclass, sqrt(total));
      }
    }

    /* Metal Quantity: Estimation */

    if (QT_VARS(QT_EST,ANAM_QT_Q) > 0)
    {
      calest.setQest(nclass-1, anam_discrete_DD->getDDStatZmoy(nclass - 1)
                     * calest.getTest(nclass - 1));
      for (int iclass = nclass - 2; iclass >= 0; iclass--)
        calest.setQest(iclass, calest.getQest(iclass + 1)
                       + anam_discrete_DD->getDDStatZmoy(iclass) *
                       (calest.getTest(iclass) - calest.getTest(iclass + 1)));
    }

    /* Metal Quantity: Standard Deviation */

    if (QT_VARS(QT_STD,ANAM_QT_Q) > 0)
    {
      for (int iclass = 0; iclass < nclass; iclass++)
      {
        double total = 0.;
        for (int ivar = 0; ivar < nclass - 1; ivar++)
        {
          double value = (ivar < nb_std) ?
              db->getArray(iech, cols_std[ivar]) : 1.;
          double prod = value * CQ(iclass, ivar + 1);
          total += prod * prod;
        }
        calest.setQstd(iclass, sqrt(total));
      }
    }

    /* Z: Estimation */

    double zestim = 0.;
    if (QT_VARS(QT_EST,ANAM_QT_Z) > 0)
    {
      zestim = anam_discrete_DD->getDDStatZmoy(
          nclass - 1) * calest.getTest(nclass-1);
      for (int iclass = 0; iclass < nclass - 1; iclass++)
        zestim += anam_discrete_DD->getDDStatZmoy(iclass)
            * (calest.getTest(iclass) - calest.getTest(iclass + 1));
    }

    /* Z: Standard Deviation */

    double zstdev = 0.;
    if (QT_VARS(QT_STD,ANAM_QT_Z) > 0)
    {
      double total = 0;
      for (int ivar = 0; ivar < nclass - 1; ivar++)
      {
        double value = (ivar < nb_std) ?
            db->getArray(iech, cols_std[ivar]) : 1.;
        double prod = value * anam_discrete_DD->getDDStatCnorm(ivar);
        total += prod * prod;
      }
      zstdev = sqrt(total);
    }

    /* Store the results */

    if (ncutmine > 0)
    {
      st_interpolate_qt_local(z_max, cutmine, calest, calcut);
      calcut.calculateBenefitGrade();
      st_recovery_local(db, iech, iptr, ncode, codes, qt_vars, zestim, zstdev,
                        calcut);
    }
    else
    {
      calest.calculateBenefitGrade();
      st_recovery_local(db, iech, iptr, ncode, codes, qt_vars, zestim, zstdev,
                        calest);
    }
  }
  return 0;
}

/*****************************************************************************/
/*!
 **  Calculate Experimental Grade-Tonnage curves from factors
 **  Case of Discrete Indicator Residuals
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  anam         Point anamorphosis
 ** \param[in]  ncutmine     Number of required cutoffs
 ** \param[in]  cutmine      Array of the requested cutoffs
 ** \param[in]  z_max        Maximum grade array (only for QT interpolation)
 ** \param[in]  flag_correct 1 if Tonnage order relationship must be corrected
 ** \param[in]  nb_est       Number of columns for factor estimation
 ** \param[in]  cols_est     Array of columns for factor estimation
 ** \param[in]  nb_std       Number of columns for factor st. dev.
 ** \param[in]  cols_std     Array of columns for factor st. dev.
 ** \param[in]  iptr         Rank for storing the results
 ** \param[in]  ncode        Number of stored results
 ** \param[in]  codes        Array of codes for stored results
 ** \param[in]  qt_vars      Array of variables to be calculated
 **
 ** \param[out] calest       Selectivity
 ** \param[out] calcut       Translated selectivity
 *****************************************************************************/
static int st_anam_factor2qt_discrete_IR(Db *db,
                                         AAnam *anam,
                                         int ncutmine,
                                         double *cutmine,
                                         double z_max,
                                         int flag_correct,
                                         int nb_est,
                                         int *cols_est,
                                         int nb_std,
                                         int *cols_std,
                                         int iptr,
                                         int ncode,
                                         int *codes,
                                         int *qt_vars,
                                         Selectivity& calest,
                                         Selectivity& calcut)
{
  AnamDiscreteIR *anam_discrete_IR = dynamic_cast<AnamDiscreteIR*>(anam);
  int nech = db->getSampleNumber();
  int ncleff = MAX(nb_est, nb_std);

  /* Calculate the Recovery Functions from the factors */

  for (int iech = 0; iech < nech; iech++)
  {
    if (st_is_sample_skipped(db, iech, nb_est, cols_est, nb_std, cols_std))
      continue;

    /* Calculate the tonnage and the recovered grade */

    double total = 0.;
    for (int ivar = 0; ivar < ncleff; ivar++)
    {
      double value = db->getArray(iech, cols_est[ivar]);
      total += value;
      calest.setTest(ivar, total * anam_discrete_IR->getIRStatT(ivar + 1));
    }

    /* Correct order relationship */

    if (flag_correct) calest.correctTonnageOrder();

    /* Tonnage: Standard Deviation */

    if (QT_VARS(QT_STD,ANAM_QT_T) > 0)
    {
      total = 0.;
      for (int ivar = 0; ivar < ncleff; ivar++)
      {
        double value = db->getArray(iech, cols_std[ivar]);
        total += value * value;
        calest.setTstd(ivar, sqrt(total) * anam_discrete_IR->getIRStatT(ivar + 1));
      }
    }

    /* Metal Quantity: Estimation */

    if (QT_VARS(QT_EST,ANAM_QT_Q) > 0)
    {
      calest.setQest(ncleff-1, anam_discrete_IR->getIRStatZ(ncleff)
                     * calest.getTest(ncleff-1));
      for (int ivar = ncleff - 2; ivar >= 0; ivar--)
      {
        calest.setQest(ivar,
                       calest.getQest(ivar+1) + anam_discrete_IR->getIRStatZ(ivar + 1)
                       * (calest.getTest(ivar) - calest.getTest(ivar + 1)));
      }
    }

    /* Metal Quantity: Standard Deviation */

    if (QT_VARS(QT_STD,ANAM_QT_Q) > 0)
    {
      for (int ivar = 0; ivar < ncleff; ivar++)
      {
        total = 0.;
        for (int jvar = 0; jvar < ivar; jvar++)
        {
          double value = db->getArray(iech, cols_std[jvar]);
          total += value * value;
        }
        double prod = anam_discrete_IR->getIRStatB(ivar + 1)
            + anam_discrete_IR->getIRStatZ(ivar + 1) * anam_discrete_IR->getIRStatT(
                ivar + 1);
        total *= prod * prod;
        for (int jvar = ivar + 1; jvar < ncleff; jvar++)
        {
          double prod = db->getArray(iech, cols_std[jvar])
              * anam_discrete_IR->getIRStatB(ivar + 1);
          total += prod * prod;
        }
        calest.setQstd(ivar, sqrt(total));
      }
    }

    /* Z: Estimation */

    double zestim = 0.;
    if (QT_VARS(QT_EST,ANAM_QT_Z) > 0)
    {
      zestim = anam_discrete_IR->getIRStatZ(ncleff) * calest.getTest(ncleff-1);
      for (int ivar = 0; ivar < ncleff - 1; ivar++)
        zestim += anam_discrete_IR->getIRStatZ(ivar + 1)
            * (calest.getTest(ivar) - calest.getTest(ivar + 1));
    }

    /* Z: Standard Deviation */

    double zstdev = 0.;
    if (QT_VARS(QT_STD,ANAM_QT_Z) > 0)
    {
      total = 0.;
      for (int ivar = 0; ivar < ncleff; ivar++)
      {
        double prod = db->getArray(iech, cols_std[ivar])
            * anam_discrete_IR->getIRStatB(ivar + 1);
        total += prod * prod;
      }
      zstdev = sqrt(total);
    }

    /* Store the results */

    if (ncutmine > 0)
    {
      st_interpolate_qt_local(z_max, cutmine, calest, calcut);
      calcut.calculateBenefitGrade();
      st_recovery_local(db, iech, iptr, ncode, codes, qt_vars, zestim, zstdev,
                        calcut);
    }
    else
    {
      calest.calculateBenefitGrade();
      st_recovery_local(db, iech, iptr, ncode, codes, qt_vars, zestim, zstdev,
                        calest);
    }
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Calculate the recoveries (z,T,Q,m,B) starting from the factors
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  anam         Point anamorphosis
 ** \param[in]  ncutmine     Number of required cutoffs
 ** \param[in]  cutmine      Array of the requested cutoffs
 ** \param[in]  z_max        Maximum grade array (only for QT interpolation)
 ** \param[in]  flag_correct 1 if Tonnage order relationship must be corrected
 ** \param[in]  ncode        Number of operators
 ** \param[in]  codes        Array of codes for stored results
 ** \param[in]  nb_est       Number of columns for factor estimation
 ** \param[in]  cols_est     Array of columns for factor estimation
 ** \param[in]  nb_std       Number of columns for factor st. dev.
 ** \param[in]  cols_std     Array of columns for factor st. dev.
 **
 ** \param[out] ncut         Actual number of cutoffs
 ** \param[out] qt_vars      Array for storage (Dimension: 2*ANAM_N_QT)
 **
 ** \remark If the argument 'zcut' is provided, the recovery curves are
 ** \remark calculated for these cutoffs. Otherwise, they are calculated
 ** \remark for the estimated cutoffs in the discrete case
 **
 *****************************************************************************/
int anam_factor2qt(Db *db,
                   AAnam *anam,
                   int ncutmine,
                   double *cutmine,
                   double z_max,
                   int flag_correct,
                   int ncode,
                   int *codes,
                   int nb_est,
                   int *cols_est,
                   int nb_std,
                   int *cols_std,
                   int *ncut,
                   int *qt_vars)
{
  int iptr = -1;
  int flag_inter = 0;
  int verbose = 0;
  Selectivity calest;
  Selectivity calcut;

  AnamHermite *anam_hermite = dynamic_cast<AnamHermite*>(anam);
  AnamDiscreteDD *anam_discrete_DD = dynamic_cast<AnamDiscreteDD*>(anam);
  AnamDiscreteIR *anam_discrete_IR = dynamic_cast<AnamDiscreteIR*>(anam);

  /* Preliminary checks */

  (*ncut) = 0;
  if (db == nullptr) return 1;
  if (anam == nullptr) return 1;
  if (nb_est <= 0 && nb_std <= 0)
  {
    messerr("The number of factors is zero");
    return 1;
  }
  int nvar = MAX(nb_est, nb_std);

  /* Get the number of initial cutoffs */

  int nmax = 0;
  switch (anam->getType().toEnum())
  {
    case EAnam::E_HERMITIAN:
      nmax = anam_hermite->getNbPoly();
      if (nvar >= nmax)
      {
        messerr(
            "Number of factors (%d) must be smaller than Number of classes (%d)",
            nvar, nmax);
        return 1;
      }
      break;

    case EAnam::E_DISCRETE_DD:
      nmax = anam_discrete_DD->getNClass();
      if (nvar >= nmax)
      {
        messerr(
            "Number of factors (%d) must be smaller than Number of classes (%d)",
            nvar, nmax);
        return 1;
      }
      if (ncutmine <= 0)
        ncutmine = nmax;
      else
        flag_inter = 1;
      break;

    case EAnam::E_DISCRETE_IR:
      nmax = anam_discrete_IR->getNCut();
      if (nvar > nmax)
      {
        messerr(
            "Number of factors (%d) must be smaller or equal to Number of cutoffs",
            nvar, nmax);
        return 1;
      }
      if (ncutmine <= 0)
        ncutmine = nmax;
      else
        flag_inter = 1;
      break;

    default:
      messerr("This Anamorphosis cannot be used in this method");
      return 1;
  }

  /* Analyzing the code requirements */

  int nvarout = st_code_analyze(verbose, ncode, codes,
                                nb_est, nb_std, ncutmine,
                                TEST, flag_inter, qt_vars);
  if (nvarout <= 0) return 1;

  /* Variable allocation */

  (*ncut) = ncutmine;
  iptr = db->addColumnsByConstant(nvarout, TEST);
  if (iptr < 0) return 1;

  /* Core allocation */

  calest = Selectivity(nmax);
  if (ncutmine > 0)
    calcut = Selectivity(ncutmine);

  /* Dispatch according to the type of Anamorphosis */

  switch (anam->getType().toEnum())
  {
    case EAnam::E_HERMITIAN:
      if (st_anam_factor2qt_hermitian(db, anam, ncutmine, cutmine, nb_est,
                                      cols_est, nb_std, cols_std, iptr, ncode,
                                      codes, qt_vars, calest)) return 1;
      break;

    case EAnam::E_DISCRETE_DD:
      if (st_anam_factor2qt_discrete_DD(db, anam, ncutmine, cutmine, z_max,
                                        flag_correct, nb_est, cols_est, nb_std,
                                        cols_std, iptr, ncode, codes, qt_vars,
                                        calest, calcut)) return 1;
      break;

    case EAnam::E_DISCRETE_IR:
      if (st_anam_factor2qt_discrete_IR(db, anam, ncutmine, cutmine, z_max,
                                        flag_correct, nb_est, cols_est, nb_std,
                                        cols_std, iptr, ncode, codes, qt_vars,
                                        calest, calcut)) return 1;
      break;

    default:
      messerr("This method is not programmed yet for this anamorphosis");
      return 1;
  }

  return 0;
}

/****************************************************************************/
/*!
 **  Interpolate the Grade-Tonnage curves
 **
 ** \param[in] verbose  Verbose flag
 ** \param[in] zcutmine Array of cutoffs
 ** \param[in] calest   Selectivity
 **
 ** \param[out] calcut  Selectivity
 **
 *****************************************************************************/
void selectivity_interpolate(int verbose,
                             double *zcutmine,
                             Selectivity& calest,
                             Selectivity& calcut)
{
  st_interpolate_qt_global(zcutmine, calest, calcut);
  calcut.calculateBenefitGrade();
  if (verbose) calcut.dumpGini();
}

/*****************************************************************************/
/*!
 **  Calculate the Uniform Conditioning
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  anam         Point anamorphosis
 ** \param[in]  att_est      Rank of the Kriging estimate
 ** \param[in]  att_var      Rank of the Variance of Kriging estimate
 ** \param[in]  ncutmine     Number of required cutoffs
 ** \param[in]  cutmine      Array of the requested cutoffs
 ** \param[in]  proba        Probability
 ** \param[in]  var_bloc     Change of support coefficient
 ** \param[in]  ncode        Number of operators
 ** \param[in]  codes        Array of codes for stored results
 ** \param[in]  verbose      Verbose option
 **
 ** \param[out]  qt_vars     Array of results
 **
 *****************************************************************************/
int uc_f(Db *db,
         AAnam *anam,
         int att_est,
         int att_var,
         int ncutmine,
         double *cutmine,
         double proba,
         double var_bloc,
         int ncode,
         int *codes,
         int verbose,
         int *qt_vars)
{
  int error, nbpoly, iptr, iptr_sV, iptr_yV, nvarout;
  double yc, sv, yv, ore, metal, mean, variance, varb, r_coef, zvstar, varv;
  double vv_min, vv_max, sv_min, sv_max, zv_min, zv_max, yv_min, yv_max;
  VectorDouble psi_hn, hn;
  VectorVectorDouble phi_b_zc;
  Selectivity calest;

  /* Initializations */

  error = 1;
  iptr = iptr_sV = iptr_yV = -1;
  vv_min = zv_min = sv_min = yv_min =  1.e30;
  vv_max = zv_max = sv_max = yv_max = -1.e30;
  AnamHermite *anam_hermite = dynamic_cast<AnamHermite*>(anam);
  anam_hermite->setFlagBound(1);

  /* Preliminary checks */

  if (db == nullptr) goto label_end;
  if (anam == nullptr) goto label_end;
  if (anam->getType() != EAnam::HERMITIAN)
  {
    messerr("Uniform Conditioning is restricted to Hermitian Anamorphosis");
    goto label_end;
  }

  if (anam_hermite->getVariance() <= var_bloc)
  {
    messerr("The coefficient of change of support (%lf)", var_bloc);
    messerr("must be smaller than the variance (%lf)",
            anam_hermite->getVariance());
    goto label_end;
  }
  if (ncutmine <= 0)
  {
    messerr("You must define some cutoff values");
    goto label_end;
  }
  nbpoly = anam_hermite->getNbPoly();
  mean = anam_hermite->getMean();
  variance = anam_hermite->getVariance();

  /* Core allocation */

  psi_hn.resize(nbpoly);
  phi_b_zc.resize(ncutmine);

  /* Memorize the punctual Hermite polynomials */

  for (int ip = 0; ip < nbpoly; ip++)
    psi_hn[ip] = anam_hermite->getPsiHn(ip);

  /* Add variables for storage */

  iptr_sV = db->addColumnsByConstant(1, TEST);
  if (iptr_sV < 0) goto label_end;
  iptr_yV = db->addColumnsByConstant(1, TEST);
  if (iptr_yV < 0) goto label_end;

  /* Analyzing the codes */

  nvarout = st_code_analyze(verbose, ncode, codes, 1, 1, ncutmine, proba, 0,
                            qt_vars);
  if (nvarout <= 0) goto label_end;
  if (QT_FLAG(ANAM_QT_Z))
  {
    messerr("The recovery option 'Z' is not available in this function");
    goto label_end;
  }
  iptr = db->addColumnsByConstant(nvarout, TEST);
  if (iptr < 0) goto label_end;

  /* Core allocation */

  calest = Selectivity(ncutmine);

  /* Transforming Point anamorphosis into Block Anamorphosis */

  if (anam_point_to_block(anam, 0, var_bloc, TEST, TEST)) goto label_end;
  r_coef = anam_hermite->getRCoef();
  varb = anam_hermite->getVariance();

  /* Transform cutmine into gaussian equivalent */

  for (int icut = 0; icut < ncutmine; icut++)
    cutmine[icut] = anam_hermite->RawToTransformValue(cutmine[icut]);

  /* Fill the array phi_b_zc */

  for (int icut = 0; icut < ncutmine; icut++)
  {
    hn = hermitePolynomials(cutmine[icut], 1., nbpoly);
    phi_b_zc[icut] = hermiteCoefMetal(cutmine[icut], hn);
  }

  /* Computing S and Y on panels */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    anam_hermite->setPsiHn(psi_hn);
    anam_hermite->calculateMeanAndVariance();
    zvstar = db->getArray(iech, att_est);
    varv = db->getArray(iech, att_var);
    if (anam_point_to_block(anam, 0, varv, TEST, TEST)) goto label_end;
    db->setArray(iech, iptr_sV, anam_hermite->getRCoef());
    db->setArray(iech, iptr_yV, anam_hermite->RawToTransformValue(zvstar));

    if (varv < vv_min) vv_min = varv;
    if (varv > vv_max) vv_max = varv;
    if (zvstar < zv_min) zv_min = zvstar;
    if (zvstar > zv_max) zv_max = zvstar;
  }

  /* Loop on the panels to compute the grade-tonnage functions */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    sv = db->getArray(iech, iptr_sV);
    yv = db->getArray(iech, iptr_yV);

    /* Loop on the cutoffs */

    for (int icut = 0; icut < ncutmine; icut++)
    {
      yc = cutmine[icut];
      ore = 1.
          - law_cdf_gaussian(
              (yc - yv * sv / r_coef) / sqrt(
                  1. - (sv / r_coef) * (sv / r_coef)));
      hn = hermitePolynomials(yv, 1., nbpoly);
      for (int ih = 0; ih < nbpoly; ih++)
        hn[ih] *= pow(sv / r_coef, (double) ih);
      matrix_product(1, nbpoly, 1, hn.data(), phi_b_zc[icut].data(), &metal);

      if (sv < sv_min) sv_min = sv;
      if (sv > sv_max) sv_max = sv;
      if (yv < yv_min) yv_min = yv;
      if (yv > yv_max) yv_max = yv;

      /* Storing the grade-tonnage functions */

      calest.setZcut(icut, yc);
      calest.setTest(icut, ore);
      calest.setQest(icut, metal);
    }

    calest.calculateBenefitGrade();
    st_recovery_local(db, iech, iptr, ncode, codes, qt_vars, TEST, TEST,
                      calest);
  }

  /* Verbose printout (optional) */

  if (verbose)
  {
    message("Uniform Conditioning on %d panels and %d cutoffs\n",
            db->getSampleNumber(true), ncutmine);
    message("- Number of Polynomials = %d\n", nbpoly);
    message("- Mean                  = %lf\n", mean);
    message("- Punctual Variance     = %lf\n", variance);
    message("- Block Variance        = %lf\n", varb);
    message("- Change of Support     = %lf\n", r_coef);
    message("- var_V in [%lf, %lf]\n", vv_min, vv_max);
    message("- S_V   in [%lf, %lf]\n", sv_min, sv_max);
    message("- Z_V   in [%lf, %lf]\n", zv_min, zv_max);
    message("- Y_V   in [%lf, %lf]\n", yv_min, yv_max);
  }

  /* Set the error return code */

  error = 0;

  label_end:
  if (iptr_sV >= 0) db->deleteColumnByUID(iptr_sV);
  if (iptr_yV >= 0) db->deleteColumnByUID(iptr_yV);
  return (error);
}

/*****************************************************************************/
/*!
 **  Correct the estimation and st. deviation of estimation if KO
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  iech         Rank of the sample
 ** \param[in]  att_est      Rank of the Kriging estimate
 ** \param[in]  att_std      Rank of the St, Deviation of Kriging estimate
 **
 ** \param[out] krigest      Kriging estimation
 ** \param[out] krigstd      Standard deviation of the estimation error
 **
 ** \remarks In SK: krigstd returns the standard deviation of the estimation error
 ** \remarks In OK: krigstd reads the square root of the estimation variance
 ** \remarks and returns the standard deviation
 **
 *****************************************************************************/
static void st_correct_from_OK(Db *db,
                               int iech,
                               int att_est,
                               int att_std,
                               int /*flag_OK*/,
                               double *krigest,
                               double *krigstd)
{
  *krigest = db->getArray(iech, att_est);
  *krigstd = db->getArray(iech, att_std);

//  if (flag_OK)
//  {
//    double var2 = 1. - (*krigstd) * (*krigstd);
//    double stdv = sqrt(var2);
//    *krigstd  = stdv;
//  }
}

/*****************************************************************************/
/*!
 **  Starting from the initial pointer, return the pointer for the estimation
 **  as well as the one for the standard deviation (if required)
 **
 ** \param[in]  iptr_init    Initial pointer
 ** \param[in]  ncutmine     Number of cutoffs
 ** \param[in]  icut         Rank of the cutoff
 ** \param[in]  flag_est     1 for computing the Estimation
 ** \param[in]  flag_std     1 for computing the St. Deviation
 **
 ** \param[out] iptr_est     Starting pointer for the Estimation
 ** \param[out] iptr_std     Starting pointer for the St. Deviation
 **
 ** \remarks If not used the output pointers are set to -1
 **
 *****************************************************************************/
static void st_get_starting_pointers(int iptr_init,
                                     int ncutmine,
                                     int icut,
                                     int flag_est,
                                     int flag_std,
                                     int *iptr_est,
                                     int *iptr_std)
{
  int iptr, ncut;

  ncut = MAX(1, ncutmine);
  iptr = iptr_init + icut;

  if (flag_est)
  {
    *iptr_est = iptr;
    iptr += ncut;
  }
  else
  {
    *iptr_est = -1;
  }
  if (flag_std)
  {
    *iptr_std = iptr;
    iptr += ncut;
  }
  else
  {
    *iptr_std = -1;
  }
}

/*****************************************************************************/
/*!
 **  Prepare the vectors of estimation and st. dev.
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  att_est      Rank of the Kriging estimate
 ** \param[in]  att_std      Rank of the St, Deviation of Kriging estimate
 ** \param[in]  flag_OK      1 if kriging is performed with OK
 **
 ** \param[out] krigest      Vector of estamations
 ** \param[out] krigstd      Vector of standard deviation
 **
 *****************************************************************************/
static void st_ce_get_vectors(Db *db,
                              int att_est,
                              int att_std,
                              int flag_OK,
                              VectorDouble &krigest,
                              VectorDouble &krigstd)
{
  int nech = db->getSampleNumber();

  krigest.resize(nech);
  krigstd.resize(nech);
  for (int iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    st_correct_from_OK(db, iech, att_est, att_std, flag_OK, &krigest[iech],
                       &krigstd[iech]);
  }
}

/*****************************************************************************/
/*!
 **  Calculate the Conditional value and variance in the Gaussian Model
 **
 ** \return Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  nbsimu       Number of Monte Carlo simulations (0: Hermite)
 ** \param[in]  phis         Array of the Polynomial expansion
 ** \param[in]  att_est      Rank of the Kriging estimate
 ** \param[in]  att_std      Rank of the St, Deviation of Kriging estimate
 ** \param[in]  flag_OK      1 if kriging is performed with OK
 ** \param[in]  flag_est     Flag for calculation of the Estimation
 ** \param[in]  flag_std     Flag for calculation of the St. Deviation
 ** \param[in]  iptr_Z       Address of the Z variable
 **
 *****************************************************************************/
static int st_ce_compute_Z(Db *db,
                           int nbsimu,
                           const VectorDouble phis,
                           int att_est,
                           int att_std,
                           int flag_OK,
                           int flag_est,
                           int flag_std,
                           int iptr_Z)

{
  VectorDouble krigest, krigstd, valest, valstd;
  int iptr_est, iptr_std;

  st_get_starting_pointers(iptr_Z, 1, 0, flag_est, flag_std, &iptr_est,
                           &iptr_std);

  st_ce_get_vectors(db, att_est, att_std, flag_OK, krigest, krigstd);

  if (nbsimu <= 0)
  {
    valest = hermiteCondExp(krigest, krigstd, phis);
    valstd = hermiteCondStd(krigest, krigstd, phis);
  }
  else
  {
    valest = MCCondExp(krigest, krigstd, phis, nbsimu);
    valstd = MCCondStd(krigest, krigstd, phis, nbsimu);
  }

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    if (iptr_est >= 0) db->setArray(iech, iptr_est, valest[iech]);
    if (iptr_std >= 0) db->setArray(iech, iptr_std, valstd[iech]);
  }

  return (0);
}

/*****************************************************************************/
/*!
 **  Calculate the Conditional value and variance in the Gaussian Model
 **
 ** \return Error return code
 **
 ** \param[in]  krigest      Estimation value
 ** \param[in]  krigstd      Standard deviation of Estimation value
 ** \param[in]  phis         Array of the Polynomial expansion
 **
 *****************************************************************************/
double ce_compute_Z2(double krigest, double krigstd, const VectorDouble &phis)
{
  VectorDouble dd;

  /* Core allocation */

  int nbpoly = static_cast<int>(phis.size());
  dd.resize(nbpoly * nbpoly);

  /* Loading the Conditional Expectation arrays */

  message("calculating DD with nbpoly = %d\n", nbpoly);
  for (int ih = 0; ih < nbpoly; ih++)
    for (int jh = 0; jh < nbpoly; jh++)
    {
      DD(ih,jh) = (ih + jh >= nbpoly) ? 0 :
                                (pow(-1., ih) * phis[ih + jh]
                                 * sqrt(ut_cnp(ih + jh, jh)));
    }

  /* Loop on the samples */

  double krigvar = krigstd * krigstd;
  double krigrho = sqrt(1. - krigvar);
  VectorDouble hh = hermitePolynomials(krigest / krigrho, 1., nbpoly);

  /* Calculating conditional variance */

  double valstd = 0.;
  for (int ih = 1; ih < nbpoly; ih++)
  {
    double factor = 0.;
    for (int jh = 0; jh < nbpoly; jh++)
      factor += hh[jh] * pow(krigrho, (double) jh) * DD(ih, jh);
    valstd += pow(krigvar, (double) ih) * factor * factor;
  }
  valstd = sqrt(valstd);

  return valstd;
}

/*****************************************************************************/
/*!
 **  Calculate the Tonnage by Conditional Expectation
 **
 ** \return Error return code
 **
 ** \param[in]  mode         1 for T (Proba. above); 2 for Proba. below
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  ncutmine     Number of required cutoffs
 ** \param[in]  nbsimu       Number of Monte Carlo simulations (0: Hermite)
 ** \param[in]  yc           Array of the requested cutoffs (gaussian scale)
 ** \param[in]  att_est      Rank of the Kriging estimate
 ** \param[in]  att_std      Rank of the St, Deviation of Kriging estimate
 ** \param[in]  flag_OK      1 if kriging is performed with OK
 ** \param[in]  flag_est     1 for computing the Estimation
 ** \param[in]  flag_std     1 for computing the St. Deviation
 ** \param[in]  iptr_T       Address of the T variable
 **
 *****************************************************************************/
static int st_ce_compute_T(int mode,
                           Db *db,
                           int ncutmine,
                           int nbsimu,
                           double *yc,
                           int att_est,
                           int att_std,
                           int flag_OK,
                           int flag_est,
                           int flag_std,
                           int iptr_T)
{
  VectorDouble krigest, krigstd, valest, valstd;
  int iptr_est, iptr_std;

  /* Loop on the samples */

  iptr_est = iptr_std = -1;
  for (int icut = 0; icut < ncutmine; icut++)
  {
    st_get_starting_pointers(iptr_T, ncutmine, icut, flag_est, flag_std,
                             &iptr_est, &iptr_std);

    st_ce_get_vectors(db, att_est, att_std, flag_OK, krigest, krigstd);

    if (nbsimu <= 0)
    {
      valest = hermiteIndicator(yc[icut], krigest, krigstd);
      valstd = hermiteIndicatorStd(yc[icut], krigest, krigstd);
    }
    else
    {
      valest = MCIndicator(yc[icut], krigest, krigstd, nbsimu);
      valstd = MCIndicatorStd(yc[icut], krigest, krigstd, nbsimu);
    }

    for (int iech = 0; iech < db->getSampleNumber(); iech++)
    {
      if (!db->isActive(iech)) continue;
      if (iptr_est >= 0)
      {
        if (mode == 1)
          db->setArray(iech, iptr_est, valest[iech]);
        else
          db->setArray(iech, iptr_est, 1. - valest[iech]);
      }
      if (iptr_std >= 0) db->setArray(iech, iptr_std, valstd[iech]);
    }
  }

  return (0);
}

/*****************************************************************************/
/*!
 **  Calculate the Quantile by Conditional Expectation
 **
 ** \return Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  anam         AAnam structure
 ** \param[in]  proba        Probability threshold
 ** \param[in]  att_est      Rank of the Kriging estimate
 ** \param[in]  att_std      Rank of the St, Deviation of Kriging estimate
 ** \param[in]  flag_OK      1 if kriging is performed with OK
 ** \param[in]  iptr_QUANT   Address of the Quantile
 **
 *****************************************************************************/
static int st_ce_compute_quant(Db *db,
                               AAnam *anam,
                               double proba,
                               int att_est,
                               int att_std,
                               int flag_OK,
                               int iptr_QUANT)
{
  double krigest, krigstd;

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    st_correct_from_OK(db, iech, att_est, att_std, flag_OK, &krigest, &krigstd);
    double y = krigest + krigstd * law_invcdf_gaussian(proba);
    db->setArray(iech, iptr_QUANT, anam->TransformToRawValue(y));
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Calculate the Metal Quantity by Conditional Expectation
 **
 ** \return Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  ncutmine     Number of required cutoffs
 ** \param[in]  nbsimu       Number of Monte Carlo simulations (0: Hermite)
 ** \param[in]  yc           Array of the requested cutoffs (gaussian scale)
 ** \param[in]  phis         Array of the Polynomial expansion
 ** \param[in]  att_est      Rank of the Kriging estimate
 ** \param[in]  att_std      Rank of the St, Deviation of Kriging estimate
 ** \param[in]  flag_OK      1 if kriging is performed with OK
 ** \param[in]  flag_est     1 for computing the Estimation
 ** \param[in]  flag_std     1 for computing the St. Deviation
 ** \param[in]  iptr_Q       Address of the Q variable
 **
 *****************************************************************************/
static int st_ce_compute_Q(Db *db,
                           int ncutmine,
                           int nbsimu,
                           double *yc,
                           VectorDouble phis,
                           int att_est,
                           int att_std,
                           int flag_OK,
                           int flag_est,
                           int flag_std,
                           int iptr_Q)
{
  VectorDouble krigest, krigstd, valest, valstd;
  int iptr_est, iptr_std;

  /* Loop on the cutoffs */

  iptr_est = iptr_std = -1;
  for (int icut = 0; icut < ncutmine; icut++)
  {
    st_get_starting_pointers(iptr_Q, ncutmine, icut, flag_est, flag_std,
                             &iptr_est, &iptr_std);

    st_ce_get_vectors(db, att_est, att_std, flag_OK, krigest, krigstd);

    if (nbsimu <= 0)
    {
      valest = hermiteMetal(yc[icut], krigest, krigstd, phis);
      valstd = hermiteMetalStd(yc[icut], krigest, krigstd, phis);
    }
    else
    {
      valest = MCMetal(yc[icut], krigest, krigstd, phis, nbsimu);
      valstd = MCMetalStd(yc[icut], krigest, krigstd, phis, nbsimu);
    }

    for (int iech = 0; iech < db->getSampleNumber(); iech++)
    {
      if (!db->isActive(iech)) continue;
      if (iptr_est >= 0) db->setArray(iech, iptr_est, valest[iech]);
      if (iptr_std >= 0) db->setArray(iech, iptr_std, valstd[iech]);
    }
  }

  return (0);
}

/*****************************************************************************/
/*!
 **  Calculate the Conventional Benefit by Conditional Expectation
 **
 ** \return Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  ncutmine     Number of required cutoffs
 ** \param[in]  cutmine      Array of the requested cutoffs
 ** \param[in]  count        Number of items (Estim + St. Dev.)
 ** \param[in]  iptr_T       Address of the Tonnage
 ** \param[in]  iptr_Q       Address of the Metal Quantity
 ** \param[in]  iptr_B       Address of the output quantity
 **
 *****************************************************************************/
static int st_ce_compute_B(Db *db,
                           int ncutmine,
                           double *cutmine,
                           int count,
                           int iptr_T,
                           int iptr_Q,
                           int iptr_B)
{
  double t, q, b;
  int jptr_T, jptr_Q, jptr_B;

  /* Loop on the samples */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;

    /* Loop on the cutoffs */

    jptr_T = iptr_T;
    jptr_Q = iptr_Q;
    jptr_B = iptr_B;
    for (int icut = 0; icut < ncutmine; icut++)
    {
      t = db->getArray(iech, jptr_T);
      q = db->getArray(iech, jptr_Q);
      b = q - t * cutmine[icut];
      db->setArray(iech, jptr_B, b);
      jptr_T += count;
      jptr_Q += count;
      jptr_B++;
    }
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Calculate the Average Grade by Conditional Expectation
 **
 ** \return Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  ncutmine     Number of required cutoffs
 ** \param[in]  count        Number of items (Estim + St. Dev.)
 ** \param[in]  iptr_T       Address of the Tonnage
 ** \param[in]  iptr_Q       Address of the Metal Quantity
 ** \param[in]  iptr_M       Address of the output quantity
 **
 *****************************************************************************/
static int st_ce_compute_M(Db *db,
                           int ncutmine,
                           int count,
                           int iptr_T,
                           int iptr_Q,
                           int iptr_M)
{
  double t, q, m;
  int jptr_T, jptr_Q, jptr_M;

  /* Loop on the samples */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;

    /* Loop on the cutoffs */

    jptr_T = iptr_T;
    jptr_Q = iptr_Q;
    jptr_M = iptr_M;
    for (int icut = 0; icut < ncutmine; icut++)
    {
      t = db->getArray(iech, jptr_T);
      q = db->getArray(iech, jptr_Q);
      m = (t > EPSILON3) ? q / t :
                           TEST;
      db->setArray(iech, jptr_M, m);
      jptr_T += count;
      jptr_Q += count;
      jptr_M++;
    }
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Convert the set of Z-cutoffs into Y-cutoffs
 **
 ** \return  Pointer to the newly allocated vector of values
 **
 ** \param[in]  anam_hermite Point Hermite anamorphosis
 ** \param[in]  ncutmine     Number of required cutoffs
 ** \param[in]  cutmine      Array of the requested cutoffs
 **
 ** \remarks The returned array must be freed by the calling function
 **
 *****************************************************************************/
static double* st_ztoy_cutoffs(AnamHermite *anam_hermite,
                               int ncutmine,
                               double *cutmine)
{
  double *yc;

  // Initializations

  yc = nullptr;
  if (ncutmine < 0) return (yc);

  // Core allocation

  yc = (double*) mem_alloc(sizeof(double) * ncutmine, 0);
  if (yc == nullptr) return (yc);

  // Loop on the cutoff values

  for (int icut = 0; icut < ncutmine; icut++)
    yc[icut] = anam_hermite->RawToTransformValue(cutmine[icut]);
  return (yc);
}

/*****************************************************************************/
/*!
 **  Calculate the Conditional Expectation
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  anam         Point anamorphosis
 ** \param[in]  att_est      Rank of the Kriging estimate
 ** \param[in]  att_std      Rank of the St, Deviation of Kriging estimate
 ** \param[in]  flag_est     1 for computing the Estimation
 ** \param[in]  flag_std     1 for computing the St. Deviation
 ** \param[in]  flag_OK      1 if kriging has ben performed in Ordinary Kriging
 ** \param[in]  ncutmine     Number of required cutoffs
 ** \param[in]  cutmine      Array of the requested cutoffs
 ** \param[in]  proba        Probability
 ** \param[in]  ncode        Number of operators
 ** \param[in]  codes        Array of codes for stored results
 ** \param[in]  nbsimu       Number of Monte Carlo simulations (0 : Hermite)
 ** \param[in]  verbose      Verbose option
 **
 ** \param[out] qt_vars      Array for storage (Dimension: 2*ANAM_N_QT)
 **
 *****************************************************************************/
int ce_f(Db *db,
         AAnam *anam,
         int att_est,
         int att_std,
         int flag_est,
         int flag_std,
         int flag_OK,
         int ncutmine,
         double *cutmine,
         double proba,
         int ncode,
         int *codes,
         int nbsimu,
         int verbose,
         int *qt_vars)
{
  int error, nbpoly, iptr_Z, iptr_T, iptr_Q, iptr_B, iptr_M, need_T, need_Q,
      count;
  int iptr_est, iptr_std, iptr_PROBA, iptr_QUANT;
  double *yc;

  /* Initializations */

  error = 1;
  count = need_T = need_Q = 0;
  iptr_Z = iptr_T = iptr_Q = iptr_B = iptr_M = iptr_PROBA = iptr_QUANT = -1;
  yc = nullptr;
  AnamHermite *anam_hermite = dynamic_cast<AnamHermite*>(anam);

  if (anam->getType() != EAnam::HERMITIAN)
  {
    messerr("The argument 'anam' must be Gaussian");
    goto label_end;
  }
  if (ncode <= 0)
  {
    messerr("You must specify at least one Recovery Function");
    goto label_end;
  }
  nbpoly = anam_hermite->getNbPoly();

  /* Analyzing the codes */

  count = flag_est + flag_std;
  if (st_code_analyze(verbose, ncode, codes, flag_est, flag_std, ncutmine,
                      proba, 0, qt_vars) <= 0) goto label_end;
  yc = st_ztoy_cutoffs(anam_hermite, ncutmine, cutmine);
  need_T = QT_FLAG(ANAM_QT_T) || QT_FLAG(ANAM_QT_B) || QT_FLAG(ANAM_QT_M) ||
  QT_FLAG(ANAM_QT_PROBA);
  need_Q = QT_FLAG(ANAM_QT_Q) || QT_FLAG(ANAM_QT_B) || QT_FLAG(ANAM_QT_M);
  if (yc == nullptr) need_T = need_Q = 0;

  /* Add the variables */

  if (QT_FLAG(ANAM_QT_Z))
  {
    iptr_Z = db->addColumnsByConstant(count, TEST);
    if (iptr_Z < 0) goto label_end;
  }
  if (need_T)
  {
    iptr_T = db->addColumnsByConstant(count * ncutmine, TEST);
    if (iptr_T < 0) goto label_end;
  }
  if (need_Q)
  {
    iptr_Q = db->addColumnsByConstant(count * ncutmine, TEST);
    if (iptr_Q < 0) goto label_end;
  }
  if (QT_FLAG(ANAM_QT_B) && need_T && need_Q)
  {
    iptr_B = db->addColumnsByConstant(ncutmine, TEST);
    if (iptr_B < 0) goto label_end;
  }
  if (QT_FLAG(ANAM_QT_M) && need_T && need_Q)
  {
    iptr_M = db->addColumnsByConstant(ncutmine, TEST);
    if (iptr_M < 0) goto label_end;
  }
  if (QT_FLAG(ANAM_QT_PROBA) && need_T)
  {
    iptr_PROBA = db->addColumnsByConstant(count * ncutmine, TEST);
    if (iptr_PROBA < 0) goto label_end;
  }
  if (QT_FLAG(ANAM_QT_QUANT))
  {
    iptr_QUANT = db->addColumnsByConstant(1, TEST);
    if (iptr_QUANT < 0) goto label_end;
  }

  /* Optional printout */

  if (verbose)
  {
    if (nbsimu > 0)
      message(
          "Conditional expectation under the gaussian model (Monte-Carlo)\n");
    else
      message("Conditional expectation under the gaussian model (Hermite)\n");
    message(" Max. degree of Hermite polynomials : %6d\n", nbpoly - 1);
    message(" Number of values                   : %6d\n",
            db->getSampleNumber(true));
    message(" Number of cutoffs                  : %6d\n", ncutmine);
    if (nbsimu > 0)
      message(" Number of Monte-Carlo simulations  : %6d\n", nbsimu);
  }

  /* Computing the estimation */

  if (QT_FLAG(ANAM_QT_Z))
  {
    if (st_ce_compute_Z(db, nbsimu, anam_hermite->getPsiHn(), att_est, att_std,
                        flag_OK, flag_est, flag_std, iptr_Z)) goto label_end;
  }

  /* Compute Conditional Expectation for Tonnage */

  if (need_T)
  {
    if (st_ce_compute_T(1, db, ncutmine, nbsimu, yc, att_est, att_std, flag_OK,
                        flag_est, flag_std, iptr_T)) goto label_end;
  }

  /* Compute Conditional Expectation for Metal Quantity */

  if (QT_FLAG(ANAM_QT_Q) && need_Q)
  {
    if (st_ce_compute_Q(db, ncutmine, nbsimu, yc, anam_hermite->getPsiHn(),
                        att_est, att_std, flag_OK, flag_est, flag_std, iptr_Q))
      goto label_end;
  }

  /* Compute Conditional Expectation for Conventional Benefit */

  if (QT_FLAG(ANAM_QT_B) && need_T && need_Q)
  {
    if (st_ce_compute_B(db, ncutmine, cutmine, count, iptr_T, iptr_Q, iptr_B))
      goto label_end;
  }

  /* Compute Conditional Expectation for Average recoveable grade */

  if (QT_FLAG(ANAM_QT_M) && need_T && need_Q)
  {
    if (st_ce_compute_M(db, ncutmine, count, iptr_T, iptr_Q, iptr_M))
      goto label_end;
  }

  /* Compute Conditional Expectation for Tonnage */

  if (QT_FLAG(ANAM_QT_PROBA) && need_T)
  {
    if (st_ce_compute_T(2, db, ncutmine, nbsimu, yc, att_est, att_std, flag_OK,
                        flag_est, flag_std, iptr_PROBA)) goto label_end;
  }

  /* Compute Conditional Expectation for Quantile */

  if (QT_FLAG(ANAM_QT_QUANT))
  {
    st_get_starting_pointers(iptr_QUANT, ncutmine, 0, flag_est, flag_std,
                             &iptr_est, &iptr_std);
    if (st_ce_compute_quant(db, anam, proba, att_est, att_std, flag_OK,
                            iptr_QUANT)) goto label_end;
  }

  /* Set the error return code */

  error = 0;

  label_end: yc = (double*) mem_free((char* ) yc);
  if (!QT_FLAG(ANAM_QT_T))
    (void) db_attribute_del_mult(db, iptr_T, count * ncutmine);
  if (!QT_FLAG(ANAM_QT_Q))
    (void) db_attribute_del_mult(db, iptr_Q, count * ncutmine);
  return (error);
}

/*****************************************************************************/
/*!
 **  Calculate the factors corresponding to an input data vector
 **
 ** \return  Error return code
 **
 ** \param[in]  anam        anamorphosis model
 ** \param[in]  db          Db structure
 ** \param[in]  ifacs       Array of factor ranks
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int anamZToFactor(AAnam *anam,
                  Db *db,
                  const VectorInt& ifacs,
                  const NamingConvention& namconv)
{
  if (anam == nullptr)
  {
    messerr("You must define the 'anam' argument");
    return 1;
  }
  if (db == nullptr)
  {
    messerr("You must define the 'db' argument");
    return 1;
  }
  int nvar = db->getVariableNumber();
  if (nvar != 1)
  {
    messerr("This function is only coded for the monovariate Db");
    return 1;
  }
  int nfact = (int) ifacs.size();
  if (nfact <= 0)
  {
    messerr("You must define the list of factors");
    return 1;
  }
  int nmax = anam->getNFactor();
  for (int ifac = 0; ifac < nfact; ifac++)
    if (ifacs[ifac] < 1 || ifacs[ifac] > nmax)
    {
      messerr("Error in the rank of the factor(%d): it should lie in [1,%d]",
              ifacs[ifac], nmax);
      return 1;
    }

  /* Create the factors */

  int iptr = db->addColumnsByConstant(nfact, TEST);
  if (iptr <= 0) return 1;

  // Loop on the samples

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    double zval = db->getVariable(iech, 0);
    if (FFFF(zval)) continue;
    VectorDouble factors = anam->z2factor(zval, ifacs);
    if (factors.empty()) continue;
    for (int ifac = 0; ifac < nfact; ifac++)
      db->setArray(iech, iptr + ifac, factors[ifac]);
  }

  /* Set the error return code */

  namconv.setNamesAndLocators(db, ELoc::Z, nfact, db, iptr, "Factor");

  return 0;
}
