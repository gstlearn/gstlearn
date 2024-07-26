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

#include "Anamorphosis/AnamDiscrete.hpp"
#include "Anamorphosis/AnamDiscreteDD.hpp"
#include "Anamorphosis/AnamDiscreteIR.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Variogram/Vario.hpp"
#include "Polynomials/Hermite.hpp"
#include "Polynomials/MonteCarlo.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Db/Db.hpp"
#include "Stats/Selectivity.hpp"

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
        r_coef = sqrt(anam->invertVariance(cvv));
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
        r_coef = sqrt(anam->invertVariance(cvv));
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
        r_coef = sqrt(anam->invertVariance(cvv));
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
 **  Correct the estimation and st. deviation of estimation if KO
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  iech         Rank of the sample
 ** \param[in]  col_est      Rank of the Kriging estimate
 ** \param[in]  col_std      Rank of the St, Deviation of Kriging estimate
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
                               int col_est,
                               int col_std,
                               bool /*flag_OK*/,
                               double *krigest,
                               double *krigstd)
{
  *krigest = db->getArray(iech, col_est);
  *krigstd = db->getArray(iech, col_std);

//  if (flag_OK)
//  {
//    double var2 = 1. - (*krigstd) * (*krigstd);
//    double stdv = sqrt(var2);
//    *krigstd  = stdv;
//  }
}

/*****************************************************************************/
/*!
 **  Prepare the vectors of estimation and st. dev.
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  col_est      Rank of the Kriging estimate
 ** \param[in]  col_std      Rank of the St, Deviation of Kriging estimate
 ** \param[in]  flag_OK      1 if kriging is performed with OK
 **
 ** \param[out] krigest      Vector of estamations
 ** \param[out] krigstd      Vector of standard deviation
 **
 *****************************************************************************/
static void st_ce_get_vectors(Db *db,
                              int col_est,
                              int col_std,
                              bool flag_OK,
                              VectorDouble &krigest,
                              VectorDouble &krigstd)
{
  int nech = db->getSampleNumber();
  krigest.resize(nech);
  krigstd.resize(nech);
  for (int iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    st_correct_from_OK(db, iech, col_est, col_std, flag_OK, &krigest[iech],
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
 ** \param[in]  anam         Anamorphosis Hermite
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  iptr0        Starting storage address
 ** \param[in]  col_est      Rank of the Kriging estimate
 ** \param[in]  col_std      Rank of the St, Deviation of Kriging estimate
 ** \param[in]  nbsimu       Number of Monte Carlo simulations (0: Hermite)
 ** \param[in]  flag_OK      1 if kriging is performed with OK
 **
 *****************************************************************************/
static int st_ce_Z(Db *db,
                   const AnamHermite *anam,
                   const Selectivity *selectivity,
                   int iptr0,
                   int col_est,
                   int col_std,
                   int nbsimu,
                   bool flag_OK)
{
  VectorDouble krigest, krigstd, valest, valstd;

  st_ce_get_vectors(db, col_est, col_std, flag_OK, krigest, krigstd);

  if (nbsimu <= 0)
  {
    valest = hermiteCondExp(krigest, krigstd, anam->getPsiHns());
    valstd = hermiteCondStd(krigest, krigstd, anam->getPsiHns());
  }
  else
  {
    valest = MCCondExp(krigest, krigstd, anam->getPsiHns(), nbsimu);
    valstd = MCCondStd(krigest, krigstd, anam->getPsiHns(), nbsimu);
  }

  int iptrEst = selectivity->getAddressQTEst(ESelectivity::Z, iptr0);
  int iptrStd = selectivity->getAddressQTStd(ESelectivity::Z, iptr0);
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    if (selectivity->isUsedEst(ESelectivity::Z))
      db->setArray(iech, iptrEst, valest[iech]);
    if (selectivity->isUsedStD(ESelectivity::Z))
      db->setArray(iech, iptrStd, valstd[iech]);
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Calculate the Tonnage by Conditional Expectation
 **
 ** \return Error return code
 **
 ** \param[in]  mode         1 for T (Proba. above); 2 for Proba. below
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  iptr0        Starting storage address
 ** \param[in]  col_est      Rank of the Kriging estimate
 ** \param[in]  col_std      Rank of the St, Deviation of Kriging estimate
 ** \param[in]  ycuts        Vector of Gaussian cutoffs
 ** \param[in]  nbsimu       Number of Monte Carlo simulations (0: Hermite)
 ** \param[in]  flag_OK      1 if kriging is performed with OK
 **
 *****************************************************************************/
static int st_ce_T(int mode,
                   Db *db,
                   const Selectivity *selectivity,
                   int iptr0,
                   int col_est,
                   int col_std,
                   const VectorDouble &ycuts,
                   int nbsimu,
                   bool flag_OK)
{
  VectorDouble krigest, krigstd, valest, valstd;

  /* Loop on the samples */

  for (int icut = 0; icut < selectivity->getNCuts(); icut++)
  {
    st_ce_get_vectors(db, col_est, col_std, flag_OK, krigest, krigstd);

    if (nbsimu <= 0)
    {
      valest = hermiteIndicator(ycuts[icut], krigest, krigstd);
      valstd = hermiteIndicatorStd(ycuts[icut], krigest, krigstd);
    }
    else
    {
      valest = MCIndicator(ycuts[icut], krigest, krigstd, nbsimu);
      valstd = MCIndicatorStd(ycuts[icut], krigest, krigstd, nbsimu);
    }

    int iptrEst = selectivity->getAddressQTEst(ESelectivity::T, iptr0, icut);
    int iptrStd = selectivity->getAddressQTStd(ESelectivity::T, iptr0, icut);
    int iptrProba = selectivity->getAddressQTEst(ESelectivity::PROP, iptr0, icut);
    for (int iech = 0; iech < db->getSampleNumber(); iech++)
    {
      if (!db->isActive(iech)) continue;
      if (mode == 1)
      {
        if (selectivity->isUsedEst(ESelectivity::T))
          db->setArray(iech, iptrEst, valest[iech]);
        if (selectivity->isUsedStD(ESelectivity::T))
          db->setArray(iech, iptrStd, valstd[iech]);
      }
      else
      {
        if (selectivity->isUsedEst(ESelectivity::PROP))
          db->setArray(iech, iptrProba, 1. - valest[iech]);
      }
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
 ** \param[in]  anam         Hermite anamorphosis
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  iptr0        Staring storage address
 ** \param[in]  col_est      Rank of the Kriging estimate
 ** \param[in]  col_std      Rank of the Kriging St, Deviation
 ** \param[in]  proba        Probability threshold
 ** \param[in]  flag_OK      1 if kriging is performed with OK
 **
 *****************************************************************************/
static int st_ce_quant(Db *db,
                       const AnamHermite *anam,
                       const Selectivity *selectivity,
                       int iptr0,
                       int col_est,
                       int col_std,
                       double proba,
                       bool flag_OK)
{
  double krigest, krigstd;

  int iptrQuant = selectivity->getAddressQTEst(ESelectivity::QUANT, iptr0);

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    st_correct_from_OK(db, iech, col_est, col_std, flag_OK, &krigest, &krigstd);
    double y = krigest + krigstd * law_invcdf_gaussian(proba);
    db->setArray(iech, iptrQuant, anam->transformToRawValue(y));
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
 ** \param[in]  anam         Hermite anamorphosis
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  iptr0        Starting storage address
 ** \param[in]  col_est      Rank of the Kriging estimate
 ** \param[in]  col_std      Rank of the St, Deviation of Kriging estimate
 ** \param[in]  ycuts        Array of the requested cutoffs (gaussian scale)
 ** \param[in]  nbsimu       Number of Monte Carlo simulations (0: Hermite)
 ** \param[in]  flag_OK      1 if kriging is performed with OK
 **
 *****************************************************************************/
static int st_ce_Q(Db *db,
                   const AnamHermite *anam,
                   const Selectivity *selectivity,
                   int iptr0,
                   int col_est,
                   int col_std,
                   const VectorDouble &ycuts,
                   int nbsimu,
                   bool flag_OK)
{
  VectorDouble krigest, krigstd, valest, valstd;

  /* Loop on the cutoffs */

  for (int icut = 0; icut < selectivity->getNCuts(); icut++)
  {
    st_ce_get_vectors(db, col_est, col_std, flag_OK, krigest, krigstd);

    if (nbsimu <= 0)
    {
      valest = hermiteMetal(ycuts[icut], krigest, krigstd, anam->getPsiHns());
      valstd = hermiteMetalStd(ycuts[icut], krigest, krigstd, anam->getPsiHns());
    }
    else
    {
      valest = MCMetal(ycuts[icut], krigest, krigstd, anam->getPsiHns(), nbsimu);
      valstd = MCMetalStd(ycuts[icut], krigest, krigstd, anam->getPsiHns(), nbsimu);
    }

    int iptrEst = selectivity->getAddressQTEst(ESelectivity::Q, iptr0, icut);
    int iptrStd = selectivity->getAddressQTStd(ESelectivity::Q, iptr0, icut);
    for (int iech = 0; iech < db->getSampleNumber(); iech++)
    {
      if (!db->isActive(iech)) continue;
      if (selectivity->isUsedEst(ESelectivity::Q))
        db->setArray(iech, iptrEst, valest[iech]);
      if (selectivity->isUsedStD(ESelectivity::Q))
        db->setArray(iech, iptrStd, valstd[iech]);
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
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  iptr0        Starting storage address
 ** \param[in]  ycuts        Array of the requested cutoffs
 **
 *****************************************************************************/
static int st_ce_B(Db *db,
                   const Selectivity *selectivity,
                   int iptr0,
                   const VectorDouble &ycuts)
{
  double t, q;

  /* Loop on the cutoffs */

  for (int icut = 0; icut < selectivity->getNCuts(); icut++)
  {
    int iptrT = selectivity->getAddressQTEst(ESelectivity::T, iptr0, icut);
    int iptrQ = selectivity->getAddressQTEst(ESelectivity::Q, iptr0, icut);
    int iptrB = selectivity->getAddressQTEst(ESelectivity::B, iptr0, icut);

    /* Loop on the samples */

    for (int iech = 0; iech < db->getSampleNumber(); iech++)
    {
      if (!db->isActive(iech)) continue;

      t = db->getArray(iech, iptrT);
      q = db->getArray(iech, iptrQ);
      db->setArray(iech, iptrB, q - t * ycuts[icut]);
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
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  iptr0        Starting storage address
 **
 *****************************************************************************/
static int st_ce_M(Db *db, const Selectivity *selectivity, int iptr0)
{
  double t, q;

  /* Loop on the cutoffs */

  for (int icut = 0; icut < selectivity->getNCuts(); icut++)
  {
    int iptrT = selectivity->getAddressQTEst(ESelectivity::T, iptr0, icut);
    int iptrQ = selectivity->getAddressQTEst(ESelectivity::Q, iptr0, icut);
    int iptrM = selectivity->getAddressQTEst(ESelectivity::M, iptr0, icut);

    /* Loop on the samples */

    for (int iech = 0; iech < db->getSampleNumber(); iech++)
    {
      if (!db->isActive(iech)) continue;

      t = db->getArray(iech, iptrT);
      q = db->getArray(iech, iptrQ);
      db->setArray(iech, iptrM, (t > EPSILON3) ? q / t : TEST);
    }
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Calculate the Conditional Expectation
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  anam         Anamorphosis structure
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  iptr0        Rank of the pointer for storage
 ** \param[in]  col_est      Rank of variable containing Kriging estimate
 ** \param[in]  col_std      Rank of Variable containing Kriging St. deviation
 ** \param[in]  flag_OK      1 if kriging has ben performed in Ordinary Kriging
 ** \param[in]  proba        Probability
 ** \param[in]  nbsimu       Number of Simulation outcomes
 **
 *****************************************************************************/
int _conditionalExpectation(Db *db,
                            AAnam *anam,
                            const Selectivity *selectivity,
                            int iptr0,
                            int col_est,
                            int col_std,
                            bool flag_OK,
                            double proba,
                            int nbsimu)
{
  AnamHermite *anam_hermite = dynamic_cast<AnamHermite*>(anam);

  /* Analyzing the codes */

  VectorDouble ycuts = anam_hermite->rawToTransformVec(selectivity->getZcut());
  int need_T = selectivity->isNeededT();

  /* Computing the estimation */

  if (selectivity->isUsed(ESelectivity::Z))
  {
    if (st_ce_Z(db, anam_hermite, selectivity, iptr0, col_est, col_std, nbsimu,
                flag_OK)) return 1;
  }

  /* Compute Conditional Expectation for Tonnage */

  if (selectivity->isUsed(ESelectivity::T))
  {
    if (st_ce_T(1, db, selectivity, iptr0, col_est, col_std, ycuts, nbsimu,
                flag_OK)) return 1;
  }

  /* Compute Conditional Expectation for Metal Quantity */

  if (selectivity->isUsed(ESelectivity::Q))
  {
    if (st_ce_Q(db, anam_hermite, selectivity, iptr0, col_est, col_std, ycuts,
                nbsimu, flag_OK)) return 1;
  }

  /* Compute Conditional Expectation for Conventional Benefit */

  if (selectivity->isUsed(ESelectivity::B))
  {
    if (st_ce_B(db, selectivity, iptr0, ycuts)) return 1;
  }

  /* Compute Conditional Expectation for Average recoveable grade */

  if (selectivity->isUsed(ESelectivity::M))
  {
    if (st_ce_M(db, selectivity, iptr0)) return 1;
  }

  /* Compute Conditional Expectation for Tonnage */

  if (selectivity->isUsed(ESelectivity::PROP) && need_T)
  {
    if (st_ce_T(2, db, selectivity, iptr0, col_est, col_std, ycuts, nbsimu,
                flag_OK)) return 1;
  }

  /* Compute Conditional Expectation for Quantile */

  if (selectivity->isUsed(ESelectivity::QUANT))
  {
    if (st_ce_quant(db, anam_hermite, selectivity, iptr0, col_est, col_std,
                    proba, flag_OK)) return 1;
  }

  if (! selectivity->isUsed(ESelectivity::T))
    (void) db_attribute_del_mult(db,
                                 selectivity->getAddressQTEst(ESelectivity::T, iptr0),
                                 selectivity->getNumberQTEst(ESelectivity::T));
  if (! selectivity->isUsed(ESelectivity::Q))
    (void) db_attribute_del_mult(db,
                                 selectivity->getAddressQTEst(ESelectivity::Q, iptr0),
                                 selectivity->getNumberQTEst(ESelectivity::Q));
  return 0;
}

/*****************************************************************************/
/*!
 **  Calculate the Uniform Conditioning
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  anam         Block Hermite anamorphosis
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  iptr0        Pointer for storage
 ** \param[in]  col_est      Rank of variable containing Kriging estimate
 ** \param[in]  col_var      Rank of Variable containing Variance of Kriging estimate
 **
 *****************************************************************************/
int _uniformConditioning(Db *db,
                         AnamHermite *anam,
                         Selectivity *selectivity,
                         int iptr0,
                         int col_est,
                         int col_var)
{
//  anam->setFlagBound(1);
  int nbpoly  = anam->getNbPoly();
  int ncuts = selectivity->getNCuts();

  /* Core allocation */

  VectorVectorDouble phi_b_zc(ncuts);

  /* Memorize the punctual Hermite polynomials */

  VectorDouble psi_hn = anam->getPsiHns();

  /* Add variables for storage */

  int iptr_sV = db->addColumnsByConstant(1, TEST);
  if (iptr_sV < 0) return 1;
  int iptr_yV = db->addColumnsByConstant(1, TEST);
  if (iptr_yV < 0) return 1;

  // Calculate the change of support coefficient

  double r_coef = anam->getRCoef();

  /* Transform zcuts into gaussian equivalent */

  VectorDouble ycuts = anam->rawToTransformVec(selectivity->getZcut());

  /* Fill the array phi_b_zc */

  for (int icut = 0; icut < ncuts; icut++)
  {
    VectorDouble hn = hermitePolynomials(ycuts[icut], 1., nbpoly);
    phi_b_zc[icut] = hermiteCoefMetal(ycuts[icut], hn);
  }

  /* Computing S and Y on panels */

  double vv_min =  1.e30;
  double zv_min =  1.e30;
  double vv_max = -1.e30;
  double zv_max = -1.e30;
  double sv_min =  1.e30;
  double yv_min =  1.e30;
  double sv_max = -1.e30;
  double yv_max = -1.e30;

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    anam->setPsiHns(psi_hn);
    anam->calculateMeanAndVariance();
    double zvstar = db->getArray(iech, col_est);
    double varv = db->getArray(iech, col_var);
    if (anam_point_to_block(anam, 0, varv, TEST, TEST)) continue;
    db->setArray(iech, iptr_sV, anam->getRCoef());
    db->setArray(iech, iptr_yV, anam->rawToTransformValue(zvstar));

    if (varv < vv_min) vv_min = varv;
    if (varv > vv_max) vv_max = varv;
    if (zvstar < zv_min) zv_min = zvstar;
    if (zvstar > zv_max) zv_max = zvstar;
  }

  /* Loop on the panels to compute the grade-tonnage functions */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    double sv = db->getArray(iech, iptr_sV);
    double yv = db->getArray(iech, iptr_yV);

    /* Loop on the cutoffs */

    for (int icut = 0; icut < ncuts; icut++)
    {
      double yc = ycuts[icut];
      double ore = 1. - law_cdf_gaussian(
              (yc - yv * sv / r_coef) / sqrt(1. - (sv / r_coef) * (sv / r_coef)));
      VectorDouble hn = hermitePolynomials(yv, 1., nbpoly);
      for (int ih = 0; ih < nbpoly; ih++)
        hn[ih] *= pow(sv / r_coef, (double) ih);
      double metal;
      matrix_product_safe(1, nbpoly, 1, hn.data(), phi_b_zc[icut].data(), &metal);

      if (sv < sv_min) sv_min = sv;
      if (sv > sv_max) sv_max = sv;
      if (yv < yv_min) yv_min = yv;
      if (yv > yv_max) yv_max = yv;

      /* Storing the grade-tonnage functions */

      selectivity->setTest(icut, ore);
      selectivity->setQest(icut, metal);
    }
    selectivity->calculateBenefitAndGrade();
    selectivity->storeInDb(db, iech, iptr0, TEST, TEST);
  }

  if (iptr_sV >= 0) db->deleteColumnByUID(iptr_sV);
  if (iptr_yV >= 0) db->deleteColumnByUID(iptr_yV);
  return 0;
}
