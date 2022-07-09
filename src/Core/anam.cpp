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
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Db/Db.hpp"

#include <math.h>
#include <Stats/SelectivityGlobal.hpp>

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

/****************************************************************************/
/*!
 **  Calculate the theoretical grade tonnage value
 **
 ** \return  Array f results (Dimension: 7 * nclass)
 **
 ** \param[in] anam         AAnam structure to be updated
 ** \param[in] zcut         Array of cutoffs
 ** \param[in] flag_correct 1 if Tonnage order relationship must be corrected
 ** \param[in] verbose      Verbose flag
 **
 ** \remark In the case of Discrete Anamorphosis, the number of classes
 ** \remark is defined by the number of cutoffs
 **
 *****************************************************************************/
SelectivityGlobal anam_selectivity(AAnam *anam,
                                   VectorDouble zcut,
                                   int flag_correct,
                                   int verbose)
{
  SelectivityGlobal calest;

  /* Dispatch according to the anamorphosis */

  if (anam->getType() == EAnam::HERMITIAN)
  {
    AnamHermite *anam_hermite = dynamic_cast<AnamHermite*>(anam);
    calest = anam_hermite->calculateSelectivity(zcut);
  }
  else if (anam->getType() == EAnam::DISCRETE_DD)
  {
    AnamDiscreteDD *anam_discrete_DD = dynamic_cast<AnamDiscreteDD*>(anam);
    calest = anam_discrete_DD->calculateSelectivity(flag_correct);
  }
  else if (anam->getType() == EAnam::DISCRETE_IR)
  {
    AnamDiscreteIR *anam_discrete_IR = dynamic_cast<AnamDiscreteIR*>(anam);
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
 **  Calculate the recoveries (z,T,Q,m,B) starting from the factors
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  anam         Point anamorphosis
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  cols_est     Array of columns for factor estimation
 ** \param[in]  cols_std     Array of columns for factor st. dev.
 ** \param[in]  z_max        Maximum grade array (only for QT interpolation)
 ** \param[in]  flag_correct 1 if Tonnage order relationship must be corrected
 ** \param[in]  verbose      Verbose flag
 **
 *****************************************************************************/
int anamFactor2QT(Db *db,
                  AAnam *anam,
                  const Selectivity* selectivity,
                  const VectorInt& cols_est,
                  const VectorInt& cols_std,
                  double z_max,
                  int flag_correct)
{
  int iptr = -1;
  AnamHermite *anam_hermite = dynamic_cast<AnamHermite*>(anam);
  AnamDiscreteDD *anam_discrete_DD = dynamic_cast<AnamDiscreteDD*>(anam);
  AnamDiscreteIR *anam_discrete_IR = dynamic_cast<AnamDiscreteIR*>(anam);
  int nb_est = (int) cols_est.size();
  int nb_std = (int) cols_std.size();

  /* Preliminary checks */

  if (db == nullptr)
  {
    messerr("You must define a Db");
    return 1;
  }
  if (anam == nullptr)
  {
    messerr("You must define an Anamorphosis");
    return 1;
  }
  if (nb_est <= 0 && nb_std <= 0)
  {
    messerr("The number of factors is zero");
    return 1;
  }

  int nvarout = selectivity->getVariableNumber();
  if (nvarout <= 0)
  {
    messerr("No recovery function is defined");
    return 1;
  }

  /* Variable allocation */

  iptr = db->addColumnsByConstant(nvarout, TEST);
  if (iptr < 0) return 1;

  /* Dispatch according to the type of Anamorphosis */

  switch (anam->getType().toEnum())
  {
    case EAnam::E_HERMITIAN:
      anam_hermite->factor2QT(db, selectivity, cols_est, cols_std, iptr);
      break;

    case EAnam::E_DISCRETE_DD:
      anam_discrete_DD->factor2QT(db, selectivity, cols_est, cols_std, iptr,
                                  z_max, flag_correct);
      break;

    case EAnam::E_DISCRETE_IR:
      anam_discrete_IR->factor2QT(db, selectivity, cols_est, cols_std, iptr,
                                  z_max, flag_correct);
      break;

    default:
      messerr("This method is not programmed yet for this anamorphosis");
      return 1;
  }

  return 0;
}

/*****************************************************************************/
/*!
 **  Calculate the Uniform Conditioning
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  anam         Point anamorphosis
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  att_est      Rank of the Kriging estimate
 ** \param[in]  att_var      Rank of the Variance of Kriging estimate
 ** \param[in]  zcuts        Array of the requested cutoffs
 ** \param[in]  var_bloc     Change of support coefficient
 ** \param[in]  verbose      Verbose option
 **
 *****************************************************************************/
int uc(Db *db,
       AAnam *anam,
       const Selectivity* selectivity,
       int att_est,
       int att_var,
       const VectorDouble& zcuts,
       double var_bloc,
       int verbose)
{
  int error, nbpoly, iptr, iptr_sV, iptr_yV, nvarout, ncuts;
  double yc, sv, yv, ore, metal, mean, variance, varb, r_coef, zvstar, varv;
  double vv_min, vv_max, sv_min, sv_max, zv_min, zv_max, yv_min, yv_max;
  VectorDouble psi_hn, hn, ycuts;
  VectorVectorDouble phi_b_zc;
  SelectivityGlobal calest;

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
  if (zcuts.empty())
  {
    messerr("You must define some cutoff values");
    goto label_end;
  }
  nbpoly = anam_hermite->getNbPoly();
  mean = anam_hermite->getMean();
  variance = anam_hermite->getVariance();
  ncuts = (int) zcuts.size();

  /* Core allocation */

  psi_hn.resize(nbpoly);
  phi_b_zc.resize(ncuts);

  /* Memorize the punctual Hermite polynomials */

  for (int ip = 0; ip < nbpoly; ip++)
    psi_hn[ip] = anam_hermite->getPsiHn(ip);

  /* Add variables for storage */

  iptr_sV = db->addColumnsByConstant(1, TEST);
  if (iptr_sV < 0) goto label_end;
  iptr_yV = db->addColumnsByConstant(1, TEST);
  if (iptr_yV < 0) goto label_end;

  /* Analyzing the codes */

  nvarout = selectivity->getVariableNumber();
  if (nvarout <= 0) goto label_end;
  if (selectivity->isUsed(ESelectivity::Z))
  {
    messerr("The recovery option 'Z' is not available in this function");
    goto label_end;
  }
  iptr = db->addColumnsByConstant(nvarout, TEST);
  if (iptr < 0) goto label_end;

  /* Core allocation */

  calest = SelectivityGlobal(ncuts);

  /* Transforming Point anamorphosis into Block Anamorphosis */

  if (anam_point_to_block(anam, 0, var_bloc, TEST, TEST)) goto label_end;
  r_coef = anam_hermite->getRCoef();
  varb = anam_hermite->getVariance();

  /* Transform cutmine into gaussian equivalent */

  ycuts = zcuts;
  for (int icut = 0; icut < ncuts; icut++)
    ycuts[icut] = anam_hermite->RawToTransformValue(zcuts[icut]);

  /* Fill the array phi_b_zc */

  for (int icut = 0; icut < ncuts; icut++)
  {
    hn = hermitePolynomials(ycuts[icut], 1., nbpoly);
    phi_b_zc[icut] = hermiteCoefMetal(ycuts[icut], hn);
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

    for (int icut = 0; icut < ncuts; icut++)
    {
      yc = ycuts[icut];
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

    calest.calculateBenefitAndGrade();
    anam_hermite->recoveryLocal(db, selectivity, iech, iptr, TEST, TEST, calest);
  }

  /* Verbose printout (optional) */

  if (verbose)
  {
    message("Uniform Conditioning on %d panels and %d cutoffs\n",
            db->getSampleNumber(true), ncuts);
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
                               bool /*flag_OK*/,
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
 ** \param[in]  anam         Anamorphosis Hermite
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  iptr0        Starting storage address
 ** \param[in]  att_est      Rank of the Kriging estimate
 ** \param[in]  att_std      Rank of the St, Deviation of Kriging estimate
 ** \param[in]  nbsimu       Number of Monte Carlo simulations (0: Hermite)
 ** \param[in]  flag_OK      1 if kriging is performed with OK
 **
 *****************************************************************************/
static int st_ce_compute_Z(Db *db,
                           const AnamHermite* anam,
                           const Selectivity* selectivity,
                           int iptr0,
                           int att_est,
                           int att_std,
                           int nbsimu,
                           bool flag_OK)
{
  VectorDouble krigest, krigstd, valest, valstd;

  st_ce_get_vectors(db, att_est, att_std, flag_OK, krigest, krigstd);

  if (nbsimu <= 0)
  {
    valest = hermiteCondExp(krigest, krigstd, anam->getPsiHn());
    valstd = hermiteCondStd(krigest, krigstd, anam->getPsiHn());
  }
  else
  {
    valest = MCCondExp(krigest, krigstd, anam->getPsiHn(), nbsimu);
    valstd = MCCondStd(krigest, krigstd, anam->getPsiHn(), nbsimu);
  }

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    if (selectivity->isUsedEst(ESelectivity::Z))
      db->setArray(iech,
                   selectivity->getAddressQTEst(ESelectivity::Z, iptr0), valest[iech]);
    if (selectivity->isUsedStD(ESelectivity::Z))
      db->setArray(iech,
                   selectivity->getAddressQTStD(ESelectivity::Z, iptr0), valstd[iech]);
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
 ** \arams[in]  iptr0        Starting storage address
 ** \param[in]  att_est      Rank of the Kriging estimate
 ** \param[in]  att_std      Rank of the St, Deviation of Kriging estimate
 ** \param[in]  ycuts        Vector of Gaussian cutoffs
 ** \param[in]  nbsimu       Number of Monte Carlo simulations (0: Hermite)
 ** \param[in]  flag_OK      1 if kriging is performed with OK
 **
 *****************************************************************************/
static int st_ce_compute_T(int mode,
                           Db *db,
                           const Selectivity* selectivity,
                           int iptr0,
                           int att_est,
                           int att_std,
                           const VectorDouble& ycuts,
                           int nbsimu,
                           bool flag_OK)
{
  VectorDouble krigest, krigstd, valest, valstd;

  /* Loop on the samples */

  for (int icut = 0; icut < selectivity->getNCuts(); icut++)
  {
    st_ce_get_vectors(db, att_est, att_std, flag_OK, krigest, krigstd);

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

    for (int iech = 0; iech < db->getSampleNumber(); iech++)
    {
      if (!db->isActive(iech)) continue;
      if (mode == 1)
      {
        if (selectivity->isUsedEst(ESelectivity::T))
          db->setArray(
              iech,
              selectivity->getAddressQTEst(ESelectivity::T, iptr0, icut),
              valest[iech]);
        if (selectivity->isUsedStD(ESelectivity::T))
           db->setArray(iech,
                        selectivity->getAddressQTStD(ESelectivity::T, iptr0, icut),
                        valstd[iech]);
      }
      else
      {
        if (selectivity->isUsedEst(ESelectivity::PROBA))
          db->setArray(
              iech,
              selectivity->getAddressQTEst(ESelectivity::PROBA, iptr0, icut),
              1. - valest[iech]);
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
 ** \param[in]  proba        Probability threshold
 ** \param[in]  att_std      Rank of the St, Deviation of Kriging estimate
 ** \param[in]  proba        Probability threshold
 ** \param[in]  flag_OK      1 if kriging is performed with OK
 **
 *****************************************************************************/
static int st_ce_compute_quant(Db *db,
                               const AnamHermite* anam,
                               const Selectivity* selectivity,
                               int iptr0,
                               int att_est,
                               int att_std,
                               double proba,
                               bool flag_OK)
{
  double krigest, krigstd;

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    st_correct_from_OK(db, iech, att_est, att_std, flag_OK, &krigest, &krigstd);
    double y = krigest + krigstd * law_invcdf_gaussian(proba);
    db->setArray(iech,
                 selectivity->getAddressQTEst(ESelectivity::QUANT, iptr0),
                 anam->TransformToRawValue(y));
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
 ** \param[in]  iptr         Starting storage address
 ** \param[in]  att_est      Rank of the Kriging estimate
 ** \param[in]  att_std      Rank of the St, Deviation of Kriging estimate
 ** \param[in]  ycuts        Array of the requested cutoffs (gaussian scale)
 ** \param[in]  nbsimu       Number of Monte Carlo simulations (0: Hermite)
 ** \param[in]  flag_OK      1 if kriging is performed with OK
 **
 *****************************************************************************/
static int st_ce_compute_Q(Db *db,
                           const AnamHermite* anam,
                           const Selectivity* selectivity,
                           int iptr0,
                           int att_est,
                           int att_std,
                           const VectorDouble& ycuts,
                           int nbsimu,
                           bool flag_OK)
{
  VectorDouble krigest, krigstd, valest, valstd;

  /* Loop on the cutoffs */

  for (int icut = 0; icut < selectivity->getNCuts(); icut++)
  {
    st_ce_get_vectors(db, att_est, att_std, flag_OK, krigest, krigstd);

    if (nbsimu <= 0)
    {
      valest = hermiteMetal(ycuts[icut], krigest, krigstd, anam->getPsiHn());
      valstd = hermiteMetalStd(ycuts[icut], krigest, krigstd, anam->getPsiHn());
    }
    else
    {
      valest = MCMetal(ycuts[icut], krigest, krigstd, anam->getPsiHn(), nbsimu);
      valstd = MCMetalStd(ycuts[icut], krigest, krigstd, anam->getPsiHn(), nbsimu);
    }

    for (int iech = 0; iech < db->getSampleNumber(); iech++)
    {
      if (!db->isActive(iech)) continue;
      if (selectivity->isUsedEst(ESelectivity::Q))
        db->setArray(iech,
            selectivity->getAddressQTEst(ESelectivity::Q, iptr0, icut),
            valest[iech]);
      if (selectivity->isUsedStD(ESelectivity::Q))
        db->setArray(iech,
                     selectivity->getAddressQTStD(ESelectivity::Q, iptr0, icut),
                     valstd[iech]);
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
static int st_ce_compute_B(Db *db,
                           const Selectivity* selectivity,
                           int iptr0,
                           const VectorDouble& ycuts)
{
  double t, q, b;

  /* Loop on the samples */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;

    /* Loop on the cutoffs */

    for (int icut = 0; icut < selectivity->getNCuts(); icut++)
    {
      t = db->getArray(iech,
                       selectivity->getAddressQTEst(ESelectivity::T, iptr0, icut));
      q = db->getArray(iech,
                       selectivity->getAddressQTEst(ESelectivity::Q, iptr0, icut));
      b = q - t * ycuts[icut];
      db->setArray(iech,
                   selectivity->getAddressQTEst(ESelectivity::B, iptr0, icut), b);
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
static int st_ce_compute_M(Db *db,
                           const Selectivity* selectivity,
                           int iptr0)
{
  double t, q, m;

  /* Loop on the samples */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;

    /* Loop on the cutoffs */

    for (int icut = 0; icut < selectivity->getNCuts(); icut++)
    {
      t = db->getArray(iech,
                       selectivity->getAddressQTEst(ESelectivity::T, iptr0, icut));
      q = db->getArray(iech,
                       selectivity->getAddressQTEst(ESelectivity::Q, iptr0, icut));
      m = (t > EPSILON3) ? q / t : TEST;
      db->setArray(iech,
                   selectivity->getAddressQTEst(ESelectivity::M, iptr0, icut), m);
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
 ** \param[in]  anam         Point anamorphosis
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  zcuts        Array of the requested cutoffs
 ** \param[in]  att_est      Rank of the Kriging estimate
 ** \param[in]  att_std      Rank of the St, Deviation of Kriging estimate
 ** \param[in]  flag_est     1 for computing the Estimation
 ** \param[in]  flag_std     1 for computing the St. Deviation
 ** \param[in]  flag_OK      1 if kriging has ben performed in Ordinary Kriging
 ** \param[in]  proba        Probability
 ** \param[in]  verbose      Verbose option
 **
 *****************************************************************************/
int ce(Db *db,
       AAnam *anam,
       const Selectivity* selectivity,
       int att_est,
       int att_std,
       bool flag_est,
       bool flag_std,
       bool flag_OK,
       double proba,
       int nbsimu,
       int verbose)
{
  AnamHermite *anam_hermite = dynamic_cast<AnamHermite*>(anam);
  int ncuts = selectivity->getNCuts();

  if (anam->getType() != EAnam::HERMITIAN)
  {
    messerr("The argument 'anam' must be Gaussian");
    return 1;
  }
  int nbpoly = anam_hermite->getNbPoly();

  /* Analyzing the codes */

  int count = 0;
  if (flag_est) count++;
  if (flag_std) count++;
  VectorDouble ycuts = anam_hermite->RawToTransformVec(selectivity->getZcut());
  int need_T = selectivity->isNeededT();
  int nvarout = selectivity->getVariableNumber();

  /* Add the variables */

  int iptr = db->addColumnsByConstant(nvarout, TEST);
  if (iptr < 0) return 1;

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
    message(" Number of cutoffs                  : %6d\n", ncuts);
    if (nbsimu > 0)
      message(" Number of Monte-Carlo simulations  : %6d\n", nbsimu);
  }

  /* Computing the estimation */

  const AnamHermite* anamH = dynamic_cast<const AnamHermite*>(anam);
  if (selectivity->isUsed(ESelectivity::Z))
  {
    if (st_ce_compute_Z(db, anamH, selectivity, iptr, att_est, att_std,
                        nbsimu, flag_OK)) return 1;
  }

  /* Compute Conditional Expectation for Tonnage */

  if (selectivity->isUsed(ESelectivity::T))
  {
    if (st_ce_compute_T(1, db, selectivity, iptr, att_est, att_std,
                        ycuts, nbsimu, flag_OK)) return 1;
  }

  /* Compute Conditional Expectation for Metal Quantity */

  if (selectivity->isUsed(ESelectivity::Q))
  {
    if (st_ce_compute_Q(db, anamH, selectivity, iptr, att_est, att_std,
                        ycuts, nbsimu, flag_OK)) return 1;
  }

  /* Compute Conditional Expectation for Conventional Benefit */

  if (selectivity->isUsed(ESelectivity::B))
  {
    if (st_ce_compute_B(db, selectivity, iptr, ycuts)) return 1;
  }

  /* Compute Conditional Expectation for Average recoveable grade */

  if (selectivity->isUsed(ESelectivity::M))
  {
    if (st_ce_compute_M(db, selectivity, iptr)) return 1;
  }

  /* Compute Conditional Expectation for Tonnage */

  if (selectivity->isUsed(ESelectivity::PROBA) && need_T)
  {
    if (st_ce_compute_T(2, db, selectivity, iptr, att_est, att_std,
                        ycuts, nbsimu, flag_OK)) return 1;
  }

  /* Compute Conditional Expectation for Quantile */

  if (selectivity->isUsed(ESelectivity::QUANT))
  {
    if (st_ce_compute_quant(db, anamH, selectivity, iptr, att_est, att_std,
                            proba, flag_OK)) return 1;
  }

  if (! selectivity->isUsed(ESelectivity::T))
    (void) db_attribute_del_mult(db,
                                 selectivity->getAddressQTEst(ESelectivity::T, iptr),
                                 selectivity->getNumberQTEst(ESelectivity::T));
  if (! selectivity->isUsed(ESelectivity::Q))
    (void) db_attribute_del_mult(db,
                                 selectivity->getAddressQTEst(ESelectivity::Q, iptr),
                                 selectivity->getNumberQTEst(ESelectivity::Q));
  return 0;
}

