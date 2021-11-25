/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "Basic/Utilities.hpp"
#include "LithoRule/RuleProp.hpp"
#include "LithoRule/PropDef.hpp"
#include "LithoRule/Rule.hpp"
#include "LithoRule/RuleShift.hpp"
#include "LithoRule/RuleShadow.hpp"
#include "geoslib_e.h"
#include "geoslib_old_f.h"

/*! \cond */
#define INDLOC(ifac1,ifac2)  ((ifac2)*propdef->nfac[0]+(ifac1))
#define PROPFIX(ifac1,ifac2) (propdef->propfix[INDLOC(ifac1,ifac2)])
#define PROPWRK(ifac1,ifac2) (propdef->propwrk[INDLOC(ifac1,ifac2)])
/*! \endcond */

/****************************************************************************/
/*!
 **  Free a Rule structure
 **
 ** \return  Pointer to the newly freed Rule structure
 **
 ** \param[in]  rule Rule structure to be freed
 **
 *****************************************************************************/
GSTLEARN_EXPORT Rule* rule_free(const Rule *rule)

{
  if (rule != nullptr) delete rule;
  return (nullptr);
}

/****************************************************************************/
/*!
 **  Locate the current proportions
 **
 ** \param[in]  propdef    PropDef structure
 ** \param[in]  ifac_ref   Conditional (first variable) facies
 **                        (Only used for EProcessOper::CONDITIONAL)
 **
 ****************************************************************************/
static int st_proportion_locate(PropDef *propdef, int ifac_ref)
{
  int ifac;

  switch (propdef->mode.toEnum())
  {
    case EProcessOper::E_COPY:
    case EProcessOper::E_MARGINAL:
      for (ifac = 0; ifac < propdef->nfaccur; ifac++)
        propdef->proploc[ifac] = propdef->propwrk[ifac];
      break;

    case EProcessOper::E_CONDITIONAL:
      for (ifac = 0; ifac < propdef->nfaccur; ifac++)
        propdef->proploc[ifac] = PROPWRK(ifac_ref - 1, ifac);
      break;

    default:
      messerr("Unknown process operation");
      break;
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Transform the proportions (from CST to WRK)
 **
 ** \return  -1 if the proportion is not defined; 0 otherwise
 **
 ** \param[in]  propdef   PropDef structure
 **
 ****************************************************************************/
static int st_proportion_transform(PropDef *propdef)

{
  double total, pp;
  int ifac1, ifac2;

  /* Dispatch */

  switch (propdef->mode.toEnum())
  {
    case EProcessOper::E_COPY:
      for (ifac1 = 0; ifac1 < propdef->nfac[0]; ifac1++)
      {
        pp = PROPFIX(ifac1, 0);
        if (FFFF(pp)) return (-1);
        PROPWRK(ifac1,0) = pp;
      }
      break;

    case EProcessOper::E_MARGINAL:
      for (ifac1 = 0; ifac1 < propdef->nfac[0]; ifac1++)
      {
        PROPWRK(ifac1,0) = 0.;
        for (ifac2 = 0; ifac2 < propdef->nfac[1]; ifac2++)
        {
          pp = PROPFIX(ifac1, ifac2);
          if (FFFF(pp)) return (-1);
          PROPWRK(ifac1,0) += pp;
        }
      }
      break;

    case EProcessOper::E_CONDITIONAL:
      for (ifac1 = 0; ifac1 < propdef->nfac[0]; ifac1++)
        for (ifac2 = 0; ifac2 < propdef->nfac[1]; ifac2++)
        {
          pp = PROPFIX(ifac1, ifac2);
          if (FFFF(pp)) return (-1);
          PROPWRK(ifac1,ifac2) = pp;
        }

      for (ifac1 = 0; ifac1 < propdef->nfac[0]; ifac1++)
      {
        total = 0.;
        for (ifac2 = 0; ifac2 < propdef->nfac[1]; ifac2++)
          total += PROPWRK(ifac1, ifac2);

        for (ifac2 = 0; ifac2 < propdef->nfac[1]; ifac2++)
          PROPWRK(ifac1,ifac2) = (total <= 0.) ? 1. / propdef->nfac[1] :
                                                 PROPWRK(ifac1,ifac2) / total;
      }
      break;

    default:
      messageAbort("This should never happen in st_proportion_transform");
      break;
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Set the method to compute Proportions
 **
 ** \param[in]  propdef  PropDef structure
 ** \param[in]  mode     Type of operation (EProcessOper)
 **
 ****************************************************************************/
GSTLEARN_EXPORT void proportion_rule_process(PropDef *propdef,
                                             const EProcessOper &mode)
{
  /* Assignments */

  propdef->mode = mode;

  /* Assign the current value for the number of facies */

  if (mode == EProcessOper::COPY || mode == EProcessOper::MARGINAL)
    propdef->nfaccur = propdef->nfac[0];
  if (mode == EProcessOper::CONDITIONAL) propdef->nfaccur = propdef->nfac[1];

  /* In the stationary case, transform the proportions (from CST to WRK) */

  if (propdef->case_stat) st_proportion_transform(propdef);

  return;
}

/****************************************************************************/
/*!
 **  Print the (non-stationary) proportions
 **
 ** \param[in]  propdef   PropDef structure
 **
 *****************************************************************************/
GSTLEARN_EXPORT void proportion_print(PropDef *propdef)

{
  if (propdef == nullptr) return;
  mestitle(0, "Proportions");

  print_matrix("Initial :", 0, 1, propdef->nfac[1], propdef->nfac[0], NULL,
               propdef->propfix.data());

  print_matrix("Working :", 0, 1, propdef->nfac[1], propdef->nfac[0], NULL,
               propdef->propwrk.data());

  print_matrix("Current :", 0, 1, propdef->nfaccur, 1, NULL,
               propdef->proploc.data());
}

/****************************************************************************/
/*!
 **  Check if the proportion has changed since the previous usage
 **  and store the current proportions for future comparison
 **
 ** \return  1 if the proportions are unchanged; 0 otherwise
 **
 ** \param[in]  propdef PropDef structure
 **
 ****************************************************************************/
static int st_proportion_changed(PropDef *propdef)

{
  /* Compare with the memory proportion array */

  int modify = !ut_vector_same(propdef->proploc, propdef->propmem);
  if (!modify) return (1);

  /* Print the proportions (optional) */

  if (debug_query("props")) proportion_print(propdef);

  for (int ifac = 0; ifac < propdef->nfaccur; ifac++)
    propdef->propmem[ifac] = propdef->proploc[ifac];

  return (0);
}

/****************************************************************************/
/*!
 **  Set the (non-stationary) proportions
 **
 ** \return  Error return code
 ** \return  - the target point does not lie within the proportion grid
 ** \return  - in conditional processing, the reference facies does not exist
 **
 ** \param[in]  propdef    PropDef structure
 ** \param[in]  db         Db input structure
 ** \param[in]  iech       Rank of the data in the input Db
 ** \param[in]  isimu      Rank of the simulation (EProcessOper::CONDITIONAL)
 ** \param[in]  nbsimu     Number of simulations
 **
 ** \param[out] jech       Rank of the auxiliary data in the input Db
 **
 ** \remark  At the end of this function, the local proportions are stored
 ** \remark  in the array proploc of the structure PropDef
 ** \remark  The argument 'isimu' is only used for
 ** \remark            propdef->mode == EProcessOper::CONDITIONAL (simbipgs)
 **
 *****************************************************************************/
static int st_proportion_define(PropDef *propdef,
                                const Db *db,
                                int iech,
                                int isimu,
                                int nbsimu,
                                int *jech)
{
  int ifac, ifac_ref;

  /* Non-stationary case : Load the proportions in propcst */

  (*jech) = 0;
  if (!propdef->case_stat)
  {
    if (propdef->case_prop_interp)
    {

      /* Case where the proportions must be interpolated */

      (*jech) = index_point_to_grid(db, iech, 1, propdef->dbprop,
                                    propdef->coor.data());
      if ((*jech) < 0)
      {
        messerr("At the data #%d, the proportion matrix is undefined",
                iech + 1);
        return (1);
      }

      /* Load the proportions (into CST) */

      for (ifac = 0; ifac < propdef->nfacprod; ifac++)
        propdef->propfix[ifac] = propdef->dbprop->getProportion(*jech, ifac);
    }
    else
    {

      /* The proportions are already available from the dbin */

      for (ifac = 0; ifac < propdef->nfacprod; ifac++)
        propdef->propfix[ifac] = db->getProportion(iech, ifac);
    }

    /* Transform proportions (from CST to WRK) */

    st_proportion_transform(propdef);
  }

  /* Locate the current proportions (from WRK to LOC) */

  ifac_ref = -1;
  if (propdef->mode == EProcessOper::CONDITIONAL)
  {
    ifac_ref = (int) db->getSimvar(ELoc::FACIES, iech, isimu, 0, 0, nbsimu, 1);
    if (ifac_ref < 1 || ifac_ref > propdef->nfac[0]) return (1);
  }
  st_proportion_locate(propdef, ifac_ref);

  return (0);
}

/****************************************************************************/
/*!
 **  Set the (non-stationary) proportions and define thresholds (for shadow only)
 **
 ** \return  Error return code
 **
 ** \param[in]  propdef    PropDef structure
 ** \param[in]  db         Db input structure
 ** \param[in]  rule       Rule structure
 ** \param[in]  facies     Facies of interest (or GV_ITEST)
 ** \param[in]  iech       Rank of the data in the input Db
 ** \param[in]  isimu      Rank of the simulation (EProcessOper::CONDITIONAL)
 ** \param[in]  nbsimu     Number of simulations (EProcessOper::CONDITIONAL)
 **
 ** \param[out] t1min      Minimum threshold for Y1
 ** \param[out] t1max      Maximum threshold for Y1
 ** \param[out] t2min      Minimum threshold for Y2
 ** \param[out] t2max      Maximum threshold for Y2
 ** \param[out] sh_dsup    Local or global upwards shift (shadow)
 ** \param[out] sh_down    Local or global downwards shift (shadow)
 **
 *****************************************************************************/
GSTLEARN_EXPORT int rule_thresh_define_shadow(PropDef *propdef,
                                              Db *db,
                                              const RuleShadow *rule,
                                              int facies,
                                              int iech,
                                              int isimu,
                                              int nbsimu,
                                              double *t1min,
                                              double *t1max,
                                              double *t2min,
                                              double *t2max,
                                              double *sh_dsup,
                                              double *sh_down)
{
  int unmodify, facloc, jech;

  /* Set the debugging information */

  debug_index(iech + 1);

  /* Processing an "unknown" facies */

  if (!IFFFF(facies) && (facies < 1 || facies > propdef->nfaccur))
  {
    *t1min = *t2min = get_rule_extreme(-1);
    *t1max = *t2max = get_rule_extreme(+1);
    return (0);
  }

  /* Define the proportions */

  if (st_proportion_define(propdef, db, iech, isimu, nbsimu, &jech))
  {
    *t1min = *t2min = get_rule_extreme(-1);
    *t1max = *t2max = get_rule_extreme(+1);
    return (0);
  }

  /* Check if the proportions have been changed */

  unmodify = st_proportion_changed(propdef);

  /* In case of Shadow, return the upwards and downwards values */

  *sh_dsup = (propdef->case_stat) ? rule->getShDsup() :
                                    propdef->proploc[1];
  *sh_down = (propdef->case_stat) ? rule->getShDown() :
                                    propdef->proploc[2];

  /* In the special cases, only the first proportion is significant */

  propdef->proploc[1] = (1 - propdef->proploc[0]) / 2;
  propdef->proploc[2] = (1 - propdef->proploc[0]) / 2;

  /* Set the proportions and translate proportions into thresholds */

  if (!unmodify)
  {
    if (rule->setProportions(propdef->proploc)) return (1);
  }

  /* Convert the proportions into thresholds */

  facloc = (IFFFF(facies)) ? 1 :
                             facies;
  VectorDouble bounds = rule->getThresh(facloc);
  *t1min = bounds[0];
  *t1max = bounds[1];
  *t2min = bounds[2];
  *t2max = bounds[3];

  return (0);
}

/****************************************************************************/
/*!
 **  Set the (non-stationary) proportions and define thresholds
 **
 ** \return  Error return code
 **
 ** \param[in]  propdef    PropDef structure
 ** \param[in]  db         Db input structure
 ** \param[in]  rule       Rule structure
 ** \param[in]  facies     Facies of interest (or ITEST) starting from 1
 ** \param[in]  iech       Rank of the data in the input Db
 ** \param[in]  isimu      Rank of the simulation (EProcessOper::CONDITIONAL)
 ** \param[in]  nbsimu     Number of simulations
 ** \param[in]  flag_check 1 if the consistency check with the actual
 **                        proportion of the current facies must be done
 **
 ** \param[out] t1min      Minimum threshold for Y1
 ** \param[out] t1max      Maximum threshold for Y1
 ** \param[out] t2min      Minimum threshold for Y2
 ** \param[out] t2max      Maximum threshold for Y2
 **
 *****************************************************************************/
GSTLEARN_EXPORT int rule_thresh_define(PropDef *propdef,
                                       Db *db,
                                       const Rule *rule,
                                       int facies,
                                       int iech,
                                       int isimu,
                                       int nbsimu,
                                       int flag_check,
                                       double *t1min,
                                       double *t1max,
                                       double *t2min,
                                       double *t2max)
{
  int unmodify, facloc, jech;

  /* Set the debugging information */

  debug_index(iech + 1);

  /* Processing an "unknown" facies */

  if (!IFFFF(facies) && (facies < 1 || facies > propdef->nfaccur))
  {
    *t1min = *t2min = get_rule_extreme(-1);
    *t1max = *t2max = get_rule_extreme(+1);
    return (0);
  }

  /* Define the proportions */

  if (st_proportion_define(propdef, db, iech, isimu, nbsimu, &jech))
  {
    *t1min = *t2min = get_rule_extreme(-1);
    *t1max = *t2max = get_rule_extreme(+1);
    return (0);
  }

  /* Check if the proportions have been changed */

  unmodify = st_proportion_changed(propdef);

  /* Check that the facies is compatible with the proportions */

  if (flag_check && !IFFFF(facies) && rule->getModeRule() == ERule::STD)
  {
    if (propdef->proploc[facies - 1] <= 0.)
    {
      messerr(
          "The presence of facies (%d) at sample (%d) is not consistent with the zero proportion",
          facies, iech + 1);
      if (!propdef->case_stat)
        messerr("Check the proportions in the cell (%d) of the Proportion Db",
                jech + 1);
      return (1);
    }
  }

  /* Set the proportions and translate proportions into thresholds */

  if (!unmodify)
  {
    if (rule->setProportions(propdef->proploc)) return (1);

    /* In the case of SHIFT, update the thresholds */

    if (rule->getModeRule() == ERule::SHIFT && 0) rule->updateShift();
  }

  /* Convert the proportions into thresholds */

  facloc = (IFFFF(facies)) ? 1 :
                             facies;
  VectorDouble bounds = rule->getThresh(facloc);
  *t1min = bounds[0];
  *t1max = bounds[1];
  *t2min = bounds[2];
  *t2max = bounds[3];

  return (0);
}

/****************************************************************************/
/*!
 **  Apply the Rule transformation to the GRFs of a Db
 **  (Shadow case)
 **
 ** \return  Error return code
 **
 ** \param[in]  db        Output Db structure
 ** \param[in]  dbprop    Db structure used for proportions (non-stationary case)
 ** \param[in]  rule      Lithotype Rule definition
 ** \param[in]  model     First Model structure (only for SHIFT)
 ** \param[in]  props     Array of proportions for the facies
 ** \param[in]  flag_stat 1 for stationary; 0 otherwise
 ** \param[in]  nfacies   Number of facies
 **
 ** \remark The input variable must be locatorized as Z or ELoc::SIMU
 ** \remark It will be changed in this function to locator ELoc::SIMU
 **
 *****************************************************************************/
GSTLEARN_EXPORT int db_rule_shadow(Db *db,
                                   Db *dbprop,
                                   RuleShadow *rule,
                                   Model *model,
                                   const VectorDouble &props,
                                   int flag_stat,
                                   int nfacies)
{
  int iptr, error, flag_used[2], nbsimu, igrf, ngrf;
  PropDef *propdef;

  /* Initializations */

  error = 1;
  nbsimu = 1;
  iptr = -1;
  propdef = nullptr;

  /* Preliminary checks */

  ngrf = rule->getGRFNumber();
  for (igrf = 0; igrf < 2; igrf++)
    flag_used[igrf] = rule->isYUsed(igrf);

  propdef = proportion_manage(1, 1, flag_stat, ngrf, 0, nfacies, 0, db, dbprop,
                              props, propdef);
  if (propdef == nullptr) goto label_end;

  /* General setting for lithotype */

  rule->particularities(db, dbprop, model, 1, flag_stat);
  proportion_rule_process(propdef, EProcessOper::COPY);

  /**********************/
  /* Add the attributes */
  /**********************/

  /* Storage of the simulations in the output file */
  iptr = db->addFields(nbsimu, 0.);
  if (iptr < 0) goto label_end;
  db->setLocatorsByAttribute(nbsimu, iptr, ELoc::FACIES);

  /* Identify the Non conditional simulations at target points */
  for (igrf = 0; igrf < 2; igrf++)
  {
    if (!flag_used[igrf]) continue;
    iptr = db_attribute_identify(db, ELoc::SIMU, igrf);
    if (iptr < 0)
    {
      iptr = db_attribute_identify(db, ELoc::Z, igrf);
      if (iptr < 0)
      {
        messerr(
            "The variable containing the simulation of the GRF %d is missing in the Db",
            igrf + 1);
        goto label_end;
      }
      db->setLocatorByAttribute(iptr, ELoc::SIMU, igrf);
    }
  }

  /* Combine the conditional simulation for each GRF */

  for (int isimu = 0; isimu < nbsimu; isimu++)
    if (rule->gaus2facResult(propdef, db, flag_used, 0, isimu, nbsimu))
      goto label_end;

  /* Set the error return flag */

  error = 0;

  label_end: propdef = proportion_manage(-1, 1, flag_stat, ngrf, 0, nfacies, 0,
                                         db, dbprop, props, propdef);
  return (error);
}

/****************************************************************************/
/*!
 **  Apply the Rule transformation to convert a set of Gaussian vectors
 **  into the corresponding Facies in a Db
 **
 ** \return  Error return code
 **
 ** \param[in]  db        Output Db structure
 ** \param[in]  ruleprop  RuleProp structure
 ** \param[in]  model     First Model structure (only for SHIFT)
 ** \param[in]  namconv   Naming convention
 **
 ** \remark The input variable must be locatorized as Z or ELoc::SIMU
 **
 *****************************************************************************/
int _db_rule(Db *db,
             const RuleProp *ruleprop,
             Model *model,
             NamingConvention namconv)
{
  if (db == nullptr)
  {
    messerr("The Db is not defined");
    return 1;
  }
  if (ruleprop == nullptr)
  {
    messerr("RuleProp must be defined");
    return 1;
  }
  int flag_stat = ruleprop->isFlagStat();
  const Rule *rule = ruleprop->getRule();
  const VectorDouble &propcst = ruleprop->getPropCst();
  const Db *dbprop = ruleprop->getDbprop();

  int error = 1;
  int iptr = -1;
  PropDef *propdef = nullptr;
  int ngrf = rule->getGRFNumber();
  VectorInt flagUsed = rule->whichGRFUsed();
  int nfacies = rule->getFaciesNumber();
  bool flagReturn = false;

  /* Preliminary checks */

  if (db->getLocatorNumber(ELoc::SIMU) != ngrf && db->getLocatorNumber(ELoc::Z)
      != ngrf)
  {
    messerr("The Rule specifies the use of %d underlying GRF(s)", ngrf);
    messerr(
        "The input 'db' should have one variable per GRF with locator 'SIMU' or 'Z'");
    goto label_end;
  }

  propdef = proportion_manage(1, 1, flag_stat, ngrf, 0, nfacies, 0, db, dbprop,
                              propcst, propdef);
  if (propdef == nullptr) goto label_end;
  if (rule->particularities(db, dbprop, model, 1, flag_stat)) goto label_end;
  proportion_rule_process(propdef, EProcessOper::COPY);

  /**********************/
  /* Add the attributes */
  /**********************/

  /* Storage of the simulations in the output file */
  iptr = db->addFields(1, 0., "Facies", ELoc::FACIES);
  if (iptr < 0) goto label_end;

  /* Identify the Non conditional simulations at target points */

  if (db->getLocatorNumber(ELoc::SIMU) != ngrf)
  {
    db->switchLocator(ELoc::Z, ELoc::SIMU);
    flagReturn = true;
  }

  /* Translate Gaussian into Facies */

  if (rule->gaus2facResult(propdef, db, flagUsed.data(), 0, 0, 1))
    goto label_end;

  // Returning to the initial locators (if the initial variable
  // had a ELoc::Z locator which has been temporarily modified into ELoc::SIMU)

  if (flagReturn) db->switchLocator(ELoc::SIMU, ELoc::Z);

  // Naming convention

  namconv.setNamesAndLocators(nullptr, VectorInt(), db, iptr);
  error = 0;

  label_end: propdef = proportion_manage(-1, 1, flag_stat, ngrf, 0, nfacies, 0,
                                         db, dbprop, propcst, propdef);
  return (error);
}

/****************************************************************************/
/*!
 **  Apply the Rule transformation to derive the bounds variables for a Db
 **  (Shadow case)
 **
 ** \return  Error return code
 **
 ** \param[in]  db        Db structure
 ** \param[in]  dbprop    Db structure used for proportions (non-stationary case)
 ** \param[in]  rule      Lithotype Rule definition
 ** \param[in]  model     First Model structure (only for SHIFT)
 ** \param[in]  props     Array of proportions for the facies
 ** \param[in]  flag_stat 1 for stationary; 0 otherwise
 ** \param[in]  nfacies   Number of facies
 **
 *****************************************************************************/
GSTLEARN_EXPORT int db_bounds_shadow(Db *db,
                                     Db *dbprop,
                                     RuleShadow *rule,
                                     Model *model,
                                     const VectorDouble &props,
                                     int flag_stat,
                                     int nfacies)
{
  int flag_used[2], ngrf, error, iptr, igrf;
  double *coor;
  PropDef *propdef;

  /* Initializations */

  error = 1;
  ngrf = 0;
  coor = nullptr;
  propdef = nullptr;

  /**********************/
  /* Preliminary checks */
  /**********************/

  /* Input Db */

  if (db == nullptr)
  {
    messerr("The Db is not defined");
    goto label_end;
  }
  if (!db->isVariableNumberComparedTo(1)) goto label_end;

  /* Rule */

  if (rule == nullptr)
  {
    messerr("The Rule is not defined");
    goto label_end;
  }
  ngrf = rule->getGRFNumber();
  for (igrf = 0; igrf < 2; igrf++)
    flag_used[igrf] = rule->isYUsed(igrf);

  /*******************/
  /* Core allocation */
  /*******************/

  coor = db_sample_alloc(db, ELoc::X);
  if (coor == nullptr) goto label_end;

  propdef = proportion_manage(1, 1, flag_stat, ngrf, 0, nfacies, 0, db, dbprop,
                              props, propdef);
  if (propdef == nullptr) goto label_end;

  /* General setting for lithotype */

  rule->particularities(db, dbprop, model, 1, flag_stat);
  proportion_rule_process(propdef, EProcessOper::COPY);

  /**********************/
  /* Add the attributes */
  /**********************/

  /* Lower bound at input data points */
  if (db_locator_attribute_add(db, ELoc::L, ngrf, 0, 0., &iptr)) goto label_end;

  /* Upper bound at input data points */
  if (db_locator_attribute_add(db, ELoc::U, ngrf, 0, 0., &iptr)) goto label_end;

  /* Calculate the thresholds and store them in the Db file */

  for (igrf = 0; igrf < ngrf; igrf++)
  {
    if (!flag_used[igrf]) continue;
    if (rule->evaluateBounds(propdef, db, db, 0, igrf, 0, 0)) goto label_end;
  }

  /* Set the error return flag */

  error = 0;

  label_end: propdef = proportion_manage(-1, 1, flag_stat, ngrf, 0, nfacies, 0,
                                         db, dbprop, props, propdef);
  coor = db_sample_free(coor);
  return (error);
}

/****************************************************************************/
/*!
 **  Apply the Rule transformation to derive the bounds variables for a Db
 **
 ** \return  Error return code
 **
 ** \param[in]  db        Db structure
 ** \param[in]  ruleprop  RuleProp structure
 ** \param[in]  model     First Model structure (only for SHIFT)
 ** \param[in]  namconv   Naming convention
 **
 *****************************************************************************/
int _db_bounds(Db *db,
               const RuleProp *ruleprop,
               Model *model,
               NamingConvention namconv)
{
  if (db == nullptr)
  {
    messerr("The Db is not defined");
    return 1;
  }
  if (ruleprop == nullptr)
  {
    messerr("RuleProp must be defined");
    return 1;
  }
  int flag_stat = ruleprop->isFlagStat();
  const Rule *rule = ruleprop->getRule();
  const VectorDouble &propcst = ruleprop->getPropCst();
  const Db *dbprop = ruleprop->getDbprop();

  int error = 1;
  int iptrl, iptru;
  PropDef *propdef = nullptr;

  VectorInt flagUsed = rule->whichGRFUsed();
  int nfacies = rule->getFaciesNumber();
  int ngrf = rule->getGRFNumber();

  /* Input Db */

  int nvar = db->getVariableNumber();
  if (!db->isVariableNumberComparedTo(1)) goto label_end;

  /* Model (for SHIFT case) */

  if (rule->checkModel(model, nvar)) goto label_end;

  /*******************/
  /* Core allocation */
  /*******************/

  propdef = proportion_manage(1, 1, flag_stat, ngrf, 0, nfacies, 0, db, dbprop,
                              propcst, propdef);
  if (propdef == nullptr) goto label_end;

  /* General setting for lithotype */

  if (rule->particularities(db, dbprop, model, 1, flag_stat)) goto label_end;
  proportion_rule_process(propdef, EProcessOper::COPY);

  /**********************/
  /* Add the attributes */
  /**********************/

  /* Lower bound at input data points */
  if (db_locator_attribute_add(db, ELoc::L, ngrf, 0, 0., &iptrl))
    goto label_end;

  /* Upper bound at input data points */
  if (db_locator_attribute_add(db, ELoc::U, ngrf, 0, 0., &iptru))
    goto label_end;

  /* Calculate the thresholds and store them in the Db file */

  for (int igrf = 0; igrf < ngrf; igrf++)
  {
    if (!flagUsed[igrf]) continue;
    if (rule->evaluateBounds(propdef, db, db, 0, igrf, 0, 0)) goto label_end;
  }

  // Naming convention

  namconv.setLocatorOutType(ELoc::L);
  namconv.setNamesAndLocators(nullptr, VectorInt(), db, iptrl, "Lower", ngrf);
  namconv.setLocatorOutType(ELoc::U);
  namconv.setNamesAndLocators(nullptr, VectorInt(), db, iptru, "Upper", ngrf);
  error = 0;

  label_end: propdef = proportion_manage(-1, 1, flag_stat, ngrf, 0, nfacies, 0,
                                         db, dbprop, propcst, propdef);
  return (error);
}

/****************************************************************************/
/*!
 **  Set memory proportion so as to provoke the update at first usage
 **
 ** \param[in]  propdef     Pointer to Propdef structure
 **
 ****************************************************************************/
GSTLEARN_EXPORT void propdef_reset(PropDef *propdef)
{
  if (propdef == nullptr) return;
  if (propdef->propmem.empty()) return;

  for (int ifac = 0; ifac < (int) propdef->propmem.size(); ifac++)
    propdef->propmem[ifac] = -1;
}

/****************************************************************************/
/*!
 **  Allocate or deallocate a proportion array
 **
 ** \return  Pointer on the returned PropDef structure
 **
 ** \param[in]  mode        1 for allocation; -1 for deallocation
 ** \param[in]  flag_facies 1 if Gibbs is used for facies
 ** \param[in]  flag_stat   1 if the proportions are stationary
 ** \param[in]  ngrf1       Number of GRFs for the first PGS
 ** \param[in]  ngrf2       Number of GRFs for the second PGS
 ** \param[in]  nfac1       Number of facies for the first PGS
 ** \param[in]  nfac2       Number of facies for the second PGS
 ** \param[in]  db          Db structure containing the data
 ** \param[in]  dbprop      Db structure containing the proportions
 **                         (only used in the non-stationary case)
 ** \param[in]  propcst     Constant set of proportions (used if flag_stat)
 ** \param[in]  proploc     PropDef structure (used for mode<0)
 **
 ****************************************************************************/
GSTLEARN_EXPORT PropDef* proportion_manage(int mode,
                                           int flag_facies,
                                           int flag_stat,
                                           int ngrf1,
                                           int ngrf2,
                                           int nfac1,
                                           int nfac2,
                                           Db *db,
                                           const Db *dbprop,
                                           const VectorDouble &propcst,
                                           PropDef *proploc)
{
  int ifac, error, nfacprod;
  const Db *db_loc;
  PropDef *propdef;

  /* Initializations */

  error = 1;
  nfacprod = nfac1;
  if (nfac2 > 0) nfacprod *= nfac2;

  /* Dispatch */

  if (mode > 0)
  {
    propdef = new PropDef;
    propdef->case_facies = flag_facies;
    propdef->case_stat = flag_stat;
    propdef->case_prop_interp = (dbprop != nullptr && is_grid(dbprop));
    propdef->ngrf[0] = ngrf1;
    propdef->ngrf[1] = ngrf2;
    propdef->nfac[0] = nfac1;
    propdef->nfac[1] = nfac2;
    propdef->nfaccur = nfac1;
    propdef->nfacprod = nfacprod;
    propdef->mode = EProcessOper::UNDEFINED;
    if (propdef->nfaccur <= 0)
    {
      messerr(" The number of facies may not be zero");
      goto label_end;
    }
    propdef->propfix.resize(nfacprod, 0.);
    propdef->propwrk.resize(nfacprod, 0.);
    propdef->proploc.resize(nfacprod, 0.);
    propdef->propmem.resize(nfacprod, 0.);

    if (flag_facies)
    {

      // Case of facies: Use of the proportions

      if (!flag_stat)
      {
        // Non-stationary case

        db_loc = (propdef->case_prop_interp) ? dbprop :
                                               db;
        if (db_loc == nullptr)
        {
          messerr("You have requested Non-stationary proportions");
          messerr("No file is provided containing Proportion variables");
          messerr("Please provide variables with 'proportion' locators");
          messerr("either in the input 'Db' or in 'dbprop'");
          goto label_end;
        }
        if (db_loc->getProportionNumber() != nfacprod)
        {
          messerr(
              "In the non-stationary case, the number of proportion variables (%d)",
              db_loc->getProportionNumber());
          messerr(
              "must be equal to the number of facies (%d) in the Lithotype Rule",
              nfacprod);
          goto label_end;
        }
        propdef->dbprop = db_loc;
        propdef->coor.resize(db_loc->getNDim());
      }
      else
      {

        // Stationary case

        double pref = 1. / (double) nfacprod;
        for (ifac = 0; ifac < nfacprod; ifac++)
        {
          propdef->propfix[ifac] = (propcst.empty()) ? pref :
                                                       propcst[ifac];
          propdef->propwrk[ifac] = (propcst.empty()) ? pref :
                                                       propcst[ifac];
          propdef->proploc[ifac] = (propcst.empty()) ? pref :
                                                       propcst[ifac];
          propdef->propmem[ifac] = (propcst.empty()) ? pref :
                                                       propcst[ifac];
        }
      }

      /* Set memory proportion so as to provoke the update at first usage */

      propdef_reset(propdef);
    }
  }
  else
  {
    propdef = proploc;
    if (propdef == nullptr) return (propdef);

    /* Deallocation */

    propdef->nfaccur = 0;
    propdef->nfacprod = 0;
    propdef->dbprop = nullptr;

    delete propdef;
    propdef = nullptr;
  }

  /* Set the error return code */

  error = 0;

  label_end: if (error)
  {
    if (propdef != nullptr) delete propdef;
    propdef = nullptr;
  }
  return (propdef);
}

/****************************************************************************/
/*!
 **  Calculate all the thresholds at each sample of a Db
 **
 ** \return  Error return code
 **
 ** \param[in]  db        Db structure
 ** \param[in]  ruleprop  RuleProp structure
 ** \param[in]  model     First Model structure (only for SHIFT)
 ** \param[in]  namconv   Naming Convention
 **
 *****************************************************************************/
int _db_threshold(Db *db,
                  const RuleProp *ruleprop,
                  Model *model,
                  NamingConvention namconv)
{
  if (db == nullptr)
  {
    messerr("The Db is not defined");
    return 1;
  }
  if (model == nullptr)
  {
    messerr("The Model is not defined");
    return 1;
  }
  if (ruleprop == nullptr)
  {
    messerr("RuleProp must be defined");
    return 1;
  }
  int flag_stat = ruleprop->isFlagStat();
  const Rule *rule = ruleprop->getRule();
  if (rule->getModeRule() != ERule::STD)
  {
    messerr("This function is only programmed for standard rule");
    return 1;
  }
  const VectorDouble &propcst = ruleprop->getPropCst();
  const Db *dbprop = ruleprop->getDbprop();

  int rank, iptr;
  double t1min, t1max, t2min, t2max;

  /* Initializations */

  int error = 1;
  int ngrf = 0;
  int nfacies = 0;
  PropDef *propdef = nullptr;

  /**********************/
  /* Preliminary checks */
  /**********************/

  ngrf = rule->getGRFNumber();
  if (rule->checkModel(model)) return 1;

  /*******************/
  /* Core allocation */
  /*******************/

  nfacies = rule->getFaciesNumber();
  propdef = proportion_manage(1, 1, flag_stat, ngrf, 0, nfacies, 0, db, dbprop,
                              propcst, propdef);
  if (propdef == nullptr) goto label_end;

  if (rule->particularities(db, dbprop, model, 1, flag_stat)) goto label_end;
  proportion_rule_process(propdef, EProcessOper::COPY);

  /**********************/
  /* Add the attributes */
  /**********************/

  iptr = db->addFields(2 * ngrf * nfacies, 0.);
  if (iptr < 0) goto label_end;

  /* Calculate the thresholds and store them in the Db file */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    rank = 0;
    for (int ifac = 0; ifac < nfacies; ifac++)
    {
      if (rule_thresh_define(propdef, db, rule, ifac + 1, iech, 0, 0, 0, &t1min,
                             &t1max, &t2min, &t2max)) goto label_end;
      db->setArray(iech, iptr + rank, t1min);
      rank++;
      db->setArray(iech, iptr + rank, t1max);
      rank++;
      if (ngrf == 1) continue;
      db->setArray(iech, iptr + rank, t2min);
      rank++;
      db->setArray(iech, iptr + rank, t2max);
      rank++;
    }
  }

  // Naming convention

  rank = 0;
  for (int ifac = 0; ifac < nfacies; ifac++)
  {
    namconv.setNamesAndLocators(
        db, iptr + rank,
        concatenateStrings("Thresh-F", toString(ifac + 1), "-Y1-Low"));
    rank++;
    namconv.setNamesAndLocators(
        db, iptr + rank,
        concatenateStrings("Thresh-F", toString(ifac + 1), "-Y1-Up"));
    rank++;
    if (ngrf == 1) continue;
    namconv.setNamesAndLocators(
        db, iptr + rank,
        concatenateStrings("Thresh-F", toString(ifac + 1), "-Y2-Low"));
    rank++;
    namconv.setNamesAndLocators(
        db, iptr + rank,
        concatenateStrings("Thresh-F", toString(ifac + 1), "-Y2-Up"));
    rank++;
  }
  error = 0;

  label_end: propdef = proportion_manage(-1, 1, flag_stat, ngrf, 0, nfacies, 0,
                                         db, dbprop, propcst, propdef);
  return (error);
}

/****************************************************************************/
/*!
 **  Combine two basic models into a bivariate model (residuals model)
 **
 ** \return  The newly Model structure
 **
 ** \param[in]  model1      First input Model
 ** \param[in]  model2      Second input Model
 ** \param[in]  rule        Rule
 **
 ** \remarks: The drift is not copied into the new model
 **
 *****************************************************************************/
GSTLEARN_EXPORT Model* model_rule_combine(const Model *model1,
                                          const Model *model2,
                                          const Rule *rule)
{
  Model *new_model;
  int ngrf;
  double rho;

  /* Initializations */

  ngrf = 0;
  new_model = nullptr;

  /* Preliminary checks */

  if (rule == nullptr)
  {
    messerr("This function requires a valid rule.");
    return (new_model);
  }
  if (model1 == nullptr)
  {
    messerr("This function requires the first model to be defined");
    return (new_model);
  }
  ngrf = rule->getGRFNumber();

  /* Case of a bivariate input model or monogaussian: simply duplicate */

  if (model1->getVariableNumber() == 2 || ngrf == 1)
  {
    new_model = model_duplicate(model1, 0., 0);
    return (new_model);
  }

  /* If model2 is not defined, consider model1 */

  if (model2 == nullptr)
  {
    if (rule->getModeRule() == ERule::SHIFT)
    {
      new_model = model_duplicate(model1, 0., 0);
      return (new_model);
    }

    model2 = model1;
  }

  /* Subsequent checks */

  if (model1->getVariableNumber() != 1 || model2->getVariableNumber() != 1)
  {
    messerr("This function can only combine monovariate models");
    return (new_model);
  }
  if (model1->getDimensionNumber() != model2->getDimensionNumber())
  {
    messerr("The two models to be combined must share the space dimension");
    return (new_model);
  }
  if (model1->isFlagLinked() || model2->isFlagLinked())
  {
    messerr("This function cannot combine models with linked drifts");
    return (new_model);
  }

  /* Calculate the correlation coefficient */

  rho = 0.;
  if (rule->getModeRule() == ERule::STD) rho = rule->getRho();

  new_model = model_combine(model1, model2, rho);
  return (new_model);
}
