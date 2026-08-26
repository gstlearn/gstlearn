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
#include "Basic/String.hpp"
#include "Db/Db.hpp"
#include "LithoRule/PropDef.hpp"
#include "LithoRule/Rule.hpp"
#include "LithoRule/RuleProp.hpp"
#include "LithoRule/RuleShadow.hpp"
#include "Model/Model.hpp"
#include "Variogram/Vario.hpp"
#include "geoslib_old_f.h"

namespace gstlrn
{
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
   ** \param[in]  flag_stat true for stationary; false otherwise
   ** \param[in]  nfacies   Number of facies
   **
   ** \remark The input variable must be locatorized as Z or ELoc::SIMU
   ** \remark It will be changed in this function to locator ELoc::SIMU
   **
   *****************************************************************************/
  Id db_rule_shadow(
    Db* db,
    Db* dbprop,
    RuleShadow* rule,
    Model* model,
    const VectorDouble& props,
    bool flag_stat,
    Id nfacies)
  {
    Id iptr, error, nbsimu, igrf, ngrf;
    PropDef propdef;

    /* Initializations */

    error = 1;
    nbsimu = 1;
    iptr = -1;

    /* Preliminary checks */

    if (rule == nullptr)
    {
      messerr("The Rule is not defined");
      return 1;
    }
    ngrf = rule->getNGRF();
    VectorBool flagUsed = rule->whichGRFUsed();
    if (propdef.define(
          true, flag_stat, {ngrf, 0}, {nfacies, 0}, db, dbprop, props))
      goto label_end;

    /* General setting for lithotype */

    rule->particularities(db, dbprop, model, 1, flag_stat);
    propdef.defineRuleMethod(EProcessOper::COPY);

    /**********************/
    /* Add the attributes */
    /**********************/

    /* Storage of the simulations in the output file */
    iptr = db->addColumnsByConstant(nbsimu, 0.);
    if (iptr < 0) goto label_end;
    db->setLocatorsByUID(nbsimu, iptr, ELoc::FACIES, 0);

    /* Identify the Non conditional simulations at target points */
    for (igrf = 0; igrf < 2; igrf++)
    {
      if (!flagUsed[igrf]) continue;
      iptr = db->getUIDByLocator(ELoc::SIMU, igrf);
      if (iptr < 0)
      {
        iptr = db->getUIDByLocator(ELoc::Z, igrf);
        if (iptr < 0)
        {
          messerr(
            "The variable containing the simulation of the GRF %d is missing "
            "in the Db",
            igrf + 1);
          goto label_end;
        }
        db->setLocatorByUID(iptr, ELoc::SIMU, igrf);
      }
    }

    /* Combine the conditional simulation for each GRF */

    for (Id isimu = 0; isimu < nbsimu; isimu++)
      if (rule->gaus2facResult(propdef, db, flagUsed, 0, isimu, nbsimu))
        goto label_end;

    /* Set the error return flag */

    error = 0;

  label_end:
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
  Id _db_rule(
    Db* db,
    const RuleProp* ruleprop,
    Model* model,
    const NamingConvention& namconv)
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
    bool flag_stat = ruleprop->isFlagStat();
    const Rule* rule = ruleprop->getRule();
    if (rule == nullptr)
    {
      messerr("The Rule is not defined in the RuleProp");
      return 1;
    }
    const VectorDouble& propcst = ruleprop->getPropCst();
    const Db* dbprop = ruleprop->getDbprop();

    Id error = 1;
    Id iptr = -1;
    PropDef propdef;
    Id ngrf = rule->getNGRF();
    VectorBool flagUsed = rule->whichGRFUsed();
    Id nfacies = rule->getNFacies();
    bool flagReturn = false;

    /* Preliminary checks */

    Id nbsimu = db->getNLoc(ELoc::SIMU);
    Id nvar = db->getNLoc(ELoc::Z);
    if (nbsimu != ngrf && nvar != ngrf)
    {
      messerr("The Rule specifies the use of %d underlying GRF(s)", ngrf);
      messerr(
        "The input 'db' should have one variable per GRF with locator 'SIMU' "
        "or 'Z'");
      goto label_end;
    }

    if (propdef.define(
          true, flag_stat, {ngrf, 0}, {nfacies, 0}, db, dbprop, propcst))
      goto label_end;
    if (rule->particularities(db, dbprop, model, 1, flag_stat)) goto label_end;
    propdef.defineRuleMethod(EProcessOper::COPY);

    /**********************/
    /* Add the attributes */
    /**********************/

    /* Storage of the simulations in the output file */
    iptr = db->addColumnsByConstant(1, 0., "Facies", ELoc::FACIES);
    if (iptr < 0) goto label_end;

    /* Identify the Non conditional simulations at target points */

    if (db->getNLoc(ELoc::SIMU) != ngrf)
    {
      db->switchLocator(ELoc::Z, ELoc::SIMU);
      flagReturn = true;
    }

    /* Translate Gaussian into Facies */

    if (rule->gaus2facResult(propdef, db, flagUsed, 0, 0, 1)) goto label_end;

    // Returning to the initial locators (if the initial variable
    // had a ELoc::Z locator which has been temporarily modified into ELoc::SIMU)

    if (flagReturn) db->switchLocator(ELoc::SIMU, ELoc::Z);

    // Naming convention

    namconv.setOutput(VectorString(), 1, db, iptr);
    error = 0;

  label_end:
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
   ** \param[in]  flag_stat true for stationary; false otherwise
   ** \param[in]  nfacies   Number of facies
   **
   *****************************************************************************/
  Id db_bounds_shadow(
    Db* db,
    Db* dbprop,
    RuleShadow* rule,
    Model* model,
    const VectorDouble& props,
    bool flag_stat,
    Id nfacies)
  {
    Id iptr;

    /* Input Db */

    if (db == nullptr)
    {
      messerr("The Db is not defined");
      return 1;
    }
    if (!db->isNVarComparedTo(1)) return 1;

    /* Rule */

    if (rule == nullptr)
    {
      messerr("The Rule is not defined");
      return 1;
    }
    Id ngrf = rule->getNGRF();
    VectorBool flagUsed = rule->whichGRFUsed();

    /*******************/
    /* Core allocation */
    /*******************/

    Id ndim = db->getNDim();
    PropDef propdef;
    VectorDouble coor(ndim);

    if (propdef.define(
          true, flag_stat, {ngrf, 0}, {nfacies, 0}, db, dbprop, props))
      return 1;

    /* General setting for lithotype */

    if (rule->particularities(db, dbprop, model, 1, flag_stat)) return 1;
    propdef.defineRuleMethod(EProcessOper::COPY);

    /**********************/
    /* Add the attributes */
    /**********************/

    /* Lower bound at input data points */
    if (db_locator_attribute_add(db, ELoc::L, ngrf, 0, 0., &iptr)) return 1;

    /* Upper bound at input data points */
    if (db_locator_attribute_add(db, ELoc::U, ngrf, 0, 0., &iptr)) return 1;

    /* Calculate the thresholds and store them in the Db file */

    for (Id igrf = 0; igrf < ngrf; igrf++)
    {
      if (!flagUsed[igrf]) continue;
      if (rule->evaluateBounds(propdef, db, db, 0, igrf, 0, 0)) return 1;
    }

    return 0;
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
  Id _db_bounds(
    Db* db,
    const RuleProp* ruleprop,
    Model* model,
    const NamingConvention& namconv)
  {
    NamingConvention nc(namconv);

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
    bool flag_stat = ruleprop->isFlagStat();
    const Rule* rule = ruleprop->getRule();
    if (rule == nullptr)
    {
      messerr("The Rule is not defined in the RuleProp");
      return 1;
    }
    const VectorDouble& propcst = ruleprop->getPropCst();
    const Db* dbprop = ruleprop->getDbprop();

    Id error = 1;
    Id iptrl, iptru;
    PropDef propdef;

    VectorBool flagUsed = rule->whichGRFUsed();
    Id nfacies = rule->getNFacies();
    Id ngrf = rule->getNGRF();

    /* Input Db */

    Id nvar = db->getNLoc(ELoc::Z);
    if (!db->isNVarComparedTo(1)) goto label_end;

    /* Model (for SHIFT case) */

    if (rule->checkModel(model, nvar)) goto label_end;

    /*******************/
    /* Core allocation */
    /*******************/

    if (propdef.define(
          true, flag_stat, {ngrf, 0}, {nfacies, 0}, db, dbprop, propcst))
      goto label_end;

    /* General setting for lithotype */

    if (rule->particularities(db, dbprop, model, 1, flag_stat)) goto label_end;
    propdef.defineRuleMethod(EProcessOper::COPY);

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

    for (Id igrf = 0; igrf < ngrf; igrf++)
    {
      if (!flagUsed[igrf]) continue;
      if (rule->evaluateBounds(propdef, db, db, 0, igrf, 0, 0)) goto label_end;
    }

    // Naming convention

    nc.setLocatorOutType(ELoc::L);
    nc.setOutput(VectorString(), ngrf, db, iptrl, "Lower", ngrf);
    nc.setLocatorOutType(ELoc::U);
    nc.setOutput(VectorString(), ngrf, db, iptru, "Upper", ngrf);
    error = 0;

  label_end:
    return (error);
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
  Id _db_threshold(
    Db* db,
    const RuleProp* ruleprop,
    Model* model,
    const NamingConvention& namconv)
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
    bool flag_stat = ruleprop->isFlagStat();
    const Rule* rule = ruleprop->getRule();
    if (rule == nullptr)
    {
      messerr("The Rule is not defined in the RuleProp");
      return 1;
    }
    if (rule->getModeRule() != ERule::STD)
    {
      messerr("This function is only programmed for standard rule");
      return 1;
    }
    const VectorDouble& propcst = ruleprop->getPropCst();
    const Db* dbprop = ruleprop->getDbprop();

    Id rank, iptr;
    double t1min, t1max, t2min, t2max;

    /* Initializations */

    Id error = 1;
    Id ngrf = 0;
    Id nfacies = 0;
    PropDef propdef;

    /**********************/
    /* Preliminary checks */
    /**********************/

    ngrf = rule->getNGRF();
    if (rule->checkModel(model)) return 1;

    /*******************/
    /* Core allocation */
    /*******************/

    nfacies = rule->getNFacies();
    if (propdef.define(
          true, flag_stat, {ngrf, 0}, {nfacies, 0}, db, dbprop, propcst))
      goto label_end;
    if (rule->particularities(db, dbprop, model, 1, flag_stat)) goto label_end;
    propdef.defineRuleMethod(EProcessOper::COPY);

    /**********************/
    /* Add the attributes */
    /**********************/

    iptr = db->addColumnsByConstant(2 * ngrf * nfacies, 0.);
    if (iptr < 0) goto label_end;

    /* Calculate the thresholds and store them in the Db file */

    for (Id iech = 0; iech < db->getNSample(); iech++)
    {
      if (!db->isActive(iech)) continue;
      rank = 0;
      for (Id ifac = 0; ifac < nfacies; ifac++)
      {
        if (propdef.ruleThreshDefine(
              db, rule, ifac + 1, iech, 0, 0, 0, &t1min, &t1max, &t2min,
              &t2max))
          goto label_end;
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
    for (Id ifac = 0; ifac < nfacies; ifac++)
    {
      namconv.setOutput(
        VectorString(), 0, db, iptr + rank,
        concatenateStrings("Thresh-F", toStr(ifac + 1), "-Y1-Low"));
      rank++;
      namconv.setOutput(
        VectorString(), 0, db, iptr + rank,
        concatenateStrings("Thresh-F", toStr(ifac + 1), "-Y1-Up"));
      rank++;
      if (ngrf == 1) continue;
      namconv.setOutput(
        VectorString(), 0, db, iptr + rank,
        concatenateStrings("Thresh-F", toStr(ifac + 1), "-Y2-Low"));
      rank++;
      namconv.setOutput(
        VectorString(), 0, db, iptr + rank,
        concatenateStrings("Thresh-F", toStr(ifac + 1), "-Y2-Up"));
      rank++;
    }
    error = 0;

  label_end:

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
  Model* model_rule_combine(
    const Model* model1,
    const Model* model2,
    const Rule* rule)
  {
    Model* new_model;
    Id ngrf;
    double rho;

    /* Initializations */

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
    ngrf = rule->getNGRF();

    /* Case of a bivariate input model or monogaussian: simply duplicate */

    if (model1->getNVar() == 2 || ngrf == 1)
    {
      new_model = model1->clone();
      return (new_model);
    }

    /* If model2 is not defined, consider model1 */

    if (model2 == nullptr)
    {
      if (rule->getModeRule() == ERule::SHIFT)
      {
        new_model = model1->clone();
        return (new_model);
      }

      model2 = model1;
    }

    /* Subsequent checks */

    if (model1->getNVar() != 1 || model2->getNVar() != 1)
    {
      messerr("This function can only combine monovariate models");
      return (new_model);
    }
    if (model1->getNDim() != model2->getNDim())
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
} // namespace gstlrn
