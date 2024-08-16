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
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "geoslib_f_private.h"

#include "Enum/EProcessOper.hpp"
#include "Enum/ERule.hpp"

#include "Calculators/CalcMigrate.hpp"
#include "Gibbs/GibbsUMultiMono.hpp"
#include "Gibbs/GibbsUPropMono.hpp"
#include "Gibbs/GibbsMMulti.hpp"
#include "Gibbs/GibbsFactory.hpp"
#include "Covariances/CovAniso.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Basic/OptDbg.hpp"
#include "LithoRule/PropDef.hpp"
#include "LithoRule/Rule.hpp"
#include "LithoRule/RuleShift.hpp"
#include "LithoRule/RuleShadow.hpp"
#include "LithoRule/RuleProp.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Simulation/SimuBoolean.hpp"
#include "Simulation/SimuSpherical.hpp"
#include "Simulation/SimuSphericalParam.hpp"
#include "Simulation/SimuRefineParam.hpp"
#include "Simulation/SimuRefine.hpp"
#include "Simulation/CalcSimuEden.hpp"
#include "Simulation/CalcSimuFFT.hpp"

#include <math.h>
#include <string.h>

/*! \cond */
#define DATA   0
#define RESULT 1

// TODO : transform this to enum
#define TYPE_GAUS   0
#define TYPE_FACIES 1
#define TYPE_PROP   2

/*! \endcond */

static double GIBBS_RHO, GIBBS_SQR;
static Modif_Categorical ModCat = { 0, { 0, 0 }, NULL, NULL };

/****************************************************************************/
/*!
 **  Initialize the global values
 **
 *****************************************************************************/
static void st_simulation_environment(void)
{
  GIBBS_RHO = 0.;
  GIBBS_SQR = 0.;
}

/****************************************************************************/
/*!
 **  Give the rank of a proportion for a given GRF and PGS
 **
 ** \return  Returned rank
 **
 ** \param[in]  propdef    PropDef structure
 ** \param[in]  ifac       Rank of the facies
 ** \param[in]  ipgs       Rank of the GS
 **
 *****************************************************************************/
static int st_facies(PropDef *propdef, int ipgs, int ifac)
{
  if (ipgs <= 0) return (ifac);
  return (propdef->nfac[0] + ifac);
}

/****************************************************************************/
/*!
 **  Transformation function for the Modification categorical case
 **
 ** \param[in]  db        Db structure
 ** \param[in]  verbose   1 for the verbose flag
 ** \param[in]  isimu     Rank of the current simulation
 ** \param[in]  nbsimu    Number of simulations
 **
 *****************************************************************************/
void simu_func_categorical_transf(Db *db, int verbose, int isimu, int nbsimu)
{
  const Rule *rule = ModCat.rule;

  rule->gaus2facResult(ModCat.propdef, db, ModCat.flag_used, ModCat.ipgs, isimu,
                       nbsimu);

  /* Optional printout */

  if (verbose)
    message("Simulation Categorical Transformation (%d/%d)\n", isimu + 1,
            nbsimu);
}

/****************************************************************************/
/*!
 **  Updating function for the Modification continuous case
 **
 ** \param[in]  db        Db structure
 ** \param[in]  verbose   1 for the verbose flag
 ** \param[in]  isimu     Rank of the current simulation
 ** \param[in]  nbsimu    Number of simulations (stored)
 **
 *****************************************************************************/
void simu_func_continuous_update(Db *db, int verbose, int isimu, int nbsimu)
{
  int iptr_simu;
  double simval;

  /* Preliminary checks */

  check_mandatory_attribute("simu_func_continuous_update", db, ELoc::SIMU);
  check_mandatory_attribute("simu_func_continuous_update", db, ELoc::Z);
  iptr_simu = Db::getSimRank(isimu, 0, 0, nbsimu, 1);

  /* Loop on the grid cells */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    simval = db->getFromLocator(ELoc::SIMU, iech, iptr_simu);
    db->updLocVariable(ELoc::Z,iech, 0, EOperator::ADD, simval);
    db->updLocVariable(ELoc::Z,iech, 1, EOperator::ADD, simval * simval);
  }

  /* Optional printout */

  if (verbose)
    message("Simulation Continuous Update (%d/%d)\n", isimu + 1, nbsimu);
}

/****************************************************************************/
/*!
 **  Updating function for the Modification categorical case
 **
 ** \param[in]  db        Db structure
 ** \param[in]  verbose   1 for the verbose flag
 ** \param[in]  isimu     Rank of the current simulation
 ** \param[in]  nbsimu    Number of simulations (stored)
 **
 *****************************************************************************/
void simu_func_categorical_update(Db *db, int verbose, int isimu, int nbsimu)
{
  int iptr_simu, facies, rank, ipgs;
  double prop;

  /* Preliminary checks */

  ipgs = ModCat.ipgs;
  check_mandatory_attribute("simu_func_categorical_update", db, ELoc::FACIES);
  check_mandatory_attribute("simu_func_categorical_update", db, ELoc::P);
  iptr_simu = Db::getSimRank(isimu, 0, ipgs, nbsimu, 1);

  /* Loop on the grid cells */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    facies = (int) db->getFromLocator(ELoc::FACIES, iech, iptr_simu) - 1;
    rank = st_facies(ModCat.propdef, ipgs, facies);
    prop = db->getLocVariable(ELoc::P,iech, rank) + 1.;
    db->setLocVariable(ELoc::P,iech, rank, prop);
  }

  /* Optional printout */

  if (verbose)
    message("Simulation Categorical Update (%d/%d)\n", isimu + 1, nbsimu);
}

/****************************************************************************/
/*!
 **  Scaling function for the Modification continuous case
 **
 ** \param[in]  db        Db structure
 ** \param[in]  verbose   1 for the verbose flag
 ** \param[in]  nbsimu    Number of simulations
 **
 *****************************************************************************/
void simu_func_continuous_scale(Db *db, int verbose, int nbsimu)
{
  double mean, stdv;

  /* Preliminary checks */

  check_mandatory_attribute("simu_func_continuous_scale", db, ELoc::Z);

  /* Loop on the grid cells */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    mean = db->getZVariable(iech, 0) / nbsimu;
    db->setLocVariable(ELoc::Z,iech, 0, mean);
    stdv = db->getZVariable(iech, 1) / nbsimu - mean * mean;
    stdv = (stdv > 0) ? sqrt(stdv) :
                        0.;
    db->setLocVariable(ELoc::Z,iech, 1, stdv);
  }

  /* Optional printout */

  if (verbose) message("Simulation Continuous Scaling (%d)\n", nbsimu);
}

/****************************************************************************/
/*!
 **  Scaling function for the Modification categorical case
 **
 ** \param[in]  db        Db structure
 ** \param[in]  verbose   1 for the verbose flag
 ** \param[in]  nbsimu    Number of simulations
 **
 *****************************************************************************/
void simu_func_categorical_scale(Db *db, int verbose, int nbsimu)
{
  int rank, nfacies, ipgs;
  double prop;
  PropDef *propdef;

  /* Preliminary checks */

  propdef = ModCat.propdef;
  ipgs = ModCat.ipgs;
  nfacies = propdef->nfac[ipgs];
  check_mandatory_attribute("simu_func_categorical_scale", db, ELoc::P);

  /* Loop on the grid cells */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    for (int ifac = 0; ifac < nfacies; ifac++)
    {
      rank = st_facies(propdef, ipgs, ifac);
      prop = db->getLocVariable(ELoc::P,iech, rank) / (double) nbsimu;
      db->setLocVariable(ELoc::P,iech, rank, prop);
    }
  }

  /* Optional printout */

  if (verbose) message("Simulation Categorical Scaling (%d)\n", nbsimu);
}

/****************************************************************************/
/*!
 **  Check for the presence of mandatory attributes
 **
 ** \param[in]  method  Name of the method
 ** \param[in]  db      Db structure
 ** \param[in]  locatorType  Mandatory attribute type
 **
 *****************************************************************************/
void check_mandatory_attribute(const char *method,
                               Db *db,
                               const ELoc& locatorType)
{
  if (get_LOCATOR_NITEM(db,locatorType) <= 0)
    messageAbort("%s : Attributes %d are mandatory",method,locatorType.getValue());
}

/****************************************************************************/
/*!
 **  Check if the field must be kept
 **
 ** \return  1 if the field must be kept; 0 otherwise
 **
 ** \param[in]  flag_gaus  1 gaussian results; otherwise facies
 ** \param[in]  flag_prop  1 for facies proportion
 ** \param[in]  file       DATA or RESULT
 ** \param[in]  type       0 for gaussian; 1 for facies; 2 for proportion
 **
 *****************************************************************************/
static int st_keep(int flag_gaus, int flag_prop, int file, int type)
{
  int keep;

  keep = 0;

  if (file == DATA)
  {

    /* Input Db */

    goto label_end;
  }
  else
  {

    /* Output Db */

    switch (type)
    {
      case 0: /* Gaussian */
        keep = (flag_gaus && !flag_prop);
        break;

      case 1: /* Facies */
        keep = (!flag_gaus && !flag_prop);
        break;

      case 2: /* Proportion */
        keep = (flag_prop);
        break;

      default:
        break;
    }
  }

  label_end: return (keep);
}

/****************************************************************************/
/*!
 **  Checks the environment for simulations by Turning Bands
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       input Db structure (optional if non conditional)
 ** \param[in]  dbout      output Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  neigh      ANeigh structure (optional if non conditional)
 **
 *****************************************************************************/
static int st_check_simtub_environment(Db *dbin,
                                       Db *dbout,
                                       Model *model,
                                       ANeigh *neigh)
{
  int nvar = 0;
  int nfex = 0;
  bool flag_cond = (dbin != nullptr);
  int ndim = dbout->getNDim();

  /**************************************************************/
  /* Check if the Space dimension is compatible with the method */
  /**************************************************************/

  if (ndim > 3)
  {
    messerr("The Turning Band Method is not a relevant simulation model");
    messerr("for this Space Dimension (%d)", ndim);
    return 1;
  }

  /*********************************/
  /* Compatibility between two Dbs */
  /*********************************/

  if (flag_cond && !dbin->hasSameDimension(dbout)) return 1;

  /**********************/
  /* Checking the model */
  /**********************/

  if (model != nullptr)
  {
    nvar = model->getVariableNumber();
    if (nvar <= 0)
    {
      messerr("The number of variables must be positive = %d",
              model->getVariableNumber());
      return 1;
    }
    if (flag_cond && dbin->getLocNumber(ELoc::Z) != nvar)
    {
      messerr("The number of variables of the Data (%d)",
              dbin->getLocNumber(ELoc::Z));
      messerr("does not match the number of variables of the Model (%d)", nvar);
      return 1;
    }
    if (model->getCovaNumber() <= 0)
    {
      messerr("The number of covariance must be positive");
      return 1;
    }

    if (model->getDimensionNumber() <= 0)
    {
      messerr("The Space Dimension must be positive = %d",
              model->getDimensionNumber());
      return 1;
    }
    if (model->getDimensionNumber() != ndim)
    {
      messerr("The Space Dimension of the Db structure (%d)", ndim);
      messerr("Does not correspond to the Space Dimension of the model (%d)",
              model->getDimensionNumber());
      return 1;
    }

    nfex = model->getExternalDriftNumber();
    if (flag_cond && nfex != 0 && ! dbout->isGrid()
        && dbin->getLocNumber(ELoc::F) != nfex)
    {
      messerr("The Model requires %d external drift(s)", model->getExternalDriftNumber());
      messerr("but the input Db refers to %d external drift variables",
              dbin->getLocNumber(ELoc::F));
      return 1;
    }
    if (nfex != 0 && dbout->getLocNumber(ELoc::F) != nfex)
    {
      messerr("The Model requires %d external drift(s)", model->getExternalDriftNumber());
      messerr("but the output Db refers to %d external drift variables",
              dbout->getLocNumber(ELoc::F));
      return 1;
    }
  }

  /*********************************/
  /* Calculate the field extension */
  /*********************************/

  VectorDouble db_mini(ndim, TEST);
  VectorDouble db_maxi(ndim, TEST);

  dbout->getExtensionInPlace(db_mini, db_maxi, true);

  if (flag_cond)
    dbin->getExtensionInPlace(db_mini, db_maxi, true);

  if (model != nullptr)
    model->setField(VH::extensionDiagonal(db_mini, db_maxi));

  /*****************************/
  /* Checking the Neighborhood */
  /*****************************/

  if (flag_cond && neigh != nullptr)
  {
    if (ndim != (int) neigh->getNDim())
    {
      messerr("The Space Dimension of the Neighborhood (%d)", (int) neigh->getNDim());
      messerr("does not correspond to the Space Dimension of the first Db (%d)",
              ndim);
      return 1;
    }
    if (neigh->getFlagXvalid() && neigh->getType() != ENeigh::MOVING)
    {
      messerr("The Cross-Validation can only be processed with Moving neighborhood");
      return 1;
    }
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Give the rank of a "variable" for a given GRF and PGS
 **
 ** \return  Returned rank
 **
 ** \param[in]  propdef    PropDef structure
 ** \param[in]  ipgs       Rank of the GS
 ** \param[in]  igrf       Rank of the Gaussian
 **
 *****************************************************************************/
int get_rank_from_propdef(PropDef *propdef, int ipgs, int igrf)
{
  if (ipgs <= 0 || propdef == nullptr) return (igrf);
  return (propdef->ngrf[0] + igrf);
}

/****************************************************************************/
/*!
 **  Suppresses the added samples
 **
 ** \param[in]  db      Db structure
 ** \param[in]  nech    initial number of samples
 **
 *****************************************************************************/
static void st_suppress_added_samples(Db *db, int nech)
{
  int iech;

  if (nech <= 0) return;
  for (iech = db->getSampleNumber() - 1; iech >= nech; iech--)
    (void) db->deleteSample(iech);
}

/****************************************************************************/
/*!
 **  Check/Show the data against facies at the closest grid node
 **
 ** \param[in]  dbin       Input Db structure
 ** \param[in]  dbout      Output Db grid structure
 ** \param[in]  flag_check 1 check the consistency between data and grid
 ** \param[in]  flag_show  1 show the data on grid
 ** \param[in]  ipgs       Rank of the PGS
 ** \param[in]  nechin     Initial number of data
 ** \param[in]  nfacies    Number of facies
 ** \param[in]  nbsimu     Number of simulations
 **
 ** \remark Attributes ELoc::FACIES are mandatory
 ** \remark Attributes ELoc::GAUSFAC are mandatory
 **
 *****************************************************************************/
static void st_check_facies_data2grid(Db *dbin,
                                      Db *dbout,
                                      int flag_check,
                                      int flag_show,
                                      int ipgs,
                                      int nechin,
                                      int nfacies,
                                      int nbsimu)
{
  int iech, jech, isimu, facdat, facres, number;
  double *coor;

  /* Initializations */

  if (! dbout->isGrid()) return;
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(dbout);
  check_mandatory_attribute("st_check_facies_data2grid", dbgrid, ELoc::FACIES);
  number = 0;
  coor = nullptr;
  if (flag_check)
    mestitle(1, "Checking facies of data against closest grid node (PGS=%d)",
             ipgs + 1);

  /* Core allocation */

  coor = db_sample_alloc(dbin, ELoc::X);
  if (coor == nullptr) goto label_end;

  /* Loop on the data */

  for (iech = 0; iech < nechin; iech++)
  {
    if (!dbin->isActive(iech)) continue;
    facdat = (int) dbin->getZVariable(iech, 0);
    if (facdat < 1 || facdat > nfacies) continue;
    jech = index_point_to_grid(dbin, iech, 0, dbgrid, coor);
    if (jech < 0) continue;

    for (isimu = 0; isimu < nbsimu; isimu++)
    {
      facres = (int) dbgrid->getSimvar(ELoc::FACIES, jech, isimu, 0, ipgs,
                                      nbsimu, 1);
      if (flag_show)
      {
        if (facdat == facres)
          dbgrid->setSimvar(ELoc::FACIES, jech, isimu, 0, ipgs, nbsimu, 1,
                           -facdat);
        else
          dbgrid->setSimvar(ELoc::FACIES, jech, isimu, 0, ipgs, nbsimu, 1, 0.);
      }

      if (facdat == facres) continue;
      number++;

      /* The data facies is different from the grid facies */

      if (flag_check)
      {
        message("Inconsistency for Simulation (%d) between :\n", isimu + 1);
        message("- Facies (%d) at Data (#%d)\n", facdat, iech + 1);
        message("- Facies (%d) at Grid (#%d)\n", facres, jech + 1);
      }
    }
  }

  label_end: if (flag_check && number <= 0) message("No problem found\n");
  db_sample_free(coor);
}

/****************************************************************************/
/*!
 **  Initialize the Gibbs internal parameters
 **
 ** \param[in]  rho        Correlation between the two underlying GRF
 **
 *****************************************************************************/
static void st_init_gibbs_params(double rho)
{
  GIBBS_RHO = rho;
  GIBBS_SQR = sqrt(1. - rho * rho);
}

/****************************************************************************/
/*!
 **  Perform the conditional or non-conditional Pluri-gaussian
 **  simulations
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure (optional)
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  ruleprop    RuleProp structure
 ** \param[in]  model1      First Model structure
 ** \param[in]  model2      Second Model structure (optional)
 ** \param[in]  neigh       ANeigh structure
 ** \param[in]  nbsimu      Number of simulations
 ** \param[in]  seed        Seed for random number generator
 ** \param[in]  flag_gaus   1 if results must be gaussian; otherwise facies
 ** \param[in]  flag_prop   1 for facies proportion
 ** \param[in]  flag_check  1 if the facies at data must be checked against
 **                         the closest simulated grid node
 ** \param[in]  flag_show   1 if the grid node which coincides with the data
 **                         should be represented with the data facies
 **                         (only if flag_cond && !flag_gaus)
 ** \param[in]  nbtuba      Number of turning bands
 ** \param[in]  gibbs_nburn Number of bootstrap iterations
 ** \param[in]  gibbs_niter Maximum number of iterations
 ** \param[in]  percent     Amount of nugget effect added to too much continous
 **                         model (expressed in percentage of the total variance)
 ** \param[in]  namconv     Naming convention
 **
 ** \remark  When conditional, the unique variable in the input Db structure
 ** \remark  should correspond to the facies index (starting from 1)
 ** \remark  The argument 'dbin' is optional: it must be defined only for
 ** \remark  conditional simulations
 **
 *****************************************************************************/
int simpgs(Db* dbin,
           Db* dbout,
           RuleProp* ruleprop,
           Model* model1,
           Model* model2,
           ANeigh* neigh,
           int nbsimu,
           int seed,
           int flag_gaus,
           int flag_prop,
           int flag_check,
           int flag_show,
           int nbtuba,
           int gibbs_nburn,
           int gibbs_niter,
           double percent,
           const NamingConvention& namconv)
{
  int iptr, icase, nfacies, flag_used[2];
  int iptr_RP, iptr_RF, iptr_DF, iptr_DN, iptr_RN, local_seed;
  Model *models[2];
  PropDef *propdef;
  std::vector<Model*> modvec;

  /* Initializations */

  int error = 1;
  int nechin = 0;
  int ngrf = 0;
  propdef = nullptr;
  models[0] = model1;
  models[1] = model2;
  bool flag_cond = (dbin != nullptr);
  iptr = iptr_RP = iptr_RF = iptr_DF = iptr_DN = iptr_RN = -1;
  nfacies = 0;
  bool verbose = false;

  if (ruleprop == nullptr)
  {
    messerr("RuleProp must be defined");
    return 1;
  }
  int flag_stat = ruleprop->isFlagStat();
  const Rule *rule = ruleprop->getRule();
  const VectorDouble &propcst = ruleprop->getPropCst();
  const Db *dbprop = ruleprop->getDbprop();

  ngrf = rule->getGRFNumber();
  if (rule->particularities(dbout, dbprop, model1, 1, flag_stat))
    goto label_end;
  if (st_check_simtub_environment(dbin, dbout, model1, neigh)) goto label_end;

  /**********************/
  /* Preliminary checks */
  /**********************/

  /* Input Db */
  if (flag_cond)
  {
    nechin = dbin->getSampleNumber();
    if (!dbin->isVariableNumberComparedTo(1)) goto label_end;
  }

  /* Output Db */
  if (flag_prop && flag_gaus)
  {
    messerr(
        "Calculating the facies proportions is incompatible with storing the Gaussian values");
    goto label_end;
  }

  /* Model */
  for (int igrf = 0; igrf < 2; igrf++)
  {
    flag_used[igrf] = rule->isYUsed(igrf);
    if (!flag_used[igrf]) continue;
    if (models[igrf] == nullptr)
    {
      messerr("The Underlying GRF #%d is needed", igrf + 1);
      messerr("No corresponding Model is provided");
      goto label_end;
    }
    if (models[igrf]->getVariableNumber() != 1)
    {
      messerr("The number of variables in the model #%d (%d) should be 1",
              igrf + 1, model1->getVariableNumber());
      goto label_end;
    }
    if (models[igrf]->stabilize(percent, true)) goto label_end;
    if (models[igrf]->standardize(true)) goto label_end;
    modvec.push_back(models[igrf]);
  }

  /* Neighborhood */
  if (flag_cond)
  {
    if (neigh->getType() != ENeigh::UNIQUE && neigh->getType() != ENeigh::BENCH)
    {
      messerr("The only authorized Neighborhoods are UNIQUE or BENCH");
      goto label_end;
    }
  }

  /* Define the environment variables for printout */

  /**********************/
  /* Add the attributes */
  /**********************/

  nfacies = rule->getFaciesNumber();

  /* Storage of the facies proportions */
  if (flag_prop)
  {
    if (db_locator_attribute_add(dbout, ELoc::P, nfacies, 0, 0., &iptr_RP))
      goto label_end;
  }

  /* Storage of the facies simulations in the Output Db */
  if (db_locator_attribute_add(dbout, ELoc::FACIES, nbsimu, 0, 0., &iptr_RF))
    goto label_end;

  /* Storage of the facies simulations in the input file */
  if (flag_cond)
  {
    if (db_locator_attribute_add(dbin, ELoc::FACIES, nbsimu, 0, 0., &iptr_DF))
      goto label_end;
  }

  if (flag_cond)
  {
    /* Gaussian transform of the facies input data */
    if (db_locator_attribute_add(dbin, ELoc::GAUSFAC, ngrf * nbsimu, 0, 0.,
                                 &iptr)) goto label_end;

    /* Non-conditional simulations at data points */
    if (db_locator_attribute_add(dbin, ELoc::SIMU, ngrf * nbsimu, 0, 0.,
                                 &iptr_DN)) goto label_end;
  }

  /* (Non-) Conditional simulations at target points */
  if (db_locator_attribute_add(dbout, ELoc::SIMU, ngrf * nbsimu, 0, 0.,
                               &iptr_RN)) goto label_end;

  if (flag_cond)
  {
    /* Lower bound at input data points */
    if (db_locator_attribute_add(dbin, ELoc::L, ngrf, 0, 0., &iptr))
      goto label_end;

    /* Upper bound at input data points */
    if (db_locator_attribute_add(dbin, ELoc::U, ngrf, 0, 0., &iptr))
      goto label_end;
  }

  propdef = proportion_manage(1, 1, flag_stat, ngrf, 0, nfacies, 0, dbin,
                              dbprop, propcst, propdef);
  if (propdef == nullptr) goto label_end;
  simu_define_func_update(simu_func_categorical_update);
  simu_define_func_scale(simu_func_categorical_scale);
  ModCat.propdef = propdef;
  ModCat.rule = rule;
  ModCat.ipgs = 0;
  ModCat.flag_used[0] = flag_used[0];
  ModCat.flag_used[1] = flag_used[1];

  /****************************************/
  /* Convert facies into gaussian at data */
  /****************************************/

  proportion_rule_process(propdef, EProcessOper::COPY);

  if (flag_cond)
  {
    int npgs = 1;
    int ipgs = 0;

    // Create the Gibbs sampler (multi-mono case)

    AGibbs *gibbs = GibbsFactory::createGibbs(dbin, modvec, rule->getRho(), false);
    gibbs->init(npgs, ngrf, gibbs_nburn, gibbs_niter, seed);

    /* Allocate the covariance matrix inverted */

    if (gibbs->covmatAlloc(verbose)) goto label_end;

    /* Allocate the Gaussian vector */

    VectorVectorDouble y = gibbs->allocY();

    /* Loop on the simulations */

    for (int isimu = 0; isimu < nbsimu; isimu++)
    {
      for (int igrf = 0; igrf < ngrf; igrf++)
        if (rule->evaluateBounds(propdef, dbin, dbout, isimu, igrf, ipgs,
                                 nbsimu)) goto label_end;

      if (gibbs->run(y, ipgs, isimu, verbose)) goto label_end;
    }
  }

  /***************************************************/
  /* Perform the conditional simulation for each GRF */
  /***************************************************/

  local_seed = seed;
  for (int igrf = 0; igrf < 2; igrf++)
  {
    if (!flag_used[igrf]) continue;
    icase = get_rank_from_propdef(propdef, 0, igrf);
    CalcSimuTurningBands situba(nbsimu, nbtuba, flag_check, local_seed);
    local_seed = 0;
    if (situba.simulate(dbin, dbout, models[igrf], neigh, icase, false,
                        VectorDouble(), MatrixSquareSymmetric(), true)) goto label_end;
  }

  /* Convert gaussian to facies at target point */

  if (!flag_gaus)
  {
    for (int isimu = 0; isimu < nbsimu; isimu++)
      simu_func_categorical_transf(dbout, 0, isimu, nbsimu);
  }

  /* Update facies proportions at target points */

  if (flag_prop)
  {
    for (int isimu = 0; isimu < nbsimu; isimu++)
      simu_func_categorical_update(dbout, 0, isimu, nbsimu);
    simu_func_categorical_scale(dbout, 0, nbsimu);
  }

  /* Check/show facies at data against facies at the closest grid node */

  if (flag_cond && !flag_gaus && (flag_check || flag_show))
    st_check_facies_data2grid(dbin, dbout, flag_check, flag_show, 0, nechin,
                              nfacies, nbsimu);

  /********************************/
  /* Free the temporary variables */
  /********************************/

  if (dbout != nullptr)
  {
    if (!st_keep(flag_gaus, flag_prop, RESULT, TYPE_PROP) && iptr_RP >= 0)
      dbout->deleteColumnsByLocator(ELoc::P);
    else
      namconv.setNamesAndLocators(NULL, VectorString(), ELoc::Z, -1, dbout, iptr_RP, "Props",
                                  nbsimu, false);

    if (!st_keep(flag_gaus, flag_prop, RESULT, TYPE_GAUS))
      dbout->deleteColumnsByLocator(ELoc::SIMU);
    else
      namconv.setNamesAndLocators(NULL, VectorString(), ELoc::Z, -1, dbout, iptr_RN, "Gaus",
                                  ngrf * nbsimu, false);

    if (!st_keep(flag_gaus, flag_prop, RESULT, TYPE_FACIES) && iptr_RF >= 0)
      dbout->deleteColumnsByLocator(ELoc::FACIES);
    else
      namconv.setNamesAndLocators(NULL, VectorString(), ELoc::Z, -1, dbout, iptr_RF, String(),
                                  nbsimu);
  }

  if (dbin != nullptr)
  {
    if (!st_keep(flag_gaus, flag_prop, DATA, TYPE_GAUS))
      dbin->deleteColumnsByLocator(ELoc::SIMU);
    else
      namconv.setNamesAndLocators(NULL, VectorString(), ELoc::Z, -1, dbin, iptr_DN, "Gaus",
                                  ngrf * nbsimu, false);

    if (!st_keep(flag_gaus, flag_prop, DATA, TYPE_FACIES) && iptr_DF >= 0)
      dbin->deleteColumnsByLocator(ELoc::FACIES);
    else
      namconv.setNamesAndLocators(NULL, VectorString(), ELoc::Z, -1, dbin, iptr_DF, String(),
                                  nbsimu, false);

    dbin->deleteColumnsByLocator(ELoc::GAUSFAC);
    dbin->deleteColumnsByLocator(ELoc::L);
    dbin->deleteColumnsByLocator(ELoc::U);
  }

  /* Set the error return flag */

  error = 0;

  label_end:
    proportion_manage(-1, 1, flag_stat, ngrf, 0, nfacies, 0,
                      dbin, dbprop, propcst, propdef);
  st_suppress_added_samples(dbin, nechin);
  return (error);
}

/****************************************************************************/
/*!
 **  Perform the conditional or non-conditional Bi Pluri-gaussian
 **  simulations
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure (optional)
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  ruleprop    Ruleprop definition
 ** \param[in]  model11     First Model structure for First Lithotype Rule
 ** \param[in]  model12     Second Model structure for First Lithotype Rule
 ** \param[in]  model21     First Model structure for Second Lithotype Rule
 ** \param[in]  model22     Second Model structure for Second Lithotype Rule
 ** \param[in]  neigh       ANeigh structure
 ** \param[in]  nbsimu      Number of simulations
 ** \param[in]  seed        Seed for random number generator
 ** \param[in]  flag_gaus   1 gaussian results; otherwise facies
 ** \param[in]  flag_prop   1 for facies proportion
 ** \param[in]  flag_check  1 if the facies at data must be checked against
 **                         the closest simulated grid node
 ** \param[in]  flag_show   1 if the grid node which coincides with the data
 **                         should be represented with the data facies
 **                         (only if flag_cond && !flag_gaus)
 ** \param[in]  nbtuba      Number of turning bands
 ** \param[in]  gibbs_nburn Number of bootstrap iterations
 ** \param[in]  gibbs_niter Maximum number of iterations
 ** \param[in]  percent     Amount of nugget effect added to too continuous
 **                         model (expressed in percentage of the total variance)
 ** \param[in]  namconv     Naming convention
 **
 ** \remark  When conditional, the two first variables in the input Db
 ** \remark  should correspond to the two facies indices (starting from 1)
 ** \remark  The argument 'dbin' is optional: it must be defined only for
 ** \remark  conditional simulations
 ** \remark  The proportions (nfac1 * nfac2) must be ordered as follows:
 ** \remark  f1af2a, f1bf2a, f1cf2a, ..., f1bf2a, f1bf2b, ..., f1nf2m
 **
 *****************************************************************************/
int simbipgs(Db *dbin,
             Db *dbout,
             RuleProp *ruleprop,
             Model *model11,
             Model *model12,
             Model *model21,
             Model *model22,
             ANeigh *neigh,
             int nbsimu,
             int seed,
             int flag_gaus,
             int flag_prop,
             int flag_check,
             int flag_show,
             int nbtuba,
             int gibbs_nburn,
             int gibbs_niter,
             double percent,
             const NamingConvention& namconv)
{
  int     iptr,iatt_z[2];
  int     npgs,flag_cond,error,icase;
  int     nfac[2],nfactot,flag_used[2][2],nechin,ngrf[2],ngrftot;
  int     iptr_RP,iptr_RF,iptr_DF,iptr_RN,iptr_DN, local_seed;
  bool    verbose;
  Rule   *rules[2];
  Model  *models[2][2];
  std::vector<Model *> modvec[2];
  PropDef *propdef;

  /* Initializations */

  error = 1;
  npgs = 2;
  nechin = 0;
  verbose = false;
  propdef = nullptr;

  if (ruleprop == nullptr)
  {
    messerr("RuleProp must be defined");
    return 1;
  }
  if (ruleprop->getRule(0)->getModeRule() != ERule::STD)
  {
    messerr("SimuBiPgs is restricted to Standard Lithotype Rule");
    return 1;
  }
  int flag_stat = ruleprop->isFlagStat();
  Rule rule1(*ruleprop->getRule(0));
  Rule rule2(*ruleprop->getRule(1));
  const VectorDouble &propcst = ruleprop->getPropCst();
  const Db *dbprop = ruleprop->getDbprop();

  nfac[0] = rule1.getFaciesNumber();
  nfac[1] = rule2.getFaciesNumber();
  rules[0] = &rule1;
  rules[1] = &rule2;
  models[0][0] = model11;
  models[0][1] = model12;
  models[1][0] = model21;
  models[1][1] = model22;
  nfactot = nfac[0] + nfac[1];
  flag_cond = (dbin != nullptr);
  iptr_RP   = iptr_RF = iptr_DF = iptr_RN = iptr_DN = 0;
  iptr      = -1;
  for (int ipgs=0; ipgs<2; ipgs++)
  {
    ngrf[ipgs] = 0;
    iatt_z[ipgs] = -1;
  }
  if (rules[0]->particularities(dbout, dbprop, model11, 1, flag_stat))
    goto label_end;
  if (rules[1]->particularities(dbout, dbprop, model21, 1, flag_stat))
    goto label_end;

  /**********************/
  /* Preliminary checks */
  /**********************/

  /* Input Db */
  if (flag_cond)
  {
    nechin = dbin->getSampleNumber();
    if (!dbin->isVariableNumberComparedTo(2)) goto label_end;
    iatt_z[0] = dbin->getUIDByLocator(ELoc::Z, 0);
    iatt_z[1] = dbin->getUIDByLocator(ELoc::Z, 1);
  }

  /* Output Db */
  if (flag_prop && flag_gaus)
  {
    messerr(
        "Calculating the facies proportions is incompatible with storing the Gaussian values");
    goto label_end;
  }

  /* Model */

  ngrftot = 0;
  for (int ipgs=0; ipgs<npgs; ipgs++)
  {
    ngrf[ipgs] = rules[ipgs]->getGRFNumber();
    ngrftot += ngrf[ipgs];

    /* Check the validity of the model */

    for (int igrf=0; igrf<2; igrf++)
    {
      flag_used[ipgs][igrf] = rules[ipgs]->isYUsed(igrf);
      if (!flag_used[ipgs][igrf]) continue;
      if (models[ipgs][igrf] == nullptr)
      {
        messerr("Variable #%d needs the underlying GRF #%d", ipgs + 1,
                igrf + 1);
        messerr("No corresponding Model is provided");
        goto label_end;
      }
      if (models[ipgs][igrf]->getVariableNumber() != 1)
      {
        messerr(
            "The number of variables in Model #%d (%d) for Variable %d should be 1",
            igrf + 1, ipgs + 1, models[ipgs][igrf]->getVariableNumber());
        goto label_end;
      }
      if (models[ipgs][igrf]->stabilize(percent, true)) goto label_end;
      if (models[ipgs][igrf]->standardize(true)) goto label_end;

      modvec[ipgs].push_back(models[ipgs][igrf]);
    }
  }

  /* Neighborhood */
  if (neigh->getType() != ENeigh::UNIQUE && neigh->getType() != ENeigh::BENCH)
  {
    messerr("The only authorized Neighborhoods are UNIQUE or BENCH");
    goto label_end;
  }

  /* Rules */

  for (int ipgs=0; ipgs<npgs; ipgs++)
  {
    // Check the Rules (only ERule::STD case is authorized)
    if (rules[ipgs]->getModeRule() != ERule::STD)
    {
      messerr("In the Bi-PGS application, only Standard Rule is authorized");
    }
  }

  /* Final checks */

  for (int ipgs=0; ipgs<2; ipgs++)
  {
    if (flag_cond)
    {
      dbin->clearLocators(ELoc::Z);
      dbin->setLocatorByUID(iatt_z[ipgs], ELoc::Z);
    }
    if (st_check_simtub_environment(dbin, dbout, models[ipgs][0], neigh))
      goto label_end;
  }

  /* Core allocation */

  propdef = proportion_manage(1, 1, flag_stat, ngrf[0], ngrf[1], nfac[0],
                              nfac[1], dbin, dbprop, propcst, propdef);
  if (propdef == nullptr) goto label_end;
  simu_define_func_update(simu_func_categorical_update);
  simu_define_func_scale(simu_func_categorical_scale);
  ModCat.propdef = propdef;

  /**********************/
  /* Add the attributes */
  /**********************/

  /* Storage of the proportions */
  if (flag_prop)
  {
    if (db_locator_attribute_add(dbout, ELoc::P, nfactot, 0, 0., &iptr_RP))
      goto label_end;
  }

  /* Storage of the facies simulations in the output file */
  if (db_locator_attribute_add(dbout, ELoc::FACIES, npgs * nbsimu, 0, 0.,
                               &iptr_RF)) goto label_end;

  /* Storage of the facies simulations in the input file */
  if (flag_cond)
  {
    if (db_locator_attribute_add(dbin, ELoc::FACIES, npgs * nbsimu, 0, 0.,
                                 &iptr_DF)) goto label_end;
  }

  /* Gaussian transform of the facies input data */
  if (flag_cond)
  {
    if (db_locator_attribute_add(dbin, ELoc::GAUSFAC, ngrftot * nbsimu, 0, 0.,
                                 &iptr)) goto label_end;
  }

  /* Non-conditional gaussian simulations at data points */
  if (flag_cond)
  {
    if (db_locator_attribute_add(dbin, ELoc::SIMU, ngrftot * nbsimu, 0, 0.,
                                 &iptr_DN)) goto label_end;
  }

  /* Non-conditional gaussian simulations at target points */
  if (db_locator_attribute_add(dbout, ELoc::SIMU, ngrftot * nbsimu, 0, 0.,
                               &iptr_RN)) goto label_end;

  if (flag_cond)
  {
    /* Lower bound at input data points */
    if (db_locator_attribute_add(dbin, ELoc::L, ngrftot, 0, 0., &iptr))
      goto label_end;

    /* Upper bound at input data points */
    if (db_locator_attribute_add(dbin, ELoc::U, ngrftot, 0, 0., &iptr))
      goto label_end;
  }

  /* Define the environment variables for printout */

  /************************/
  /* Main loop on the PGS */
  /************************/

  for (int ipgs=0; ipgs<npgs; ipgs++)
  {
    if (flag_cond)
    {
      dbin->clearLocators(ELoc::Z);
      dbin->setLocatorByUID(iatt_z[ipgs], ELoc::Z);
    }

    if (ipgs == 0)
      proportion_rule_process(propdef, EProcessOper::MARGINAL);
    else
      proportion_rule_process(propdef, EProcessOper::CONDITIONAL);

    ModCat.rule = rules[ipgs];
    ModCat.ipgs = ipgs;
    ModCat.flag_used[0] = flag_used[ipgs][0];
    ModCat.flag_used[1] = flag_used[ipgs][1];

    /****************************************/
    /* Convert facies into gaussian at data */
    /****************************************/

    if (flag_cond)
    {

      // Create the Gibbs sampler

      AGibbs *gibbs = GibbsFactory::createGibbs(dbin, modvec[ipgs],
                                                rules[ipgs]->getRho(), false);
      gibbs->init(npgs, ngrf[ipgs], gibbs_nburn, gibbs_niter, seed);

      /* Allocate the covariance matrix inverted */

      if (gibbs->covmatAlloc(verbose)) goto label_end;

      // Core allocation

      VectorVectorDouble y = gibbs->allocY();

      /* Loop on the simulations */

      for (int isimu = 0; isimu < nbsimu; isimu++)
      {

        /* Update the proportions */

        for (int igrf = 0; igrf < ngrf[ipgs]; igrf++)
          if (rules[ipgs]->evaluateBounds(propdef, dbin, dbout, isimu, igrf,
                                          ipgs, nbsimu)) goto label_end;

        if (gibbs->run(y, ipgs, isimu, verbose)) goto label_end;
      }

      /* Convert gaussian to facies on data point */

      for (int isimu = 0; isimu < nbsimu; isimu++)
      {
        if (rules[ipgs]->gaus2facData(propdef, dbin, dbout, flag_used[ipgs],
                                      ipgs, isimu, nbsimu)) goto label_end;
      }
    }

    /***************************************************/
    /* Perform the conditional simulation for each GRF */
    /***************************************************/

    /* Define the environment variables for printout */

    local_seed = seed;
    for (int igrf = 0; igrf < 2; igrf++)
    {
      if (!flag_used[ipgs][igrf]) continue;
      icase = get_rank_from_propdef(propdef, ipgs, igrf);
      CalcSimuTurningBands situba(nbsimu, nbtuba, flag_check, local_seed);
      local_seed = 0;
      if (situba.simulate(dbin, dbout, models[ipgs][igrf], neigh, icase, false,
                          VectorDouble(), MatrixSquareSymmetric(), true)) goto label_end;
    }

    /* Convert gaussian to facies at target point */

    if (!flag_gaus) for (int isimu = 0; isimu < nbsimu; isimu++)
      simu_func_categorical_transf(dbout, 0, isimu, nbsimu);

    /* Update facies proportions at target points */

    if (flag_prop)
    {
      for (int isimu = 0; isimu < nbsimu; isimu++)
        simu_func_categorical_update(dbout, 0, isimu, nbsimu);
      simu_func_categorical_scale(dbout, 0, nbsimu);
    }

    /* Check/show facies at data against facies at the closest grid node */

    if (flag_cond && !flag_gaus && (flag_check || flag_show))
      st_check_facies_data2grid(dbin, dbout, flag_check, flag_show, ipgs,
                                nechin, nfac[ipgs], nbsimu);
  }

  /********************************/
  /* Free the temporary variables */
  /********************************/

  if (dbout != nullptr)
  {
    if (!st_keep(flag_gaus, flag_prop, RESULT, TYPE_PROP) && iptr_RP >= 0)
      dbout->deleteColumnsByLocator(ELoc::P);
    else
      namconv.setNamesAndLocators(NULL, VectorString(), ELoc::Z, -1, dbout, iptr_RP, "Props",
                                  nfactot, false);

    if (!st_keep(flag_gaus, flag_prop, RESULT, TYPE_GAUS) && iptr_RN >= 0)
      dbout->deleteColumnsByLocator(ELoc::SIMU);
    else
      namconv.setNamesAndLocators(NULL, VectorString(), ELoc::Z, -1, dbout, iptr_RN, "Gaus",
                                  ngrftot * nbsimu, false);

    if (!st_keep(flag_gaus, flag_prop, RESULT, TYPE_FACIES) && iptr_RF >= 0)
      dbout->deleteColumnsByLocator(ELoc::FACIES);
    else
      namconv.setNamesAndLocators(NULL, VectorString(), ELoc::Z, -1, dbout, iptr_RF, String(),
                                  npgs * nbsimu);
  }

  if (dbin != nullptr)
  {
    if (!st_keep(flag_gaus, flag_prop, DATA, TYPE_GAUS) && iptr_DN >= 0)
      dbin->deleteColumnsByLocator(ELoc::SIMU);
    else
      namconv.setNamesAndLocators(NULL, VectorString(), ELoc::Z, -1, dbin, iptr_DN, "Gaus",
                                  ngrftot * nbsimu, false);

    if (!st_keep(flag_gaus, flag_prop, DATA, TYPE_FACIES) && iptr_DF >= 0)
      dbin->deleteColumnsByLocator(ELoc::FACIES);
    else
      namconv.setNamesAndLocators(NULL, VectorString(), ELoc::Z, -1, dbin, iptr_DF, String(),
                                  npgs * nbsimu, false);

    dbin->deleteColumnsByLocator(ELoc::GAUSFAC);
    dbin->deleteColumnsByLocator(ELoc::L);
    dbin->deleteColumnsByLocator(ELoc::U);
  }

  /* Set the error return flag */

  error = 0;

label_end:
  st_suppress_added_samples(dbin,nechin);
  proportion_manage(-1, 1, flag_stat, ngrf[0], ngrf[1], nfac[0],
                    nfac[1], dbin, dbprop, propcst, propdef);
  return (error);
}

/****************************************************************************/
/*!
 **  Convert series of simulations to conditional expectation and variance
 **
 ** \return  Error return code
 **
 ** \param[in]  db        Db structure
 ** \param[in]  locatorType    Type of pointer containing the simulations
 ** \param[in]  nbsimu    Number of simulations
 ** \param[in]  nvar      Number of variables
 **
 ** \param[out] iptr_ce_arg   Pointer to the Conditional Expectation attributes
 ** \param[out] iptr_cstd_arg Pointer to the Conditional St. Dev. attributes
 **
 *****************************************************************************/
int db_simulations_to_ce(Db *db,
                         const ELoc& locatorType,
                         int nbsimu,
                         int nvar,
                         int *iptr_ce_arg,
                         int *iptr_cstd_arg)
{
  int error, iptr_ce, iptr_cstd, iptr_nb, nech;
  double value, count, mean, var;

  // Initializations

  error = 1;
  iptr_ce = iptr_cstd = iptr_nb = -1;
  if (db == nullptr) goto label_end;
  nech = db->getSampleNumber();
  if (nbsimu <= 0 || nvar <= 0 || nech <= 0) return (1);

  // Allocate the new attributes:

  iptr_ce = db->addColumnsByConstant(nvar, 0.);
  if (iptr_ce < 0) goto label_end;
  iptr_cstd = db->addColumnsByConstant(nvar, 0.);
  if (iptr_cstd < 0) goto label_end;
  iptr_nb = db->addColumnsByConstant(nvar, 0.);
  if (iptr_nb < 0) goto label_end;

  // Loop on the simulations

  for (int isimu = 0; isimu < nbsimu; isimu++)
  {
    // Loop on the samples

    for (int iech = 0; iech < nech; iech++)
    {
      if (!db->isActive(iech)) continue;

      // Loop on the variables

      for (int ivar = 0; ivar < nvar; ivar++)
      {
        // Arguments 'simu' and 'nvar' are interchanged to keep correct order
        value = db->getSimvar(locatorType, iech, ivar, isimu, 0, nvar, nbsimu);
        if (FFFF(value)) continue;
        db->updArray(iech, iptr_ce + ivar, EOperator::ADD, value);
        db->updArray(iech, iptr_cstd + ivar, EOperator::ADD, value * value);
        db->updArray(iech, iptr_nb + ivar, EOperator::ADD, 1.);
      }
    }
  }

  // Scale the conditional expectation and variance

  for (int iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;

    for (int ivar = 0; ivar < nvar; ivar++)
    {
      count = db->getArray(iech, iptr_nb + ivar);
      if (count <= 0)
      {
        db->setArray(iech, iptr_ce + ivar, TEST);
        db->setArray(iech, iptr_cstd + ivar, TEST);
      }
      else
      {
        mean = db->getArray(iech, iptr_ce + ivar) / count;
        db->setArray(iech, iptr_ce + ivar, mean);
        var = db->getArray(iech, iptr_cstd + ivar) / count - mean * mean;
        var = (var > 0.) ? sqrt(var) : 0.;
        db->setArray(iech, iptr_cstd + ivar, var);
      }
    }
  }

  // Set the error return code

  error = 0;

  label_end:
  db->deleteColumnsByUIDRange(iptr_nb, nvar);
  if (error)
  {
    db->deleteColumnsByUIDRange(iptr_ce, nvar);
    db->deleteColumnsByUIDRange(iptr_cstd, nvar);
    *iptr_ce_arg = -1;
    *iptr_cstd_arg = -1;
  }
  else
  {
    *iptr_ce_arg = iptr_ce;
    *iptr_cstd_arg = iptr_cstd;
  }
  return (error);
}

/****************************************************************************/
/*!
 **  Perform the Gibbs sampler
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Db structure
 ** \param[in]  model       Model structure
 ** \param[in]  nbsimu      Number of simulations
 ** \param[in]  seed        Seed for random number generator
 ** \param[in]  gibbs_nburn Initial number of iterations for bootstrapping
 ** \param[in]  gibbs_niter Maximum number of iterations
 ** \param[in]  flag_moving      1 for Moving
 ** \param[in]  flag_norm   1 if the Model must be normalized
 ** \param[in]  flag_multi_mono  1 for the Multi_mono algorithm
 ** \param[in]  flag_propagation 1 for the propagation algorithm
 ** \param[in]  flag_sym_neigh Deprecated argument
 ** \param[in]  gibbs_optstats   0: No stats - 1: Print - 2: Save Neutral file
 ** \param[in]  percent     Amount of nugget effect added to too continuous
 **                         model (expressed in percentage of total variance)
 ** \param[in]  flag_ce     1 if the conditional expectation
 **                         should be returned instead of simulations
 ** \param[in]  flag_cstd   1 if the conditional standard deviation
 **                         should be returned instead of simulations
 ** \param[in]  verbose     Verbose flag
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int gibbs_sampler(Db *dbin,
                  Model *model,
                  int nbsimu,
                  int seed,
                  int gibbs_nburn,
                  int gibbs_niter,
                  bool flag_moving,
                  bool flag_norm,
                  bool flag_multi_mono,
                  bool flag_propagation,
                  bool flag_sym_neigh,
                  int gibbs_optstats,
                  double percent,
                  bool flag_ce,
                  bool flag_cstd,
                  bool verbose,
                  const NamingConvention &namconv)
{
  DECLARE_UNUSED(flag_sym_neigh);
  int error, iptr, npgs, nvar, iptr_ce, iptr_cstd;
  PropDef *propdef;

  /* Initializations */

  error = 1;
  npgs = 1;
  nvar = 0;
  iptr_ce = iptr_cstd = -1;
  propdef = nullptr;

  /**********************/
  /* Preliminary checks */
  /**********************/

  /* Db */

  if (flag_propagation)
  {
    if (dbin->getIntervalNumber() > 0)
    {
      messerr("The propagation algorithm is incompatible with bounds");
      goto label_end;
    }
  }

  /* Model */

  if (model == nullptr)
  {
    messerr("No Model is provided");
    goto label_end;
  }
  nvar = model->getVariableNumber();
  if (!flag_propagation)
  {
    if (model->stabilize(percent, true)) goto label_end;
  }
  if (flag_norm)
  {
    if (model->standardize(true)) goto label_end;
  }

  /*******************/
  /* Core allocation */
  /*******************/

  propdef = proportion_manage(1, 0, 1, 1, 0, nvar, 0, dbin, NULL,
                              VectorDouble(), propdef);
  if (propdef == nullptr) goto label_end;

  /**********************/
  /* Add the attributes */
  /**********************/

  if (db_locator_attribute_add(dbin, ELoc::GAUSFAC, nbsimu * nvar, 0, 0.,
                               &iptr)) goto label_end;

  /*****************/
  /* Gibbs sampler */
  /*****************/

  {
    AGibbs *gibbs;
    if (!flag_multi_mono)
    {
      gibbs = GibbsFactory::createGibbs(dbin, model, flag_moving);
    }
    else
    {
      std::vector<Model*> modvec;
      modvec.push_back(model);
      gibbs = GibbsFactory::createGibbs(dbin, modvec, 0., flag_propagation);
    }
    if (gibbs == nullptr) goto label_end;
    gibbs->setOptionStats(gibbs_optstats);
    gibbs->init(npgs, nvar, gibbs_nburn, gibbs_niter, seed);

    // Allocate the Gaussian vector

    VectorVectorDouble y = gibbs->allocY();

    /* Allocate the covariance matrix inverted */

    if (gibbs->covmatAlloc(verbose)) goto label_end;

    // Invoke the Gibbs calculator

    for (int isimu = 0; isimu < nbsimu; isimu++)
      if (gibbs->run(y, 0, isimu)) goto label_end;
  }

  /* Convert the simulations */

  if (flag_ce || flag_cstd)
  {
    if (db_simulations_to_ce(dbin, ELoc::GAUSFAC, nbsimu, nvar, &iptr_ce,
                             &iptr_cstd)) goto label_end;

    // We release the attributes dedicated to simulations on Dbout

    if (!flag_ce)
    {
      dbin->deleteColumnsByUIDRange(iptr_ce, nvar);
      iptr_ce = -1;
    }
    if (!flag_cstd)
    {
      dbin->deleteColumnsByUIDRange(iptr_cstd, nvar);
      iptr_cstd = -1;
    }
    dbin->deleteColumnsByLocator(ELoc::GAUSFAC);
  }

  /* Set the error return flag */

  error = 0;
  namconv.setNamesAndLocators(dbin, VectorString(), ELoc::UNKNOWN, nvar, dbin, iptr, String(),
                              nbsimu);

label_end:
  proportion_manage(-1, 0, 1, 1, 0, nvar, 0, dbin, NULL, VectorDouble(), propdef);
  return (error);
}

/****************************************************************************/
/*!
 **  Perform a set of valid conditional or non-conditional simulations
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin         Input Db structure (optional)
 ** \param[in]  dbout        Output Db structure
 ** \param[in]  model        Model structure
 ** \param[in]  neigh        ANeigh structure (optional)
 ** \param[in]  seed         Seed for random number generator
 ** \param[in]  nbtuba       Number of turning bands
 ** \param[in]  nbsimu_min   Minimum number of simulations
 ** \param[in]  nbsimu_quant Additional quantum of simulations
 ** \param[in]  niter_max    Maximum number of iterations
 ** \param[in]  cols         Vector of column indices
 ** \param[in]  func_valid   Testing function
 **
 ** \remarks  This function calls a simulation outcome.
 ** \remarks  It provides a return code:
 ** \remarks  -1: the simulation outcome is not valid and the process must be
 ** \remarks      interrupted gently with no error (case of non-convergence)
 ** \remarks   0: the simulation outcome is not valid
 ** \remarks   1: the simulation outcome is valid and must be kept
 ** \remarks   2: the simulation outcome is valid and its returned version
 ** \remarks      must be kept (this is usefull if func_valid() has modified it).
 ** \remarks
 ** \remarks  The func_valid prototype has the following arguments:
 ** \li       ndim    Space dimension
 ** \li       nx      Array of number of grid meshes along all space direction
 ** \li       dx      Array of grid mesh along all space direction
 ** \li       x0      Array of grid origin along all space direction
 ** \li       nonval  Value corresponding to a missing grid value
 ** \li       percent Percentage of the simulations already validated
 ** \li       tab     Array containing the simulated outcome
 **
 ** \remarks  The following lines give an example of func_valid() which considers
 ** \remarks  a simulation outcome as valid if more than 50% of the valid samples
 ** \remarks  have a positive value.
 **
 ** \code
 **   int func_valid(int flag_grid,int ndim,int nech,
 **                  int *nx,double *dx,double *x0,
 **                  double nonval, double percent, double *tab)
 **  {
 **    double ratio;
 **    int i,npositive,nvalid;
 **
 **    for (i=0; i<nech; i++)
 **      {
 **         if (tab[i] == nonval) continue;
 **         nvalid++;
 **         if (tab[i] > 10.) npositive++;
 **      }
 **      ratio = (nvalid > 0) ? npositive / nvalid : 0.;
 **      return(ratio > 0.5);
 **   }
 ** \endcode
 **
 *****************************************************************************/
int simtub_constraints(Db* dbin,
                       Db* dbout,
                       Model* model,
                       ANeigh* neigh,
                       int seed,
                       int nbtuba,
                       int nbsimu_min,
                       int nbsimu_quant,
                       int niter_max,
                       VectorInt& cols,
                       int (*func_valid)(int flag_grid,
                                         int nDim,
                                         int nech,
                                         int* nx,
                                         double* dx,
                                         double* x0,
                                         double nonval,
                                         double percent,
                                         VectorDouble& tab))
{
  int *nx, iatt, retval, nbtest;
  int error, nbsimu, nvalid, isimu, ndim, iter, nech, flag_grid, i;
  double *dx, *x0, percent;
  VectorDouble tab;

  /* Initializations */

  error = 1;
  law_set_random_seed(seed);
  nx = nullptr;
  dx = x0 = nullptr;
  cols.clear();

  /* Preliminary check */

  flag_grid = dbout->isGrid();
  ndim = dbout->getNDim();
  nech = dbout->getSampleNumber();
  tab.resize(dbout->getSampleNumber());
  if (flag_grid)
  {
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(dbout);
    nx = (int*) mem_alloc(sizeof(int) * ndim, 0);
    if (nx == nullptr) goto label_end;
    dx = (double*) mem_alloc(sizeof(double) * ndim, 0);
    if (dx == nullptr) goto label_end;
    x0 = (double*) mem_alloc(sizeof(double) * ndim, 0);
    if (x0 == nullptr) goto label_end;

    for (i = 0; i < ndim; i++)
    {
      nx[i] = dbgrid->getNX(i);
      dx[i] = dbgrid->getDX(i);
      x0[i] = dbgrid->getX0(i);
    }
  }

  /* Implicit loop on the simulations */

  iatt = dbout->getColumnNumber();
  nvalid = iter = nbtest = 0;
  nbsimu = nbsimu_min + nbsimu_quant;
  while (nvalid < nbsimu_min && iter < niter_max)
  {

    /* Performing the simulations */

    iter++;
    nbtest += nbsimu;
    if (simtub(dbin, dbout, model, neigh, nbsimu, 0, nbtuba, 0)) goto label_end;

    /* Check if the simulated outcomes are valid */

    for (isimu = 0; isimu < nbsimu; isimu++, iatt++)
    {

      /* Load the target simulation into the interface buffer */

      if (db_vector_get_att_sel(dbout, iatt, tab.data())) goto label_end;

      /* Check if the simulation is valid */

      percent = 100. * nvalid / nbsimu_min;
      retval = func_valid(flag_grid, ndim, nech, nx, dx, x0, TEST, percent, tab);
      if (retval == 0)
      {

        /* Delete the current simulation */

        dbout->deleteColumnByUID(iatt);

        /* Interrupt the loop (if requested) */

        if (retval < 0)
        {
          error = 0;
          goto label_end;
        }
      }
      else
      {

        /* The current simulation is accepted */

        if (retval > 1)
        {

          /* Update the vector (optional) */

          dbout->setColumnByUID(tab, iatt);
        }
        cols.push_back(iatt);
        nvalid++;
      }
    }

    /* Optional printout */

    if (OptDbg::query(EDbg::CONVERGE))
      message("Iteration #%2d - Simulations %3d tested, %2d valid\n", iter,
              nbtest, nvalid);

    /* Define the number of simulations for the next batch */

    nbsimu = nbsimu_quant;
  }

  /* Set the error return code */

  error = 0;

  label_end:
  mem_free((char* ) nx);
  mem_free((char* ) dx);
  mem_free((char* ) x0);
  return (error);
}

/****************************************************************************/
/*!
 **  Mask the grid nodes whose value is already too large
 **
 ** \return Number of cells left to be filled
 **
 ** \param[in]  dbout     Output Db structure
 ** \param[in]  seuil     Threshold
 ** \param[in]  scale     Scaling factor for the new simulation
 ** \param[in]  iptrv     Pointer to the max-stable outcome
 ** \param[in]  iptrs     Pointer to the current selection
 **
 *****************************************************************************/
static int st_maxstable_mask(Db *dbout,
                             double seuil,
                             double scale,
                             int iptrv,
                             int iptrs)
{
  int iech, number;
  double valsim;

  for (iech = number = 0; iech < dbout->getSampleNumber(); iech++)
  {
    if (!dbout->isActive(iech)) continue;
    valsim = dbout->getArray(iech, iptrv);
    if (valsim > seuil / scale)
      dbout->setArray(iech, iptrs, 0.);
    else
      number++;
  }
  return (number);
}

/****************************************************************************/
/*!
 **  Combine the simulations of the max-stable process
 **
 ** \param[in]  dbout     Output Db structure
 ** \param[in]  scale     Scaling factor for the new simulation
 ** \param[in]  iter0     Rank of the current iteration
 ** \param[in]  iptrg     Pointer to the newly simulated outcome
 ** \param[in]  iptrv     Pointer to the max-stable outcome
 ** \param[in]  iptrr     Pointer to the max-stable rank outcome
 **
 ** \param[in,out] last   Rank of Iteration where the last grid node is covered
 **
 *****************************************************************************/
static void st_maxstable_combine(Db *dbout,
                                 double scale,
                                 int iter0,
                                 int iptrg,
                                 int iptrv,
                                 int iptrr,
                                 int *last)
{
  int iech;
  double valsim, valold;

  for (iech = 0; iech < dbout->getSampleNumber(); iech++)
  {
    if (!dbout->isActive(iech)) continue;
    valold = dbout->getArray(iech, iptrv);
    valsim = dbout->getArray(iech, iptrg) / scale;
    if (valsim > valold)
    {
      dbout->setArray(iech, iptrv, valsim);
      dbout->setArray(iech, iptrr, iter0);
      (*last) = iter0;
    }
  }
}

/****************************************************************************/
/*!
 **  Perform the non-conditional simulation of the Max-Stable Model
 **
 ** \return  Error return code
 **
 ** \param[in]  dbout     Output Db structure
 ** \param[in]  model     Model structure
 ** \param[in]  ratio     Ratio modifying the range at each iteration
 ** \param[in]  seed      Seed for random number generator
 ** \param[in]  nbtuba    Number of turning bands
 ** \param[in]  flag_simu 1 if the simulation must be stored
 ** \param[in]  flag_rank 1 if the iteration rank must be stored
 ** \param[in]  verbose   Verbose flag
 **
 ** \remarks This function uses a threshold that can be defined using
 ** \remarks keypair mechanism with keyword "MaxStableThresh".
 **
 *****************************************************************************/
int simmaxstable(Db *dbout,
                 Model *model,
                 double ratio,
                 int seed,
                 int nbtuba,
                 int flag_simu,
                 int flag_rank,
                 int verbose)
{
  double tpois, seuil;
  int error, iptrg, iptrv, iptrr, iptrs, niter, nleft, icov, last;
  static double seuil_ref = 5.;

  /* Initializations */

  error = 1;
  iptrv = iptrg = iptrs = iptrr = -1;
  law_set_random_seed(seed);
  if (st_check_simtub_environment(NULL, dbout, model, NULL)) goto label_end;
  seuil = get_keypone("MaxStableThresh", seuil_ref);

  /* Preliminary checks */

  if (model->getVariableNumber() != 1)
  {
    messerr("This feature is limited to the monovariate case");
    goto label_end;
  }
  if (!flag_simu && !flag_rank)
  {
    messerr("You must choose 'flag_simu' or  'flag_rank' or both");
    goto label_end;
  }

  /* Define the environment variables for printout */

  st_simulation_environment();

  /* Add the attributes for storing the results */

  iptrv = dbout->addColumnsByConstant(1, 0.);
  if (iptrv < 0) goto label_end;
  iptrr = dbout->addColumnsByConstant(1, 0.);
  if (iptrr < 0) goto label_end;
  if (db_locator_attribute_add(dbout, ELoc::SEL, 1, 0, 0., &iptrs))
    goto label_end;
  if (db_locator_attribute_add(dbout, ELoc::SIMU, 1, 0, 0., &iptrg))
    goto label_end;

  /* Implicit loop on the simulations */

  if (verbose)
  {
    message("Total number of cells = %d\n", dbout->getSampleNumber());
    message("Maximum simulation value = %lf\n", seuil);
  }

  tpois = 0.;
  niter = last = 0;
  while (1)
  {
    niter++;
    tpois -= log(law_uniform(0., 1.));

    /* Mask the nodes that cannot be accessed anymore */

    nleft = st_maxstable_mask(dbout, seuil, tpois, iptrv, iptrs);
    if (nleft <= 0) break;

    /* Processing the Turning Bands algorithm */

    {
      CalcSimuTurningBands situba(1, nbtuba, false, seed);
      if (situba.simulate(nullptr, dbout, model, nullptr, 0)) goto label_end;
    }

    /* Combine the newly simulated outcome to the background */

    st_maxstable_combine(dbout, tpois, niter, iptrg, iptrv, iptrr, &last);

    if (verbose)
      message("Iteration #%2d - Scale = %6.3lf - Nb. cells left = %d\n", niter,
              tpois, nleft);

    /* Update the model for next iteration */

    for (icov = 0; icov < model->getCovaNumber(); icov++)
      model->setRangeIsotropic(icov, model->getRange(icov) * ratio);
  }

  if (verbose)
  {
    message("Number of iterations = %d\n", niter);
    message("Rank of the last covering iteration = %d\n", last);
  }

  /* Set the error return flag */

  error = 0;

  label_end: if (iptrs >= 0) dbout->deleteColumnByUID(iptrs);
  if (iptrg >= 0) dbout->deleteColumnByUID(iptrg);
  if (!flag_rank && iptrr >= 0) dbout->deleteColumnByUID(iptrr);
  if (!flag_simu && iptrv >= 0) dbout->deleteColumnByUID(iptrv);
  return (error);
}

/****************************************************************************/
/*!
 **  Calculate the quantile for a given array
 **
 ** \return  Quantile value
 **
 ** \param[in]  dbout     Output Db structure
 ** \param[in]  proba     Probability
 **
 ** \param[out] sort      Sorting array
 **
 *****************************************************************************/
static double st_quantile(Db *dbout, double proba, double *sort)
{
  int iech, nech, nval, rank;

  /* Initializations */

  nech = dbout->getSampleNumber();

  /* Load the non-masked simulated values */

  for (iech = nval = 0; iech < nech; iech++)
  {
    if (!dbout->isActive(iech)) continue;
    sort[nval++] = dbout->getSimvar(ELoc::SIMU, iech, 0, 0, 0, 1, 1);
  }

  /* Sorting the array */

  ut_sort_double(0, nval, NULL, sort);

  /* Calculate the quantile */

  rank = (int) (proba * (double) nval);
  return (sort[rank]);
}

/****************************************************************************/
/*!
 **  Perform the non-conditional simulation of the Orthogonal Residual Model
 **
 ** \return  Error return code
 **
 ** \param[in]  dbout     Output Db structure
 ** \param[in]  model     Model structure
 ** \param[in]  ncut      Number of cutoffs
 ** \param[in]  zcut      Array of cutoffs
 ** \param[in]  wcut      Array of weights
 ** \param[in]  seed      Seed for random number generator
 ** \param[in]  nbtuba    Number of turning bands
 ** \param[in]  verbose   Verbose flag
 **
 *****************************************************************************/
int simRI(Db *dbout,
          Model *model,
          int ncut,
          double *zcut,
          double *wcut,
          int seed,
          int nbtuba,
          int verbose)
{
  double *pres, *pton, *sort, cumul, simval, proba, seuil;
  int icut, error, iptrg, iptrs, nech, iech, count, total;

  /* Initializations */

  error = 1;
  iptrg = iptrs = -1;
  pres = pton = sort = nullptr;
  nech = dbout->getSampleNumber();
  law_set_random_seed(seed);
  if (st_check_simtub_environment(NULL, dbout, model, NULL)) goto label_end;

  /* Preliminary checks */

  if (model->getVariableNumber() != 1)
  {
    messerr("This feature is limited to the monovariate case");
    goto label_end;
  }

  /* Define the environment variables for printout */

  st_simulation_environment();

  /* Add the attributes for storing the results */

  sort = (double*) mem_alloc(sizeof(double) * nech, 0);
  if (sort == nullptr) goto label_end;
  pton = (double*) mem_alloc(sizeof(double) * ncut, 0);
  if (pton == nullptr) goto label_end;
  pres = (double*) mem_alloc(sizeof(double) * (ncut - 1), 0);
  if (pres == nullptr) goto label_end;
  if (db_locator_attribute_add(dbout, ELoc::SEL, 1, 0, 0., &iptrs))
    goto label_end;
  if (db_locator_attribute_add(dbout, ELoc::SIMU, 1, 0, 0., &iptrg))
    goto label_end;

  /* Preliminary calculations */

  cumul = 0.;
  for (icut = 0; icut < ncut; icut++)
  {
    if (icut > 0 && zcut[icut] <= zcut[icut - 1])
    {
      messerr("The cutoff values must be ordered increasingly");
      goto label_end;
    }
    if (wcut[icut] < 0)
    {
      messerr("The weight of class (%d) cannot be negative", icut + 1);
      goto label_end;
    }
    cumul += wcut[icut];
  }
  if (cumul <= 0.)
  {
    messerr("The sum of weights cannot be negative or null");
    goto label_end;
  }
  for (icut = 0; icut < ncut; icut++)
    wcut[icut] /= cumul;
  pton[0] = 1.;
  for (icut = 1; icut < ncut; icut++)
    pton[icut] = pton[icut - 1] - wcut[icut];
  for (icut = 0; icut < ncut - 1; icut++)
    pres[icut] = pton[icut + 1] / pton[icut];

  /* Set the mask to the whole set of grid nodes */

  for (iech = 0; iech < nech; iech++)
    dbout->setLocVariable(ELoc::SEL,iech, 0, 1.);

  /* Loop on the cutoff classes */

  total = 0;
  for (icut = 0; icut < ncut; icut++)
  {

    /* Simulation in the non-masked part of the grid */

    {
      CalcSimuTurningBands situba(1, nbtuba, false, seed);
      if (situba.simulate(nullptr, dbout, model, nullptr, 0)) goto label_end;
    }

    /* Look for the quantile */

    proba = 1. - pres[icut];
    seuil = (icut < ncut - 1) ? st_quantile(dbout, proba, sort) : TEST;

    /* Update the current selection */

    for (iech = count = 0; iech < nech; iech++)
    {
      if (!dbout->getSelection(iech)) continue;
      simval = dbout->getSimvar(ELoc::SIMU, iech, 0, 0, 0, 1, 1);
      if (!FFFF(seuil) && simval >= seuil) continue;
      dbout->setSimvar(ELoc::SIMU, iech, 0, 0, 0, 1, 1, (double) (icut + 1));
      dbout->setLocVariable(ELoc::SEL,iech, 0, 0.);
      count++;
    }
    total += count;
    if (verbose)
      message("Level %3d - Proba=%lf - Affected=%7d - Total=%7d\n", icut + 1,
              proba, count, total);
  }

  /* Perform the final coding */

  for (iech = 0; iech < nech; iech++)
  {
    icut = (int) dbout->getSimvar(ELoc::SIMU, iech, 0, 0, 0, 1, 1);
    if (icut < 1 || icut > ncut)
      dbout->setSimvar(ELoc::SIMU, iech, 0, 0, 0, 1, 1, TEST);
    else
      dbout->setSimvar(ELoc::SIMU, iech, 0, 0, 0, 1, 1, zcut[icut - 1]);
  }

  /* Set the error return flag */

  error = 0;

  label_end:
  mem_free((char* ) sort);
  mem_free((char* ) pton);
  mem_free((char* ) pres);
  if (iptrs >= 0) dbout->deleteColumnByUID(iptrs);
  return (error);
}

/****************************************************************************/
/*!
 **  Perform the conditional Pluri-gaussian simulations using spde
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure (optional)
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  ruleprop    RuleProp definition
 ** \param[in]  model1      First Model structure
 ** \param[in]  model2      Second Model structure (optional)
 ** \param[in]  triswitch   Meshing option
 ** \param[in]  gext        Array of domain dilation
 ** \param[in]  flag_gaus   1 if results must be gaussian; otherwise facies
 ** \param[in]  flag_prop   1 for facies proportion
 ** \param[in]  flag_check  1 if the facies at data must be checked against
 **                         the closest simulated grid node
 ** \param[in]  flag_show   1 if the grid node which coincides with the data
 **                         should be represented with the data facies
 ** \param[in]  nfacies     Number of facies
 ** \param[in]  seed        Seed for random number generator
 ** \param[in]  nbsimu      Number of simulations
 ** \param[in]  gibbs_nburn Number of iterations (Burning step)
 ** \param[in]  gibbs_niter Maximum number of iterations
 ** \param[in]  ngibbs_int  Number of iterations internal to Gibbs (SPDE)
 ** \param[in]  verbose     Verbose flag
 ** \param[in]  percent     Amount of nugget effect added to too continous
 **                         model (expressed in percentage of the total variance)
 **
 ** \remark  When conditional, the unique variable in the input Db structure
 ** \remark  should correspond to the facies index (starting from 1)
 **
 *****************************************************************************/
int simpgs_spde(Db* dbin,
                Db* dbout,
                RuleProp* ruleprop,
                Model* model1,
                Model* model2,
                const String& triswitch,
                const VectorDouble& gext,
                int flag_gaus,
                int flag_prop,
                int flag_check,
                int flag_show,
                int nfacies,
                int seed,
                int nbsimu,
                int gibbs_nburn,
                int gibbs_niter,
                int ngibbs_int,
                int verbose,
                double percent)
{
  int iptr, ngrf, igrf, nechin, error, flag_used[2], flag_cond;
  int iptr_RF, iptr_RP;
  Model *models[2];
  PropDef *propdef;
  SPDE_Option s_option;

  /* Initializations */

  error = 1;
  nechin = 0;
  ngrf = 0;
  propdef = nullptr;
  models[0] = model1;
  models[1] = model2;
  iptr_RF = iptr_RP = 0;
  iptr = -1;
  flag_cond = (dbin != nullptr);
  law_set_random_seed(seed);

  if (ruleprop == nullptr)
  {
    messerr("RuleProp must be defined");
    return 1;
  }
  int flag_stat = ruleprop->isFlagStat();
  const Rule *rule = ruleprop->getRule();
  const VectorDouble &propcst = ruleprop->getPropCst();
  const Db *dbprop = ruleprop->getDbprop();

  if (rule->getModeRule() == ERule::SHADOW)
  {
    messerr("The 'Shadow' rule is not authorized");
    goto label_end;
  }
  if (rule->particularities(dbout, dbprop, model1, 1, flag_stat))
    goto label_end;

  if (isGlobalFlagEigen())
  {
    messerr("This method is not coded (yet) for Eigen matrices");
    messerr("Please use setGlobalFlagEigen(false)");
    goto label_end;
  }

  /**********************/
  /* Preliminary checks */
  /**********************/

  /* Input Db */
  if (flag_cond)
  {
    nechin = dbin->getSampleNumber();
    if (!dbin->isVariableNumberComparedTo(1)) goto label_end;
  }

  /* Output Db */
  if (dbout == nullptr)
  {
    messerr("'dbout' is compulsory");
    goto label_end;
  }
  if (flag_prop && flag_gaus)
  {
    messerr(
        "Calculating the facies proportions is incompatible with storing the Gaussian values");
    goto label_end;
  }

  /* Model */
  for (igrf = 0; igrf < 2; igrf++)
  {
    flag_used[igrf] = rule->isYUsed(igrf);
    if (!flag_used[igrf]) continue;
    ngrf++;
    if (models[igrf] == nullptr)
    {
      messerr("The Underlying GRF #%d is needed", igrf + 1);
      messerr("No corresponding Model is provided");
      goto label_end;
    }
    if (models[igrf]->getVariableNumber() != 1)
    {
      messerr("The number of variables in the model #%d (%d) should be 1",
              igrf + 1, model1->getVariableNumber());
      goto label_end;
    }
    if (models[igrf]->stabilize(percent, true)) goto label_end;
    if (models[igrf]->standardize(true)) goto label_end;
  }
  if (spde_check(dbin, dbout, model1, model2, verbose, gext, true, true, true,
                 false, false, true, flag_prop)) goto label_end;
  s_option = spde_option_alloc();
  spde_option_update(s_option, triswitch);

  /* Define the environment variables for printout */

  /**********************/
  /* Add the attributes */
  /**********************/

  /* Storage of the facies proportions */
  if (flag_prop)
  {
    if (db_locator_attribute_add(dbout, ELoc::P, nfacies, 0, 0., &iptr_RP))
      goto label_end;
  }

  /* Storage of the simulations in the Output Db */
  if (db_locator_attribute_add(dbout, ELoc::SIMU, nbsimu * ngrf, 0, 0.,
                               &iptr_RF)) goto label_end;

  if (flag_cond)
  {
    /* Lower bound at input data points */
    if (db_locator_attribute_add(dbin, ELoc::L, ngrf, 0, 0., &iptr))
      goto label_end;

    /* Upper bound at input data points */
    if (db_locator_attribute_add(dbin, ELoc::U, ngrf, 0, 0., &iptr))
      goto label_end;
  }

  /* Storage of the facies simulations in the Output Db */
  if (db_locator_attribute_add(dbout, ELoc::FACIES, nbsimu, 0, 0., &iptr_RF))
    goto label_end;

  propdef = proportion_manage(1, 1, flag_stat, ngrf, 0, nfacies, 0, dbin,
                              dbprop, propcst, propdef);
  if (propdef == nullptr) goto label_end;
  if (!flag_gaus) simu_define_func_transf(simu_func_categorical_transf);
  simu_define_func_update(simu_func_categorical_update);
  simu_define_func_scale(simu_func_categorical_scale);
  ModCat.propdef = propdef;
  ModCat.rule = rule;
  ModCat.ipgs = 0;
  ModCat.flag_used[0] = flag_used[0];
  ModCat.flag_used[1] = flag_used[1];

  /****************************************/
  /* Convert facies into gaussian at data */
  /****************************************/

  proportion_rule_process(propdef, EProcessOper::COPY);

  /* Initialize the Gibbs calculations */

  st_init_gibbs_params(rule->getRho());

  if (flag_cond)
    for (igrf = 0; igrf < 2; igrf++)
    {
      if (!flag_used[igrf]) continue;
      for (int isimu = 0; isimu < nbsimu; isimu++)
      {
        if (rule->evaluateBounds(propdef, dbin, dbout, isimu, igrf, 0, nbsimu))
          goto label_end;
      }
    }

  /***************************************************/
  /* Perform the conditional simulation for each GRF */
  /***************************************************/

  if (spde_prepar(dbin, dbout, gext, s_option)) goto label_end;
  if (spde_process(dbin, dbout, s_option, nbsimu, gibbs_nburn, gibbs_niter,
                   ngibbs_int)) goto label_end;

  /* Check/show facies at data against facies at the closest grid node */

  if (flag_cond && !flag_gaus && (flag_check || flag_show))
    st_check_facies_data2grid(dbin, dbout, flag_check, flag_show, 0, nechin,
                              nfacies, nbsimu);

  /********************************/
  /* Free the temporary variables */
  /********************************/

  if (!st_keep(flag_gaus, flag_prop, RESULT, TYPE_PROP) && iptr_RP >= 0)
    dbout->deleteColumnsByLocator(ELoc::P);

  if (!st_keep(flag_gaus, flag_prop, RESULT, TYPE_FACIES) && iptr_RF >= 0)
    dbout->deleteColumnsByLocator(ELoc::FACIES);

  if (!st_keep(flag_gaus, flag_prop, RESULT, TYPE_GAUS))
    dbout->deleteColumnsByLocator(ELoc::SIMU);

  if (dbin != nullptr)
  {
    dbin->deleteColumnsByLocator(ELoc::L);
    dbin->deleteColumnsByLocator(ELoc::U);
  }

  /* Set the error return flag */

  error = 0;

  label_end:
  proportion_manage(-1, 1, flag_stat, ngrf, 0, nfacies, 0,
                    dbin, dbprop, propcst, propdef);
  st_suppress_added_samples(dbin, nechin);
  return (error);
}

/****************************************************************************/
/*!
 **  Perform the conditional simulations under inequality constraints
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure (optional)
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  model       Model structure
 ** \param[in]  seed        Seed for random number generator
 ** \param[in]  nbsimu      Number of simulations
 ** \param[in]  nbtuba      Number of turning bands
 ** \param[in]  gibbs_nburn Initial number of iterations for bootstrapping
 ** \param[in]  gibbs_niter Maximum number of iterations
 ** \param[in]  flag_check  1 to check the proximity in Gaussian scale
 ** \param[in]  flag_ce     1 if the conditional expectation
 **                         should be returned instead of simulations
 ** \param[in]  flag_cstd   1 if the conditional standard deviation
 **                         should be returned instead of simulations
 ** \param[in]  verbose     Verbose flag
 **
 ** \remarks The Neighborhood does not have to be defined as this method
 ** \remarks only functions using a Unique Neighborhood. For consistency
 ** \remarks it is generated internally.
 **
 *****************************************************************************/
int simcond(Db *dbin,
            Db *dbout,
            Model *model,
            int seed,
            int nbsimu,
            int nbtuba,
            int gibbs_nburn,
            int gibbs_niter,
            int flag_check,
            int flag_ce,
            int flag_cstd,
            int verbose)
{
  PropDef *propdef;
  int nvar, error, iptr, iptr_ce, iptr_cstd;

  /* Initializations */

  error = 1;
  bool flag_ext_created = false;
  bool flag_nostat_created = false;
  nvar = model->getVariableNumber();
  iptr = -1;
  propdef = nullptr;

  /* Preliminary checks */

  NeighUnique* neighU = NeighUnique::create(false);
  law_set_random_seed(seed);
  if (st_check_simtub_environment(dbin, dbout, model, NULL)) goto label_end;
  if (manageExternalInformation(1, ELoc::F, dbin, dbout, &flag_ext_created)) goto label_end;
  if (manageExternalInformation(1, ELoc::NOSTAT, dbin, dbout, &flag_nostat_created)) goto label_end;

  /* Limitations */

  if (nvar > 1)
  {
    messerr("This method is restricted to the monovariate case");
    goto label_end;
  }
  if (dbin->getIntervalNumber() <= 0)
  {
    messerr("No bound is defined: use 'simtub' instead");
    goto label_end;
  }

  /* Core allocation */

  propdef = proportion_manage(1, 0, 1, 1, 0, nvar, 0, dbin, NULL,
                              VectorDouble(), propdef);
  if (propdef == nullptr) goto label_end;

  /* Define the environment variables for printout */

  st_simulation_environment();

  /* Add the attributes for storing the results */

  if (db_locator_attribute_add(dbin, ELoc::GAUSFAC, nbsimu, 0, 0., &iptr))
    goto label_end;
  if (db_locator_attribute_add(dbin, ELoc::SIMU, nvar * nbsimu, 0, 0., &iptr))
    goto label_end;
  if (db_locator_attribute_add(dbout, ELoc::SIMU, nvar * nbsimu, 0, 0., &iptr))
    goto label_end;

  /*****************/
  /* Gibbs sampler */
  /*****************/

  {
    AGibbs *gibbs = GibbsFactory::createGibbs(dbin, model, false);
    gibbs->init(1, 1, gibbs_nburn, gibbs_niter, seed);

    /* Allocate the covariance matrix inverted */

    if (gibbs->covmatAlloc(verbose)) goto label_end;

    /* Allocate the Gaussian Vector */

    VectorVectorDouble y = gibbs->allocY();

    /* Loop on the simulations */

    for (int isimu = 0; isimu < nbsimu; isimu++)
    {
      if (gibbs->run(y, 0, isimu, verbose)) goto label_end;
    }
  }

  /* Processing the Turning Bands algorithm */

  {
    CalcSimuTurningBands situba(nbsimu, nbtuba, flag_check, seed);
    if (situba.simulate(dbin, dbout, model, neighU, 0, false,
                        VectorDouble(), MatrixSquareSymmetric(), false, true)) goto label_end;
  }

  /* Free the temporary variables not used anymore */

  dbin->deleteColumnsByLocator(ELoc::GAUSFAC);
  dbin->deleteColumnsByLocator(ELoc::SIMU);

  /* Convert the simulations */

  if (flag_ce || flag_cstd)
  {
    if (db_simulations_to_ce(dbout, ELoc::SIMU, nbsimu, nvar, &iptr_ce,
                             &iptr_cstd)) goto label_end;

    // We release the attributes dedicated to simulations on Dbout

    dbout->deleteColumnsByLocator(ELoc::SIMU);
    if (!flag_ce)
    {
      dbout->deleteColumnsByUIDRange(iptr_ce, nvar);
      iptr_ce = -1;
    }
    if (!flag_cstd)
    {
      dbout->deleteColumnsByUIDRange(iptr_cstd, nvar);
      iptr_cstd = -1;
    }
  }

  /* Set the error return flag */

  error = 0;

  label_end:
  delete neighU;
  (void) manageExternalInformation(-1, ELoc::F, dbin, dbout, &flag_ext_created);
  (void) manageExternalInformation(-1, ELoc::NOSTAT, dbin, dbout, &flag_nostat_created);
  return (error);
}

/*****************************************************************************/
/*!
 **  Simulates the random function on the sphere
 **
 ** \param[in]  db          Data base containing the coordinates of target points
 **                         These coordinates must be expressed in long/lat
 ** \param[in]  model       Model (defined in Euclidean space) to be used
 ** \param[in]  sphepar     SimuSphericalParam structure
 ** \param[in]  seed        Seed for random number generation
 ** \param[in]  verbose     Verbose flag
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int simsph(DbGrid *db,
           Model *model,
           const SimuSphericalParam& sphepar,
           int seed,
           bool verbose,
           const NamingConvention& namconv)
{
  int flag_sphere;

  /* Preliminary checks */

  flag_sphere = isDefaultSpaceSphere();
  if (!flag_sphere)
  {
    messerr("The Spherical Simulation is restricted to Spherical coordinates");
    return 1;
  }
  if (db->getNDim() != 2)
  {
    messerr("The Simulation on Sphere is restricted to 2-D case");
    return 1;
  }
  for (int icova = 0; icova < model->getCovaNumber(); icova++)
  {
    if (model->getCova(icova)->getFlagAniso())
    {
      messerr("Only Isotropic Models may be used for Spherical Simulations");
      return 1;
    }
  }

  /* Create the new variable in the Data base */

  int iptr = db->addColumnsByConstant(1, 0., String(), ELoc::SIMU);

  SimuSpherical simsphe(1, seed);
  if (simsphe.simulate(db, model, sphepar, iptr, verbose)) return 1;

  namconv.setNamesAndLocators(db, VectorString(), ELoc::UNKNOWN, 1, db, iptr, "Simu");
  return 0;
}

/*****************************************************************************/
/*!
 **  Simulates the random function on the sphere
 **
 ** \return The Vector simulated values
 **
 ** \param[in]  mesh        MeshSpherical object
 ** \param[in]  model       Model (defined in Euclidean space) to be used
 ** \param[in]  sphepar     SimuSphericalParam structure
 ** \param[in]  seed        Seed for random number generation
 ** \param[in]  verbose     Verbose flag
 **
 *****************************************************************************/
VectorDouble simsph_mesh(MeshSpherical *mesh,
                         Model *model,
                         const SimuSphericalParam& sphepar,
                         int seed,
                         int verbose)
{
  VectorDouble simu;

  bool flag_sphere = isDefaultSpaceSphere();
  if (!flag_sphere)
  {
    messerr("The Spherical Simulation is restricted to Spherical coordinates");
    return simu;
  }
  for (int icova = 0; icova < model->getCovaNumber(); icova++)
  {
    if (model->getCova(icova)->getFlagAniso())
    {
      messerr("Only Isotropic Models may be used for Spherical Simulations");
      return simu;
    }
  }

  SimuSpherical simsphe(1, seed);
  simu = simsphe.simulate_mesh(mesh, model, sphepar, verbose);

  return simu;
}

/****************************************************************************/
/*!
 **  Refine the simulation
 **
 ** \return  Newly refined Grid.
 **
 ** \param[in]  dbin       Input grid Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  param      SimuRefineParam structure
 ** \param[in]  seed       Seed for the random number generator
 **
 ** \remark For each dimension of the space, if N stands for the number of
 ** \remark nodes in the input grid, the number of nodes of the output grid
 ** \remark will be (N-1) * 2^p + 1 where p is the param.getNmult()
 **
 *****************************************************************************/
DbGrid* simfine(DbGrid *dbin,
                Model *model,
                const SimuRefineParam& param,
                int seed)
{
  /* Preliminary check */

  if (!dbin->isVariableNumberComparedTo(1)) return nullptr;

  /* Patch the model with maximum dimension for OK */

  model->setField(dbin->getExtensionDiagonal());

  // Perform the simulation

  SimuRefine simfine(1, seed);
  DbGrid* dbout = simfine.simulate(dbin, model, param);

  return dbout;
}

/****************************************************************************/
/*!
 **  Check if the sample belongs to the time slice
 **
 ** \return  Rank of the time slice (or -1)
 **
 ** \param[in]  date     Date attached to a sample
 ** \param[in]  ntime    Number of time intervals
 ** \param[in]  time0    Origin of the first time interval
 ** \param[in]  dtime    Time interval
 **
 *****************************************************************************/
static int st_getTimeInterval(double date,
                              int ntime,
                              double time0,
                              double dtime)
{
  for (int itime = 0; itime < ntime; itime++)
  {
    double time_deb = time0 + dtime * itime;
    double time_fin = time0 + dtime * (itime + 1);
    if (date >= time_deb && date < time_fin) return (itime);
  }
  return (-1);
}

static int st_getFACIES(const DbGrid* dbgrid, int nfacies, int indFacies, int iech)
{
  int ifacies = (int) dbgrid->getArray(iech, indFacies);
  if (ifacies < 0 || ifacies > nfacies || IFFFF(ifacies)) ifacies = 0;
  return (ifacies);
}

static double st_getPORO(const DbGrid* dbgrid, int indPoro, int iech)
{
  if (indPoro <= 0) return (1);
  double poro = dbgrid->getArray(iech, indPoro);
  if (FFFF(poro) || poro < 0.) poro = 0.;
  return (poro);
}

static double st_getDATE(const DbGrid* dbgrid, int indDate, int iech)
{
  double date;

  if (indDate <= 0) return (0);
  date = dbgrid->getArray(iech, indDate);
  if (FFFF(date)) return (0);
  date = MAX(1., date);
  return (date);
}

static int st_getFLUID(const DbGrid* dbgrid, int nfluids, int indFluid, int iech)
{
  int ifluid = (int) dbgrid->getArray(iech, indFluid);
  if (ifluid < 0 || ifluid > nfluids || IFFFF(ifluid)) ifluid = 0;
  return (ifluid);
}

/*****************************************************************************/
/*!
**  Extract time charts from the fluid propagation block
**
** \return  The returned matrix
**
** \param[in]  dbgrid        Db grid structure
** \param[in]  name_facies   Name of variable containing Facies
** \param[in]  name_fluid    Name of variable containing Fluid
** \param[in]  name_poro     Name of variable containing Porosity (optional)
** \param[in]  name_date     Name of variable containing Date
** \param[in]  nfacies       number of facies (facies 0 excluded)
** \param[in]  nfluids       number of fluids
** \param[in]  facies0       Value of the target facies
** \param[in]  fluid0        Value of the target fluid
** \param[in]  ntime         Number of Time intervals
** \param[in]  time0         Starting time
** \param[in]  dtime         Time interval
** \param[in]  verbose       1 for a verbose option
**
*****************************************************************************/
MatrixRectangular fluid_extract(DbGrid *dbgrid,
                                const String& name_facies,
                                const String& name_fluid,
                                const String& name_poro,
                                const String& name_date,
                                int nfacies,
                                int nfluids,
                                int facies0,
                                int fluid0,
                                int ntime,
                                double time0,
                                double dtime,
                                bool verbose)
{
  MatrixRectangular tab;
  if (! dbgrid->isGrid())
  {
    messerr("The Fluid Propagation is restricted to regular grid");
    return tab;
  }
  if (dbgrid->getNDim() > 3)
  {
    messerr("Fluid propagation is limited to 3-D space (maximum)");
    return tab;
  }
  if (ntime < 0 || time0 < 0 || dtime <= 0)
  {
    messerr("Error in Time Interval Definition");
    messerr("Origin=%lf - Step=%lf - Number=%d",time0,dtime,ntime);
    return tab;
  }

  /* Define global variables */

  int ind_facies = dbgrid->getUID(name_facies);
  int ind_fluid  = dbgrid->getUID(name_fluid);
  int ind_date   = dbgrid->getUID(name_date);
  if (ind_facies < 0)
  {
    messerr("Variable 'Facies' must be provided");
    return 1;
  }
  if (ind_fluid < 0)
  {
    messerr("Variable 'Fluid' must be provided");
    return 1;
  }
  if (ind_date < 0)
  {
    messerr("Variable 'Date' must be provided");
    return 1;
  }
  int ind_poro = dbgrid->getUID(name_poro);

  /* Initialize the array */

  tab = MatrixRectangular(ntime, 4);
  int nxyz    = dbgrid->getSampleNumber();
  for (int itime = 0; itime < ntime; itime++)
  {
    tab.setValue(itime, 0, time0 + dtime * itime);
    tab.setValue(itime, 1, time0 + dtime * (itime + 1));
    tab.setValue(itime, 2, 0.);
    tab.setValue(itime, 3, 0.);
  }

  /* Loop on the blocks */

  double totnum = 0.;
  double totvol = 0.;
  double locnum = 0.;
  double locvol = 0.;
  double datmax = 0;
  for (int iech = 0; iech < nxyz; iech++)
  {

    if (st_getFACIES(dbgrid, nfacies, ind_facies, iech) != facies0) continue;
    if (st_getFLUID(dbgrid, nfluids, ind_fluid, iech) != fluid0) continue;
    double volume = st_getPORO(dbgrid, ind_poro, iech);
    double date = st_getDATE(dbgrid, ind_date, iech);
    if (date > datmax) datmax = date;

    totnum += 1;
    totvol += volume;
    int itime = st_getTimeInterval(date, ntime, time0, dtime);
    if (itime < 0) continue;
    locnum += 1;
    locvol += volume;

    tab.setValue(itime, 2, tab.getValue(itime, 2) + 1);
    tab.setValue(itime, 3, tab.getValue(itime, 3) + volume);
  }

  /* Final printout */

  if (verbose)
  {
    mestitle(1, "Extraction for Fluid(%d) and Facies(%d)", facies0, fluid0);
    message("Time slices: From %lf to %lf by step of %lf\n",
            time0, time0 + dtime * ntime, dtime);
    message("Total Number of Cells               = %d\n", nxyz);
    message("Maximum Date                        = %lf\n", datmax);
    message("Total Number of Invaded Cells       = %lf\n", totnum);
    message("Total Volume of Invaded Cells       = %lf\n", totvol);
    message("Total Number of Cells in Time Slice = %lf\n", locnum);
    message("Total Volume of Cells in Time Slice = %lf\n", locvol);
  }
  return tab;
}
