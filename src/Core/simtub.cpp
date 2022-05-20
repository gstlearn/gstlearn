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
#include "geoslib_f_private.h"
#include "Gibbs/GibbsUMultiMono.hpp"
#include "Gibbs/GibbsUPropMono.hpp"
#include "Gibbs/GibbsMMulti.hpp"
#include "Gibbs/GibbsFactory.hpp"
#include "Morpho/Morpho.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/ECov.hpp"
#include "Basic/MathFunc.hpp"
#include "Basic/String.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/EJustify.hpp"
#include "Basic/OptCustom.hpp"
#include "LithoRule/PropDef.hpp"
#include "LithoRule/Rule.hpp"
#include "LithoRule/RuleShift.hpp"
#include "LithoRule/RuleShadow.hpp"
#include "LithoRule/RuleProp.hpp"
#include "LithoRule/EProcessOper.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Simulation/SimuTurningBands.hpp"
#include "Simulation/SimuBoolean.hpp"

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

static int FLAG_DGM;
static double GIBBS_RHO, GIBBS_SQR;
static double R_COEFF;
static Modif_Categorical ModCat = { 0, { 0, 0 }, NULL, NULL };

/****************************************************************************/
/*!
 **  Initialize the global values
 **
 *****************************************************************************/
static void st_simulation_environment(void)
{
  FLAG_DGM = 0;
  GIBBS_RHO = 0.;
  GIBBS_SQR = 0.;
  R_COEFF = 0.;
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
  if (ipgs <= 0)
    return (ifac);
  else
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
  iptr_simu = db->getSimvarRank(isimu, 0, 0, nbsimu, 1);

  /* Loop on the grid cells */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    simval = get_LOCATOR_ITEM(db, ELoc::SIMU, iptr_simu, iech);
    db->updVariable(iech, 0, 0, simval);
    db->updVariable(iech, 1, 0, simval * simval);
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
  iptr_simu = db->getSimvarRank(isimu, 0, ipgs, nbsimu, 1);

  /* Loop on the grid cells */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    facies = (int) get_LOCATOR_ITEM(db, ELoc::FACIES, iptr_simu, iech) - 1;
    rank = st_facies(ModCat.propdef, ipgs, facies);
    prop = db->getProportion(iech, rank) + 1.;
    db->setProportion(iech, rank, prop);
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
    mean = db->getVariable(iech, 0) / nbsimu;
    db->setVariable(iech, 0, mean);
    stdv = db->getVariable(iech, 1) / nbsimu - mean * mean;
    stdv = (stdv > 0) ? sqrt(stdv) :
                        0.;
    db->setVariable(iech, 1, stdv);
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
      prop = db->getProportion(iech, rank) / (double) nbsimu;
      db->setProportion(iech, rank, prop);
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
                               const ELoc &locatorType)
{
  if (get_LOCATOR_NITEM(db,locatorType) <= 0)
    messageAbort("%s : Attributes %d are mandatory",method,locatorType.getValue());
  return;
}

/****************************************************************************/
/*!
 **  Check if the field must be kept
 **
 ** \return  1 if the field must be kept; 0 otherwise
 **
 ** \param[in]  flag_gaus  1 gaussian results; otherwise facies
 ** \param[in]  flag_modif 1 for facies proportion
 ** \param[in]  file       DATA or RESULT
 ** \param[in]  type       0 for gaussian; 1 for facies; 2 for proportion
 **
 *****************************************************************************/
static int st_keep(int flag_gaus, int flag_modif, int file, int type)
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
        keep = (flag_gaus && !flag_modif);
        break;

      case 1: /* Facies */
        keep = (!flag_gaus && !flag_modif);
        break;

      case 2: /* Proportion */
        keep = (flag_modif);
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
 ** \param[in]  neighparam ANeighParam structure (optional if non conditional)
 **
 *****************************************************************************/
static int st_check_simtub_environment(Db *dbin,
                                       Db *dbout,
                                       Model *model,
                                       ANeighParam *neighparam)
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
    if (flag_cond && dbin->getVariableNumber() != nvar)
    {
      messerr("The number of variables of the Data (%d)",
              dbin->getVariableNumber());
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

    nfex = model_nfex(model);
    if (flag_cond && nfex != 0 && !is_grid(dbout)
        && dbin->getExternalDriftNumber() != nfex)
    {
      messerr("The Model requires %d external drift(s)", model_nfex(model));
      messerr("but the input Db refers to %d external drift variables",
              dbin->getExternalDriftNumber());
      return 1;
    }
    if (nfex != 0 && dbout->getExternalDriftNumber() != nfex)
    {
      messerr("The Model requires %d external drift(s)", model_nfex(model));
      messerr("but the output Db refers to %d external drift variables",
              dbout->getExternalDriftNumber());
      return 1;
    }
  }

  /*********************************/
  /* Calculate the field extension */
  /*********************************/

  VectorDouble db_mini(ndim);
  VectorDouble db_maxi(ndim);

  db_extension(dbout, db_mini, db_maxi, false);

  if (flag_cond)
    db_extension(dbin, db_mini, db_maxi, true);

  if (model != nullptr)
    model->setField(ut_vector_extension_diagonal(db_mini, db_maxi));

  /*****************************/
  /* Checking the Neighborhood */
  /*****************************/

  if (flag_cond && neighparam != nullptr)
  {
    if (neighparam->getNDim() != ndim)
    {
      messerr("The Space Dimension of the Neighborhood (%d)", neighparam->getNDim());
      messerr("does not correspond to the Space Dimension of the first Db (%d)",
              ndim);
      return 1;
    }
    if (neighparam->getFlagXvalid() && neighparam->getType() != ENeigh::MOVING)
    {
      messerr(
          "The Cross-Validation can only be processed with Moving neighborhood");
      return 1;
    }
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Perform the conditional or non-conditional simulation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       Input Db structure (optional)
 ** \param[in]  dbout      Output Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure (optional)
 ** \param[in]  nbsimu     Number of simulations
 ** \param[in]  seed       Seed for random number generator
 ** \param[in]  nbtuba     Number of turning bands
 ** \param[in]  flag_check 1 to check the proximity in Gaussian scale
 ** \param[in]  namconv    Naming convention
 **
 ** \remark  The arguments 'dbout' and 'neigh' are optional: they must
 ** \remark  be defined only for conditional simulations
 **
 *****************************************************************************/
int simtub(Db *dbin,
           Db *dbout,
           Model *model,
           ANeighParam *neighparam,
           int nbsimu,
           int seed,
           int nbtuba,
           int flag_check,
           const NamingConvention& namconv)
{
  SimuTurningBands situba;
  int flag_cond, nvar, error, iext, inostat, iptr_in, iptr_out;

  /* Initializations */

  error = 1;
  nvar = model->getVariableNumber();
  iptr_in = iptr_out = -1;
  flag_cond = (dbin != nullptr);
  if (st_check_simtub_environment(dbin, dbout, model, neighparam)) goto label_end;
  if (manage_external_info(1, ELoc::F, dbin, dbout, &iext)) goto label_end;
  if (manage_external_info(1, ELoc::NOSTAT, dbin, dbout, &inostat))
    goto label_end;

  /* Define the environment variables for printout */

  st_simulation_environment();

  /* Add the attributes for storing the results */

  if (flag_cond)
  {
    if (db_locator_attribute_add(dbin, ELoc::SIMU, nvar * nbsimu, 0, 0.,
                                 &iptr_in)) goto label_end;
  }
  if (db_locator_attribute_add(dbout, ELoc::SIMU, nvar * nbsimu, 0, 0.,
                               &iptr_out)) goto label_end;

  // Processing the Turning Bands algorithm

  situba = SimuTurningBands(nbsimu, nbtuba, model, seed);
  if (situba.simulate(dbin, dbout, neighparam, 0)) goto label_end;

  // Check the simulation at data location

  if (flag_check) situba.checkGaussianData2Grid(dbin, dbout, model);

  /* Free the temporary variables */

  if (flag_cond) dbin->deleteColumnsByLocator(ELoc::SIMU);

  /* Set the error return flag */

  error = 0;
  namconv.setNamesAndLocators(dbin, ELoc::Z, model->getVariableNumber(), dbout,
                              iptr_out, String(), nbsimu);

  label_end:
  (void) manage_external_info(-1, ELoc::F, dbin, dbout, &iext);
  (void) manage_external_info(-1, ELoc::NOSTAT, dbin, dbout, &inostat);
  return (error);
}

/*****************************************************************************/
/*!
 **  Performs the boolean simulation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin          Db structure containing the data (optional)
 ** \param[in]  dbout         DbGrid structure containing the simulated grid
 ** \param[in]  tokens        Tokens structure
 ** \param[in]  boolparam     SimuBooleanParam structure
 ** \param[in]  seed          Seed for the random number generator
 ** \param[in]  flag_simu     Store the boolean simulation
 ** \param[in]  flag_rank     Store the object rank
 ** \param[in]  verbose       1 for a verbose output
 ** \param[in]  namconv       Naming convention
 **
 *****************************************************************************/
int simbool(Db* dbin,
            DbGrid* dbout,
            ModelBoolean* tokens,
            const SimuBooleanParam& boolparam,
            int seed,
            bool flag_simu,
            bool flag_rank,
            bool verbose,
            const NamingConvention& namconv)
{
  int iptr_cover = -1;
  if (dbin != nullptr)
  {
    if (dbin->getVariableNumber() != 1)
    {
      messerr("Conditional Boolean simulation needs 1 variable");
      return 1;
    }
    iptr_cover = dbin->addColumnsByConstant(1, 0.,"Cover",ELoc::Z,1);
    if (iptr_cover < 0) return 1;
  }

  /* Add the attributes for storing the simulation */

  int iptr_simu = -1;
  if (flag_simu)
  {
    iptr_simu = dbout->addColumnsByConstant(1, boolparam.getBackground());
    if (iptr_simu < 0) return 1;
  }
  int iptr_rank = -1;
  if (flag_rank)
  {
    iptr_rank = dbout->addColumnsByConstant(1, TEST);
    if (iptr_rank < 0) return 1;
  }

  SimuBoolean simubool(1, seed);
  if (simubool.simulate(dbin, dbout, tokens, boolparam,
                        iptr_simu, iptr_rank, verbose)) return 1;

  namconv.setNamesAndLocators(dbin, ELoc::Z, 1, dbout, iptr_simu, "Facies", 1,
                              false);
  namconv.setNamesAndLocators(dbin, ELoc::Z, 1, dbout, iptr_rank, "Rank", 1,
                              false);
  return 0;
}

/****************************************************************************/
/*!
 **  Perform the conditional or non-conditional block simulation
 **  in the scope of the Discrete Gaussian Model
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       Input Db structure (optional)
 ** \param[in]  dbout      Output Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure (optional)
 ** \param[in]  rval       Change of support coefficient
 ** \param[in]  seed       Seed for random number generator
 ** \param[in]  nbsimu     Number of simulations
 ** \param[in]  nbtuba     Number of turning bands
 ** \param[in]  flag_check 1 to check the proximity in Gaussian scale
 **
 ** \remark  The arguments 'dbout' and 'neigh' are optional: they must
 ** \remark  be defined only for conditional simulations
 **
 *****************************************************************************/
int simdgm(Db *dbin,
           Db *dbout,
           Model *model,
           ANeighParam *neighparam,
           double rval,
           int seed,
           int nbsimu,
           int nbtuba,
           int flag_check)
{
  SimuTurningBands situba;
  int flag_cond, nvar, error, iext, inostat, iptr;

  /* Initializations */

  error = 1;
  nvar = model->getVariableNumber();
  iptr = -1;
  flag_cond = (dbin != nullptr);
  if (st_check_simtub_environment(dbin, dbout, model, neighparam)) goto label_end;
  if (manage_external_info(1, ELoc::F, dbin, dbout, &iext)) goto label_end;
  if (manage_external_info(1, ELoc::NOSTAT, dbin, dbout, &inostat)) goto label_end;

  /* Define the environment variables for printout */

  st_simulation_environment();
  FLAG_DGM = 1;
  R_COEFF = rval;

  /* Add the attributes for storing the results */

  if (flag_cond)
  {
    if (db_locator_attribute_add(dbin, ELoc::SIMU, nvar * nbsimu, 0, 0., &iptr))
      goto label_end;
  }
  if (db_locator_attribute_add(dbout, ELoc::SIMU, nvar * nbsimu, 0, 0., &iptr))
    goto label_end;

  // Processing the Turning Bands algorithm

  situba = SimuTurningBands(nbsimu, nbtuba, model,seed);
  if (situba.simulate(dbin, dbout, neighparam, 0)) goto label_end;

  // Check the simulations at data locations

  if (flag_check) situba.checkGaussianData2Grid(dbin, dbout, model);

  /* Free the temporary variables */

  if (flag_cond) dbin->deleteColumnsByLocator(ELoc::SIMU);

  /* Set the error return flag */

  error = 0;

  label_end: (void) manage_external_info(-1, ELoc::F, dbin, dbout, &iext);
  (void) manage_external_info(-1, ELoc::NOSTAT, dbin, dbout, &inostat);
  return (error);
}

/****************************************************************************/
/*!
 **  Perform the conditional or non-conditional simulation
 **  with Bayesian Drift
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       Input Db structure (optional)
 ** \param[in]  dbout      Output Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure (optional)
 ** \param[in]  nbsimu     Number of simulations
 ** \param[in]  seed       Seed for random number generator
 ** \param[in]  dmean      Array giving the prior means for the drift terms
 ** \param[in]  dcov       Array containing the prior covariance matrix
 **                        for the drift terms
 ** \param[in]  nbtuba     Number of turning bands
 ** \param[in]  flag_check 1 to check the proximity in Gaussian scale
 ** \param[in]  namconv    Naming convention
 **
 ** \remark  The arguments 'dbout' and 'neigh' are optional: they must
 ** \remark  be defined only for conditional simulations
 **
 *****************************************************************************/
int simbayes(Db *dbin,
             Db *dbout,
             Model *model,
             ANeighParam *neighparam,
             int nbsimu,
             int seed,
             const VectorDouble& dmean,
             const VectorDouble& dcov,
             int nbtuba,
             int flag_check,
             const NamingConvention& namconv)
{
  SimuTurningBands situba;
  int flag_cond, nvar, error, iptr_in, iptr_out;

  /* Initializations */

  error = 1;
  nvar = model->getVariableNumber();
  iptr_in = iptr_out = -1;
  flag_cond = (dbin != nullptr);
  if (st_check_simtub_environment(dbin, dbout, model, neighparam)) goto label_end;

  /* Define the environment variables for printout */

  st_simulation_environment();

  /* Add the attributes for storing the results */

  if (flag_cond)
  {
    if (db_locator_attribute_add(dbin, ELoc::SIMU, nvar * nbsimu, 0, 0, &iptr_in))
      goto label_end;
  }
  if (db_locator_attribute_add(dbout, ELoc::SIMU, nvar * nbsimu, 0, 0, &iptr_out))
    goto label_end;

  // Processing the Turning Bands algorithm

  situba = SimuTurningBands(nbsimu, nbtuba, model, seed);
  if (situba.simulate(dbin, dbout, neighparam, 0, true, dmean, dcov))
    goto label_end;

  // Check simulations at data locations

  if (flag_check) situba.checkGaussianData2Grid(dbin, dbout, model);

  /* Free the temporary variables */

  if (flag_cond) dbin->deleteColumnsByLocator(ELoc::SIMU);

  /* Set the error return flag */

  error = 0;
  namconv.setNamesAndLocators(dbin, ELoc::Z, model->getVariableNumber(), dbout,
                              iptr_out, String(), nbsimu);

  label_end:
  return (error);
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
  if (ipgs <= 0 || propdef == nullptr)
    return (igrf);
  else
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
  return;
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
    facdat = (int) dbin->getVariable(iech, 0);
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
  coor = db_sample_free(coor);
  return;
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
 ** \param[in]  neighparam  Neighborhood structure
 ** \param[in]  nbsimu      Number of simulations
 ** \param[in]  seed        Seed for random number generator
 ** \param[in]  flag_gaus   1 if results must be gaussian; otherwise facies
 ** \param[in]  flag_modif  1 for facies proportion
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
int simpgs(Db *dbin,
           Db *dbout,
           RuleProp *ruleprop,
           Model *model1,
           Model *model2,
           ANeighParam *neighparam,
           int nbsimu,
           int seed,
           int flag_gaus,
           int flag_modif,
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
  SimuTurningBands situba;
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
  iptr_RP = iptr_RF = iptr_DF = iptr_DN = iptr_RN = nfacies = 0;
  iptr = -1;
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
  if (st_check_simtub_environment(dbin, dbout, model1, neighparam)) goto label_end;

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
  if (flag_modif && flag_gaus)
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
    if (model_stabilize(models[igrf], 1, percent)) goto label_end;
    if (model_normalize(models[igrf], 1)) goto label_end;
    modvec.push_back(models[igrf]);
  }

  /* Neighborhood */
  if (flag_cond)
  {
    if (neighparam->getType() != ENeigh::UNIQUE && neighparam->getType()
        != ENeigh::BENCH)
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
  if (flag_modif)
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
    situba = SimuTurningBands(nbsimu, nbtuba, models[igrf], local_seed);
    local_seed = 0;
    if (situba.simulate(dbin, dbout, neighparam, icase, false, VectorDouble(),
                        VectorDouble(), true)) goto label_end;
    if (flag_check) situba.checkGaussianData2Grid(dbin, dbout, models[igrf]);
  }

  /* Convert gaussian to facies at target point */

  if (!flag_gaus)
  {
    for (int isimu = 0; isimu < nbsimu; isimu++)
      simu_func_categorical_transf(dbout, 0, isimu, nbsimu);
  }

  /* Update facies proportions at target points */

  if (flag_modif)
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
    if (!st_keep(flag_gaus, flag_modif, RESULT, TYPE_PROP) && iptr_RP)
      dbout->deleteColumnsByLocator(ELoc::P);
    else
      namconv.setNamesAndLocators(NULL, ELoc::Z, -1, dbout, iptr_RP, "Props",
                                  nbsimu, false);

    if (!st_keep(flag_gaus, flag_modif, RESULT, TYPE_GAUS))
      dbout->deleteColumnsByLocator(ELoc::SIMU);
    else
      namconv.setNamesAndLocators(NULL, ELoc::Z, -1, dbout, iptr_RN, "Gaus",
                                  ngrf * nbsimu, false);

    if (!st_keep(flag_gaus, flag_modif, RESULT, TYPE_FACIES) && iptr_RF)
      dbout->deleteColumnsByLocator(ELoc::FACIES);
    else
      namconv.setNamesAndLocators(NULL, ELoc::Z, -1, dbout, iptr_RF, String(),
                                  nbsimu);
  }

  if (dbin != nullptr)
  {
    if (!st_keep(flag_gaus, flag_modif, DATA, TYPE_GAUS))
      dbin->deleteColumnsByLocator(ELoc::SIMU);
    else
      namconv.setNamesAndLocators(NULL, ELoc::Z, -1, dbin, iptr_DN, "Gaus",
                                  ngrf * nbsimu, false);

    if (!st_keep(flag_gaus, flag_modif, DATA, TYPE_FACIES) && iptr_DF)
      dbin->deleteColumnsByLocator(ELoc::FACIES);
    else
      namconv.setNamesAndLocators(NULL, ELoc::Z, -1, dbin, iptr_DF, String(),
                                  nbsimu, false);

    dbin->deleteColumnsByLocator(ELoc::GAUSFAC);
    dbin->deleteColumnsByLocator(ELoc::L);
    dbin->deleteColumnsByLocator(ELoc::U);
  }

  /* Set the error return flag */

  error = 0;

  label_end:
    propdef = proportion_manage(-1, 1, flag_stat, ngrf, 0, nfacies, 0,
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
 ** \param[in]  neighparam  Neighborhood structure
 ** \param[in]  nbsimu      Number of simulations
 ** \param[in]  seed        Seed for random number generator
 ** \param[in]  flag_gaus   1 gaussian results; otherwise facies
 ** \param[in]  flag_modif  1 for facies proportion
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
             ANeighParam *neighparam,
             int nbsimu,
             int seed,
             int flag_gaus,
             int flag_modif,
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
  SimuTurningBands situba;
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
    iatt_z[0] = db_attribute_identify(dbin, ELoc::Z, 0);
    iatt_z[1] = db_attribute_identify(dbin, ELoc::Z, 1);
  }

  /* Output Db */
  if (flag_modif && flag_gaus)
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
      if (model_stabilize(models[ipgs][igrf], 1, percent)) goto label_end;
      if (model_normalize(models[ipgs][igrf], 1)) goto label_end;

      modvec[ipgs].push_back(models[ipgs][igrf]);
    }
  }

  /* Neighborhood */
  if (neighparam->getType() != ENeigh::UNIQUE && neighparam->getType() != ENeigh::BENCH)
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
    if (st_check_simtub_environment(dbin, dbout, models[ipgs][0], neighparam))
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
  if (flag_modif)
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
      situba = SimuTurningBands(nbsimu, nbtuba, models[ipgs][igrf], local_seed);
      local_seed = 0;
      if (situba.simulate(dbin, dbout, neighparam, icase, false, VectorDouble(),
                          VectorDouble(), true)) goto label_end;
      if (flag_check) situba.checkGaussianData2Grid(dbin, dbout, models[ipgs][igrf]);
    }

    /* Convert gaussian to facies at target point */

    if (!flag_gaus) for (int isimu = 0; isimu < nbsimu; isimu++)
      simu_func_categorical_transf(dbout, 0, isimu, nbsimu);

    /* Update facies proportions at target points */

    if (flag_modif)
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
    if (!st_keep(flag_gaus, flag_modif, RESULT, TYPE_PROP) && iptr_RP)
      dbout->deleteColumnsByLocator(ELoc::P);
    else
      namconv.setNamesAndLocators(NULL, ELoc::Z, -1, dbout, iptr_RP, "Props",
                                  nfactot, false);

    if (!st_keep(flag_gaus, flag_modif, RESULT, TYPE_GAUS) && iptr_RN)
      dbout->deleteColumnsByLocator(ELoc::SIMU);
    else
      namconv.setNamesAndLocators(NULL, ELoc::Z, -1, dbout, iptr_RN, "Gaus",
                                  ngrftot * nbsimu, false);

    if (!st_keep(flag_gaus, flag_modif, RESULT, TYPE_FACIES) && iptr_RF)
      dbout->deleteColumnsByLocator(ELoc::FACIES);
    else
      namconv.setNamesAndLocators(NULL, ELoc::Z, -1, dbout, iptr_RF, String(),
                                  npgs * nbsimu);
  }

  if (dbin != nullptr)
  {
    if (!st_keep(flag_gaus, flag_modif, DATA, TYPE_GAUS) && iptr_DN)
      dbin->deleteColumnsByLocator(ELoc::SIMU);
    else
      namconv.setNamesAndLocators(NULL, ELoc::Z, -1, dbin, iptr_DN, "Gaus",
                                  ngrftot * nbsimu, false);

    if (!st_keep(flag_gaus, flag_modif, DATA, TYPE_FACIES) && iptr_DF)
      dbin->deleteColumnsByLocator(ELoc::FACIES);
    else
      namconv.setNamesAndLocators(NULL, ELoc::Z, -1, dbin, iptr_DF, String(),
                                  npgs * nbsimu, false);

    dbin->deleteColumnsByLocator(ELoc::GAUSFAC);
    dbin->deleteColumnsByLocator(ELoc::L);
    dbin->deleteColumnsByLocator(ELoc::U);
  }

  /* Set the error return flag */

  error = 0;

label_end:
  st_suppress_added_samples(dbin,nechin);
  propdef = proportion_manage(-1, 1, flag_stat, ngrf[0], ngrf[1], nfac[0],
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
                         const ELoc &locatorType,
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
        db->updArray(iech, iptr_ce + ivar, 0, value);
        db->updArray(iech, iptr_cstd + ivar, 0, value * value);
        db->updArray(iech, iptr_nb + ivar, 0, 1.);
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

  label_end: (void) db_attribute_del_mult(db, iptr_nb, nvar);
  iptr_nb = -1;
  if (error)
  {
    (void) db_attribute_del_mult(db, iptr_ce, nvar);
    (void) db_attribute_del_mult(db, iptr_cstd, nvar);
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
 ** \param[in]  neighparam  Neigh structure (optional)
 ** \param[in]  nbsimu      Number of simulations
 ** \param[in]  seed        Seed for random number generator
 ** \param[in]  gibbs_nburn Initial number of iterations for bootstrapping
 ** \param[in]  gibbs_niter Maximum number of iterations
 ** \param[in]  flag_norm   1 if the Model must be normalized
 ** \param[in]  flag_multi_mono  1 for the Multi_mono algorithm
 ** \param[in]  flag_propagation 1 for the propagation algorithm
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
                  ANeighParam *neighparam,
                  int nbsimu,
                  int seed,
                  int gibbs_nburn,
                  int gibbs_niter,
                  bool flag_norm,
                  bool flag_multi_mono,
                  bool flag_propagation,
                  bool /*flag_sym_neigh*/,
                  int gibbs_optstats,
                  double percent,
                  bool flag_ce,
                  bool flag_cstd,
                  bool verbose,
                  const NamingConvention& namconv)
{
  int error, iptr, npgs, nvar, iptr_ce, iptr_cstd;
  PropDef *propdef;

  /* Initializations */

  error = 1;
  npgs = 1;
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
    if (model_stabilize(model, 1, percent)) goto label_end;
  }
  if (flag_norm)
  {
    if (model_normalize(model, 1)) goto label_end;
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
      gibbs = GibbsFactory::createGibbs(dbin, model, neighparam);
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
      if (gibbs->run(y, 0, isimu, verbose)) goto label_end;
  }

  /* Convert the simulations */

  if (flag_ce || flag_cstd)
  {
    if (db_simulations_to_ce(dbin, ELoc::GAUSFAC, nbsimu, nvar, &iptr_ce,
                             &iptr_cstd)) goto label_end;

    // We release the attributes dedicated to simulations on Dbout

    if (!flag_ce)
    {
      (void) db_attribute_del_mult(dbin, iptr_ce, nvar);
      iptr_ce = -1;
    }
    if (!flag_cstd)
    {
      (void) db_attribute_del_mult(dbin, iptr_cstd, nvar);
      iptr_cstd = -1;
    }
    dbin->deleteColumnsByLocator(ELoc::GAUSFAC);
  }

  /* Set the error return flag */

  error = 0;
  namconv.setNamesAndLocators(dbin, ELoc::UNKNOWN, nvar, dbin, iptr, String(),
                              nbsimu);

  label_end: propdef = proportion_manage(-1, 0, 1, 1, 0,
                                         model->getVariableNumber(), 0, dbin,
                                         NULL,
                                         VectorDouble(), propdef);
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
 ** \param[in]  neighparam   ANeighParam structure (optional)
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
int simtub_constraints(Db *dbin,
                       Db *dbout,
                       Model *model,
                       ANeighParam *neighparam,
                       int seed,
                       int nbtuba,
                       int nbsimu_min,
                       int nbsimu_quant,
                       int niter_max,
                       VectorInt &cols,
                       int (*func_valid)(int flag_grid,
                                         int nDim,
                                         int nech,
                                         int *nx,
                                         double *dx,
                                         double *x0,
                                         double nonval,
                                         double percent,
                                         VectorDouble &tab))
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

  flag_grid = is_grid(dbout);
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
    if (simtub(dbin, dbout, model, neighparam, nbsimu, 0, nbtuba, 0)) goto label_end;

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

  label_end: nx = (int*) mem_free((char* ) nx);
  dx = (double*) mem_free((char* ) dx);
  x0 = (double*) mem_free((char* ) x0);
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
  return;
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
  SimuTurningBands situba;
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
  niter = nleft = last = 0;
  while (1)
  {
    niter++;
    tpois -= log(law_uniform(0., 1.));

    /* Mask the nodes that cannot be accessed anymore */

    nleft = st_maxstable_mask(dbout, seuil, tpois, iptrv, iptrs);
    if (nleft <= 0) break;

    /* Processing the Turning Bands algorithm */

    situba = SimuTurningBands(1, nbtuba, model, seed);
    if (situba.simulate(nullptr, dbout, nullptr, 0)) goto label_end;

    /* Combine the newly simulated outcome to the background */

    st_maxstable_combine(dbout, tpois, niter, iptrg, iptrv, iptrr, &last);

    if (verbose)
      message("Iteration #%2d - Scale = %6.3lf - Nb. cells left = %d\n", niter,
              tpois, nleft);

    /* Update the model for next iteration */

    for (icov = 0; icov < model->getCovaNumber(); icov++)
      model->getCova(icov)->setRange(model->getCova(icov)->getRange() * ratio);
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
  SimuTurningBands situba;
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
    dbout->setSelection(iech, 1);

  /* Loop on the cutoff classes */

  total = 0;
  for (icut = 0; icut < ncut; icut++)
  {

    /* Simulation in the non-masked part of the grid */

    situba = SimuTurningBands(1, nbtuba, model, seed);
    if (situba.simulate(nullptr, dbout, nullptr, 0)) goto label_end;

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
      dbout->setSelection(iech, 0);
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

  label_end: sort = (double*) mem_free((char* ) sort);
  pton = (double*) mem_free((char* ) pton);
  pres = (double*) mem_free((char* ) pres);
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
 ** \param[in]  flag_modif  1 for facies proportion
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
int simpgs_spde(Db *dbin,
                Db *dbout,
                RuleProp *ruleprop,
                Model *model1,
                Model *model2,
                const String &triswitch,
                const VectorDouble &gext,
                int flag_gaus,
                int flag_modif,
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
  if (flag_modif && flag_gaus)
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
    if (model_stabilize(models[igrf], 1, percent)) goto label_end;
    if (model_normalize(models[igrf], 1)) goto label_end;
  }
  if (spde_check(dbin, dbout, model1, model2, verbose, gext, 1, 1, 1, 0, 0, 1,
                 flag_modif)) goto label_end;
  s_option = spde_option_alloc();
  spde_option_update(s_option, triswitch);

  /* Define the environment variables for printout */

  /**********************/
  /* Add the attributes */
  /**********************/

  /* Storage of the facies proportions */
  if (flag_modif)
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

  if (!st_keep(flag_gaus, flag_modif, RESULT, TYPE_PROP) && iptr_RP)
    dbout->deleteColumnsByLocator(ELoc::P);

  if (!st_keep(flag_gaus, flag_modif, RESULT, TYPE_FACIES) && iptr_RF)
    dbout->deleteColumnsByLocator(ELoc::FACIES);

  if (!st_keep(flag_gaus, flag_modif, RESULT, TYPE_GAUS))
    dbout->deleteColumnsByLocator(ELoc::SIMU);

  dbin->deleteColumnsByLocator(ELoc::L);
  dbin->deleteColumnsByLocator(ELoc::U);

  /* Set the error return flag */

  error = 0;

  label_end: propdef = proportion_manage(-1, 1, flag_stat, ngrf, 0, nfacies, 0,
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
  SimuTurningBands situba;
  PropDef *propdef;
  ANeighParam *neighparam = nullptr;
  int nvar, error, iext, inostat, iptr, iptr_ce, iptr_cstd, ndim;

  /* Initializations */

  error = 1;
  nvar = model->getVariableNumber();
  ndim = model->getDimensionNumber();
  iptr = -1;
  propdef = nullptr;

  /* Preliminary checks */

  neighparam = NeighUnique::create(ndim,false);
  law_set_random_seed(seed);
  if (st_check_simtub_environment(dbin, dbout, model, NULL)) goto label_end;
  if (manage_external_info(1, ELoc::F, dbin, dbout, &iext)) goto label_end;
  if (manage_external_info(1, ELoc::NOSTAT, dbin, dbout, &inostat))
    goto label_end;

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

  situba = SimuTurningBands(nbsimu, nbtuba, model, seed);
  if (situba.simulate(dbin, dbout, neighparam, 0, false, VectorDouble(),
                      VectorDouble(), false, true)) goto label_end;
  if (flag_check) situba.checkGaussianData2Grid(dbin, dbout, model);

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
      (void) db_attribute_del_mult(dbout, iptr_ce, nvar);
      iptr_ce = -1;
    }
    if (!flag_cstd)
    {
      (void) db_attribute_del_mult(dbout, iptr_cstd, nvar);
      iptr_cstd = -1;
    }
  }

  /* Set the error return flag */

  error = 0;

  label_end:
  delete neighparam;
  (void) manage_external_info(-1, ELoc::F, dbin, dbout, &iext);
  (void) manage_external_info(-1, ELoc::NOSTAT, dbin, dbout, &inostat);
  return (error);
}
