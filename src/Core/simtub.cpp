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
#include "../../include/Gibbs/GibbsMultiMono.hpp"
#include "../../include/Gibbs/GibbsUPropMono.hpp"
#include "Gibbs/GibbsFactory.hpp"
#include "Morpho/Morpho.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Covariances/CovAniso.hpp"
#include "Basic/MathFunc.hpp"
#include "LithoRule/PropDef.hpp"
#include "LithoRule/Rule.hpp"
#include "LithoRule/RuleShift.hpp"
#include "LithoRule/RuleShadow.hpp"
#include "geoslib_e.h"

/*! \cond */
#define DATA   0
#define RESULT 1

#define AD(ivar,jvar)           (ivar) + nvar * (jvar)
#define SEEDS(ivar,is,ib,isimu)  situba->seeds[(ivar)+nvar*((is)+ncova*((ib)+nbtuba*(isimu)))]
#define AIC(icov,ivar,jvar)      aic[(icov)*nvar*nvar + AD(ivar,jvar)]
#define IND(iech,ivar)           ((iech) + (ivar) * nech)

#define GAUS   0
#define FACIES 1
#define PROP   2

#define EPS 1e-5

/*! \endcond */

static int QUANT = 1000;
static int MES_IGRF,MES_NGRF,MES_IPGS,MES_NPGS,FLAG_DGM;
static double GIBBS_RHO,GIBBS_SQR;
static double R_COEFF;
static Modif_Categorical ModCat = {0,{0,0},NULL,NULL};

/****************************************************************************/
/*!
**  Initialize the global values
**
*****************************************************************************/
static void st_simulation_environment(void)
{
  MES_IGRF = 0;
  MES_NGRF = 0;
  MES_IPGS = 0;
  MES_NPGS = 0;
  FLAG_DGM = 0;

  GIBBS_RHO = 0.;
  GIBBS_SQR = 0.;
  R_COEFF   = 0.;

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
static int st_facies(PropDef *propdef,
                     int    ipgs,
                     int    ifac)
{
  if (ipgs <= 0)
    return(ifac);
  else
    return(propdef->nfac[0] + ifac);
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
GEOSLIB_API void simu_func_categorical_transf(Db  *db,
                                              int  verbose,
                                              int  isimu,
                                              int  nbsimu)
{
  const Rule* rule = ModCat.rule;

  rule->gaus2facResult(ModCat.propdef, db, ModCat.flag_used, ModCat.ipgs,
                       isimu, nbsimu);

  /* Optional printout */

  if (verbose) 
    message("Simulation Categorical Transformation (%d/%d)\n",isimu+1,nbsimu);
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
GEOSLIB_API void simu_func_continuous_update(Db  *db,
                                             int  verbose,
                                             int  isimu,
                                             int  nbsimu)
{
  int    iptr_simu;
  double simval;

  /* Preliminary checks */

  check_mandatory_attribute("simu_func_continuous_update",db,LOC_SIMU);
  check_mandatory_attribute("simu_func_continuous_update",db,LOC_Z);
  iptr_simu = db->getSimvarRank(isimu,0,0,nbsimu,1);

  /* Loop on the grid cells */

  for (int iech=0; iech<db->getSampleNumber(); iech++)
  {
    if (! db->isActive(iech)) continue;
    simval = get_LOCATOR_ITEM(db,LOC_SIMU,iptr_simu,iech);
    db->updVariable(iech,0,0,simval);
    db->updVariable(iech,1,0,simval * simval);
  }

  /* Optional printout */

  if (verbose) 
    message("Simulation Continuous Update (%d/%d)\n",isimu+1,nbsimu);
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
GEOSLIB_API void simu_func_categorical_update(Db  *db,
                                              int  verbose,
                                              int  isimu,
                                              int  nbsimu)
{
  int iptr_simu,facies,rank,ipgs;
  double prop;

  /* Preliminary checks */

  ipgs = ModCat.ipgs;
  check_mandatory_attribute("simu_func_categorical_update",db,LOC_FACIES);
  check_mandatory_attribute("simu_func_categorical_update",db,LOC_P);
  iptr_simu = db->getSimvarRank(isimu,0,ipgs,nbsimu,1);

  /* Loop on the grid cells */

  for (int iech=0; iech<db->getSampleNumber(); iech++)
  {
    if (! db->isActive(iech)) continue;
    facies = (int) get_LOCATOR_ITEM(db,LOC_FACIES,iptr_simu,iech) - 1;
    rank   = st_facies(ModCat.propdef,ipgs,facies);
    prop   = db->getProportion(iech,rank) + 1.;
    db->setProportion(iech,rank,prop);
  }

  /* Optional printout */

  if (verbose) 
    message("Simulation Categorical Update (%d/%d)\n",isimu+1,nbsimu);
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
GEOSLIB_API void simu_func_continuous_scale(Db  *db,
                                            int  verbose,
                                            int  nbsimu)
{
  double mean,stdv;

  /* Preliminary checks */

  check_mandatory_attribute("simu_func_continuous_scale",db,LOC_Z);

  /* Loop on the grid cells */

  for (int iech=0; iech<db->getSampleNumber(); iech++)
  {
    if (! db->isActive(iech)) continue;
    mean = db->getVariable(iech,0) / nbsimu;
    db->setVariable(iech,0,mean);
    stdv = db->getVariable(iech,1) / nbsimu - mean * mean;
    stdv = (stdv > 0) ? sqrt(stdv) : 0.;
    db->setVariable(iech,1,stdv);
  }

  /* Optional printout */

  if (verbose) 
    message("Simulation Continuous Scaling (%d)\n",nbsimu);
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
GEOSLIB_API void simu_func_categorical_scale(Db  *db,
                                             int  verbose,
                                             int  nbsimu)
{
  int    rank,nfacies,ipgs;
  double prop;
  PropDef *propdef;

  /* Preliminary checks */

  propdef = ModCat.propdef;
  ipgs    = ModCat.ipgs;
  nfacies = propdef->nfac[ipgs];
  check_mandatory_attribute("simu_func_categorical_scale",db,LOC_P);

  /* Loop on the grid cells */

  for (int iech=0; iech<db->getSampleNumber(); iech++)
  {
    if (! db->isActive(iech)) continue;
    for (int ifac=0; ifac<nfacies; ifac++)
    {
      rank = st_facies(propdef,ipgs,ifac);
      prop = db->getProportion(iech,rank) / (double) nbsimu;
      db->setProportion(iech,rank,prop);
    }
  }

  /* Optional printout */

  if (verbose) 
    message("Simulation Categorical Scaling (%d)\n",nbsimu);
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
GEOSLIB_API void check_mandatory_attribute(const char *method,
                                           Db* db,
                                           ENUM_LOCS locatorType)
{
  if (get_LOCATOR_NITEM(db,locatorType) <= 0)
    messageAbort("%s : Attributes %d are mandatory",method,locatorType);
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
static int st_keep(int flag_gaus,
                   int flag_modif,
                   int file,
                   int type)
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
      case 0:                        /* Gaussian */
        keep = (flag_gaus && ! flag_modif);
        break;
        
      case 1:                        /* Facies */
        keep = (! flag_gaus && ! flag_modif);
        break;
        
      case 2:                        /* Proportion */
        keep = (flag_modif);
        break;
    }
  }

label_end:
  return(keep);
}

/****************************************************************************/
/*!
**  Checks the environment for simulations by Turning Bands
**
** \return  Error return code
**
** \param[in]  dbin   input Db structure (optional if non conditional)
** \param[in]  dbout  output Db structure
** \param[in]  model  Model structure
** \param[in]  neigh  Neigh structure (optional if non conditional)
**
*****************************************************************************/
static int st_check_simtub_environment(Db    *dbin,
                                       Db    *dbout,
                                       Model *model,
                                       Neigh *neigh)
{
  double *dbin_mini,*dbin_maxi,*dbout_mini,*dbout_maxi;
  int error,ndim,nvar,nfex,flag_cond;
  /* Initializations */

  error = 1;
  nvar = ndim = nfex = 0;
  dbin_mini = dbin_maxi = dbout_mini = dbout_maxi = (double *) NULL;
  flag_cond = (dbin != (Db *) NULL) && (neigh != (Neigh *) NULL);
  ndim = dbout->getNDim();

  /**************************************************************/
  /* Check if the Space dimension is compatible with the method */
  /**************************************************************/

  if (ndim > 3)
  {
    messerr("The Turning Band Method is not a relevant simulation model");
    messerr("for this Space Dimension (%d)",ndim);
    goto label_end;
  }

  /*********************************/
  /* Compatibility between two Dbs */
  /*********************************/

  if (flag_cond && ! dbin->hasSameDimension(dbout)) goto label_end;

  /**********************/
  /* Checking the model */
  /**********************/

  if (model != (Model *) NULL)
  {
    nvar = model->getVariableNumber();
    if (nvar <= 0)
    {
      messerr("The number of variables must be positive = %d",model->getVariableNumber());
      goto label_end;
    }
    if (flag_cond && dbin->getVariableNumber() != nvar)
    {
      messerr("The number of variables of the Data (%d)",dbin->getVariableNumber());
      messerr("does not match the number of variables of the Model (%d)",
              nvar);
      goto label_end;
    }
    if (model->getCovaNumber() <= 0)
    {
      messerr("The number of covariance must be positive");
      goto label_end;
    }

    if (model->getDimensionNumber() <= 0)
    {
      messerr("The Space Dimension must be positive = %d",model->getDimensionNumber());
      goto label_end;
    }
    if (model->getDimensionNumber() != ndim)
    {
      messerr("The Space Dimension of the Db structure (%d)",ndim);
      messerr("Does not correspond to the Space Dimension of the model (%d)",
              model->getDimensionNumber());
      goto label_end;
    }

    nfex = model_nfex(model);
    if (flag_cond && nfex != 0 && 
        ! is_grid(dbout) && dbin->getExternalDriftNumber() != nfex)
    {
      messerr("The Model requires %d external drift(s)",
              model_nfex(model));
      messerr("but the input Db refers to %d external drift variables",
              dbin->getExternalDriftNumber());
      goto label_end;
    }
    if (nfex != 0 && dbout->getExternalDriftNumber() != nfex)
    {
      messerr("The Model requires %d external drift(s)",
              model_nfex(model));
      messerr("but the output Db refers to %d external drift variables",
              dbout->getExternalDriftNumber());
      goto label_end;
    }
  }

  /*********************************/
  /* Calculate the field extension */
  /*********************************/

  /* Input Db structure */

  if (flag_cond)
  {
    dbin_mini = db_sample_alloc(dbin,LOC_X);
    if (dbin_mini == (double *) NULL) goto label_end;
    dbin_maxi = db_sample_alloc(dbin,LOC_X);
    if (dbin_maxi == (double *) NULL) goto label_end;
    if (db_extension(dbin,
                     dbin_mini,dbin_maxi,(double *) NULL)) goto label_end;
  }

  /* Output Db structure */

  dbout_mini = db_sample_alloc(dbout,LOC_X);
  if (dbout_mini == (double *) NULL) goto label_end;
  dbout_maxi = db_sample_alloc(dbout,LOC_X);
  if (dbout_maxi == (double *) NULL) goto label_end;
  if (db_extension(dbout,
                   dbout_mini,dbout_maxi,(double *) NULL)) goto label_end;

  if (model != (Model *) NULL)
    model->setField(ut_merge_extension(ndim,dbin_mini,dbin_maxi,
                                      dbout_mini,dbout_maxi));

  dbin_mini  = db_sample_free(dbin_mini);
  dbin_maxi  = db_sample_free(dbin_maxi);
  dbout_mini = db_sample_free(dbout_mini);
  dbout_maxi = db_sample_free(dbout_maxi);

  /*****************************/
  /* Checking the Neighborhood */
  /*****************************/

  if (flag_cond && neigh != (Neigh *) NULL)
  {
    if (neigh->getNDim() != ndim)
    {
      messerr("The Space Dimension of the Neighborhood (%d)",neigh->getNDim());
      messerr("does not correspond to the Space Dimension of the first Db (%d)",
              ndim);
      goto label_end;
    }
    if (neigh->getFlagXvalid() && neigh->getType() != ENeigh::MOVING)
    {
      messerr("The Cross-Validation can only be processed with Moving neighborhood");
      goto label_end;
    }
  }

  /* Set the error return code */

  error = 0;

label_end:
  return(error);
}

/****************************************************************************/
/*!
**  Calculate the Poisson intensity for the generation
**  of the Wiener-Levy along the line
**
** \param[in]  situba Situba structure
**
** \remark  The average number of points per band is calculated:
** \remark  - so as too have in average one band between two target
** \remark    points (either data or target)
** \remark  - to have a number of Poisson points per band lying within
** \remark    [nmini; nmaxi]
**
*****************************************************************************/
static void st_density(Situba *situba)
  
{
  double scale;
  int nbtuba,naverage,npoint;
  int nmini  = 5;
  int nmaxi  = 5000;

  nbtuba   = situba->nbtuba;
  npoint   = situba->nb_points_simu;
  naverage = npoint / nbtuba;
  if (naverage < nmini) naverage = nmini;
  if (naverage > nmaxi) naverage = nmaxi;
  scale = situba->field / naverage;
  situba->theta = 1. / scale;
  return;
}

/*****************************************************************************/
/*!
**  Generates the process constituted by independent gaussian
**  variables along a 1D Poisson process. The process consists in
**  the integration(s) of the previous process Perform the core
**  allocation
**
** \return  The array of gaussian values
**
** \param[in]  nt     number for exponential intervals
** \param[in]  type   degree of the IRF -> number of integrations
**
** \param[out] v0     Wiener-Levy process
** \param[out] v1     First integration of the Wiener-Levy process
** \param[out] v2     Second integration of the Wiener-Levy process
**
** \remark  This procedure allocates memory that should be freed
**
*****************************************************************************/
static void st_irf_process_alloc(int nt,
                                 int type,
                                 double **v0,
                                 double **v1,
                                 double **v2)
{
  int i;

  /* Initializations */

  (*v0) = (*v1) = (*v2) = (double *) NULL;

  /* Generations of the independent gaussian variables */

  for (i=0; i<nt; i++) (void) law_gaussian();

  /* Core allocation */

  (*v0) = (double *) mem_alloc(sizeof(double) * nt,0);
  if ((*v0) == (double *) NULL) return;
  if (type == COV_LINEAR || type == COV_ORDER1_GC) return;

  (*v1) = (double *) mem_alloc(sizeof(double) * nt,0);
  if ((*v1) == (double *) NULL) return;
  if (type == COV_ORDER3_GC) return;

  (*v2) = (double *) mem_alloc(sizeof(double) * nt,0);
  if ((*v2) == (double *) NULL) return;
  if (type == COV_ORDER5_GC) return;

  return;
}

/*****************************************************************************/
/*!
**  Generates the process constituted by independent gaussian
**  variables along a 1D Poisson process. The process consists in
**  the integration(s) of the previous process
**
** \return  The array of gaussian values
**
** \param[in]  nt     number for exponential intervals
** \param[in]  type   degree of the IRF -> number of integrations
** \param[in]  t      Array giving the coordinates of the Poisson points
**
** \param[out] v0     Wiener-Levy process
** \param[out] v1     First integration of the Wiener-Levy process
** \param[out] v2     Second integration of the Wiener-Levy process
**
*****************************************************************************/
static void st_irf_process(int     nt,
                           int     type,
                           double *t,
                           double *v0,
                           double *v1,
                           double *v2)
{
  double delta;
  int    i;

  /* Generation of the Wiener-Levy process */

  v0[0] = 0.;
  for (i=1; i<nt; i++)
    v0[i] = v0[i-1] + law_gaussian();
  if (type == COV_LINEAR || type == COV_ORDER1_GC) return;

  /* First integration of the Wiener-Levy process */

  v1[0] = 0.;
  for (i=1; i<nt; i++)
  {
    delta = t[i] - t[i-1];
    v1[i] = v1[i-1] + v0[i-1] * delta;
  }
  if (type == COV_ORDER3_GC) return;

  /* Second integration of the Wiener-Levy process */

  v2[0] = 0.;
  for (i=1; i<nt; i++)
  {
    delta = t[i] - t[i-1];
    v2[i] = v2[i-1] + v1[i-1] * delta + v0[i-1] * delta * delta / 2.;
  }
  if (type == COV_ORDER5_GC) return;

  return;
}

/*****************************************************************************/
/*!
**  Sample the Wiener-Levy (integrated) process
**
** \param[in]  type   type of polynomial generalized covariance
** \param[in]  nt0    Rank of the Poisson point
** \param[in]  t0     starting time 
** \param[in]  t      Poisson point process 
** \param[in]  v0     Wiener-Levy process
** \param[in]  v1     First integration of the Wiener-Levy process
** \param[in]  v2     Second integration of the Wiener-Levy process
**
*****************************************************************************/
static double st_irf_process_sample(int     type,
                                    int     nt0,
                                    double  t0,
                                    double *t,
                                    double *v0,
                                    double *v1,
                                    double *v2)
{
  double value,delta;

  /* Initializations */

  delta = t0 - t[nt0];

  /* Wiener-Levy process */

  value = v0[nt0];
  if (type == COV_LINEAR || type == COV_ORDER1_GC) return(value);

  /* First integration of the Wiener-Levy process */

  value = v1[nt0] + v0[nt0] * delta;
  if (type == COV_ORDER3_GC) return(value);

  /* Second integration of the Wiener-Levy process */

  value = v2[nt0] + v1[nt0] * delta + v0[nt0] * delta * delta / 2.;
  if (type == COV_ORDER5_GC) return(value);

  return(TEST);
}

/****************************************************************************/
/*!
**  Calculate the correction factor for IRF_k models
**
** \return  Correction factor
**
** \param[in]  type    type of polynomial generalized covariance
** \param[in]  theta1  Equal to inverse of theta value
** \param[in]  scale   Range of the model
**
*****************************************************************************/
static double st_irf_correc(int    type,
                            double theta1,
                            double scale)
{
  double correc;

  /* Initializations */

  correc = 1.;
  
  /* Dispatch */

  switch (type)
  {
    case COV_LINEAR:
    case COV_ORDER1_GC:
      correc = sqrt((4. * theta1) / scale);
      break;

    case COV_ORDER3_GC:
      correc = sqrt((48. * theta1) / scale) / scale;
      break;

    case COV_ORDER5_GC:
      correc = sqrt((1440. * theta1) / scale) / scale / scale;
      break;
  }

  return(correc);
}
/****************************************************************************/
/*!
**  Calculate the scale for 1D process for the K-Bessel model (param<0.5) 
**
** \return  Scale parameter of the 1D process to simulate (param<0.5)
**
** \param[in]  param       Third parameter of the K-Bessel covariance model
** \param[in]  scale       Scale parameter of the model
**
*****************************************************************************/
static double st_compute_scale_Kb(double param,
                                  double scale)
                         
{
  double scale_out;

  scale_out = scale * sqrt(law_beta1(param,0.5 - param));

  return(scale_out);
}
/****************************************************************************/
/*!
**  Calculate the scale for 1D process for the stable model 
**
** \return  Scale parameter of the 1D process to simulate 
**
** \param[in]  alpha       Third parameter of the stable covariance model
** \param[in]  scale       Scale parameter of the model
**
*****************************************************************************/
static double st_compute_scale(double alpha,
                               double scale)
                         
{
  double scale_out;
  
  if(alpha<1)
    scale_out = scale / law_stable_standard_abgd(alpha);
  else
    scale_out = scale / sqrt(law_stable_standard_abgd(alpha/2.));
			     
  return(scale_out);
}
/*****************************************************************************/
/*!
**  Prepare a spectral method
**
** \param[in]  type   type of method for generating the period
** \param[in]  scale  scale factor
** \param[in]  param  factor for the period
**
** \param[out] omega  resulting period
** \param[out] phi    uniform phase lying within [0,2 PI]
**
*****************************************************************************/
static void st_spectral(int     type,
                        double  scale,
                        double  param,
                        double *omega,
                        double *phi)
{
  double val,period;
  int    i;

  period = 0.;
  switch (type)
  {
    case COV_GAUSSIAN:
      for (i=0; i<3; i++)
      {
        val = law_gaussian();
        period += val * val;
      }
      period = sqrt(param * period) / scale;
      break;
	  
    case COV_STABLE:
      for (i=0; i<3; i++)
      {
        val = law_gaussian();
        period += val * val;
      }
      scale  = st_compute_scale(param,scale);
      period = sqrt(2. * period) / scale;
      break;
	  
    case COV_SINCARD:
      val = (law_uniform(0.,1.) >= 0.5) ? 1 : -1;
      period = val / scale;
      break;
      
    case COV_BESSEL_J:
      val = law_beta1(1.5,param - 0.5);
      period = sqrt(val) / scale;
      break;

    case COV_BESSEL_K:
      param =  sqrt(2. * law_gamma(param));
      for (i=0; i<3; i++)
      {
        val = law_gaussian();
        period += val * val;
      }
      period = sqrt(period) / (param * scale);
      break;
	  
    default:
      break;
  }

  *omega = period;
  *phi   = 2. * GV_PI * law_uniform(0.,1.);
  
  return;
} 

/*****************************************************************************/
/*!
**  Generate the 1D stochastic process  
**
** \param[in]  ib      Rank of the turning band
** \param[in]  scale   scale factor
** \param[in]  alpha   power of the variogram h^alpha
**
** \param[out] omega   period = 2piR 
**                     where R is a second kind beta variable with parameters 
**                     1-alpha/2 and alpha/2
** \param[out] phi     uniform phase lying within [0,2 PI]
** \param[out] theta_3 value of theta_alpha,3(R) 
** \param[out] correc0 value to substract from Y_alpha,2 in order to avoid 
**                     numerical problems when R is very small
**
**  \remark Y_alpha_1=theta_alpha_1(R)cos(2pi R.x+phi) 
**  \remark used to simulate a GRF with a power semi-variogram h^alpha
**  \remark according to the method proposed in 
**  \remark Emery, X. and Lantuejoul, C. (2008) 
**  \remark A spectral approach to simulating intrinsec random fields with power
**  \remark and spline generalized covariance. 
**  \remark In Computational Geosciences 12:121-132
**
*****************************************************************************/
static void st_set_power_1D(int     ib,
                            double  scale,
                            double  alpha,
                            double *omega,
                            double *phi,
                            double *theta_3,
                            double *correc0)
{  
  double R,theta_1;
  static double twoPI,log3s2,log1s2,logap1,logap1s2,logap3s2,as2,coeff,coeff3;
  static double alpha_mem = -1.;

  if (ib == 0 || alpha != alpha_mem)
  {
    as2       = alpha / 2.;
    twoPI     = 2. * GV_PI;
    log3s2    = loggamma(1.5);
    log1s2    = loggamma(0.5);
    logap1    = loggamma(1.0 + alpha);
    logap1s2  = loggamma(0.5 + as2);
    logap3s2  = loggamma(1.5 + as2);
    coeff     = 2. * sqrt(exp(logap1) / pow(twoPI,alpha)) * pow(scale,as2);
    coeff3    = sqrt(exp(log3s2 + logap1s2 - log1s2 - logap3s2));
    alpha_mem = alpha;
  }

  R        = law_beta2(1-as2,as2);
  theta_1  = coeff * sqrt((R+1.) / pow(R,as2+1));
  *correc0 = cos(*phi);
  *omega   = twoPI * R;
  *phi     = twoPI * law_uniform(0.,1.);
  *theta_3 = theta_1 / coeff3;
}

/*****************************************************************************/
/*!
**  Generate the 1D stochastic process for spline covariance
**
** \param[in]  ib      Rank of the turning band
** \param[in]  scale   scale factor
** \param[in]  k       power of the variogram h^(2k)log(h)
**
** \param[out] omega   period = 2piR 
**                     where R is a second kind beta variable with parameters 
**                     1/2 and 1/2
** \param[out] phi     uniform phase lying within [0,2 PI]
** \param[out] xi_3    value of xi_2k,3(R) 
** \param[out] correc0 value to substract from S_2k,3 in order to avoid 
**                     numerical problems when R is very small
**
**  \remark Compute the random elements 
**  \remark used to simulate a GRF with a power semi-variogram GC h^(2k)log(h)
**  \remark according to the method proposed in 
**  \remark Emery, X. and Lantuejoul, C. (2008) 
**  \remark A spectral approach to simulating intrinsec random fields with power
**  \remark and spline generalized covariance. 
**  \remark In Computational Geosciences 12:121-132
**
*****************************************************************************/
static void st_set_spline_1D(int     ib,
                             double  scale,
                             int     k,
                             double *omega,
                             double *phi,
                             double *xi_3,
                             double *correc0)
{  
  double R;
  int twokm1;
  static double twoPI, twokp1s2, log3s2, logkp3s2, logkp1, coeff;
  static int twok = 0;
  static int k_mem = -1;

  if (ib == 0 || k != k_mem)
  {
    twok      = 2  * k;
    twokm1    = twok - 1;
    twokp1s2  = twok + 1./2;
    twoPI     = 2. * GV_PI;
    log3s2    = loggamma(1.5);
    logkp3s2  = loggamma(k + 1.5);
    logkp1    = loggamma(k + 1.);
    coeff     = sqrt(2 * exp(logkp3s2 + logkp1 - log3s2) / pow(GV_PI,twokm1));
    k_mem     = k;
  }

  R        = law_beta2(1./2,1./2);
  * xi_3   = coeff * sqrt((R+1.) / pow(R,twokp1s2));
  *correc0 = cos(*phi);
  *omega   = twoPI * R * pow(scale,twok);
  *phi     = twoPI * law_uniform(0.,1.);

}

/*****************************************************************************/
/*!
**  Generate a migration process - Perform the core allocation
**
** \param[in]  tmin  minimum value
** \param[in]  tmax  maximum value
** \param[in]  scale scale of the exponential
**
** \param[out] number count of values generated
**
** \remark  This procedure allocates memory that should be freed
**
*****************************************************************************/
static double *st_migration_alloc(double  tmin,
                                  double  tmax,
                                  double  scale,
                                  int    *number)
{
  double *tab,step,delta;
  int     count,n_quant, i;
 
  /* Initializations */
  
  *number = count = 0;

  delta = tmax - tmin;
  if (scale < delta * EPS)
  {
    step  = delta * EPS;
    count = (int) ceil(delta / EPS); 
    tab   = (double *) mem_alloc(sizeof(double) * count,0);
    if (tab == (double *) NULL) return(tab);
    for (i=0; i<count; i++) tab[i] = tmin + i * step;
  }
  else
  {
    n_quant = count = 1;
    tab = (double *) mem_alloc(sizeof(double) * n_quant * QUANT,0);
    if (tab == (double *) NULL) return(tab);
    tab[0] = tmin + scale * log(law_uniform(0.,1.));
    tab[1] = tmin - scale * log(law_uniform(0.,1.));
      
    while (tab[count] <= tmax)
    {
      count++;
      if (count >= n_quant * QUANT)
      {
        n_quant++;
        tab = (double *)
          mem_realloc((char *) tab, sizeof(double)*QUANT*n_quant,0);
        if (tab == (double *) NULL) return(tab);
      }
      tab[count] = tab[count-1] - scale * log(law_uniform(0.,1.));
    }
      
    /* Resize the array */
      
    count++;
    tab = (double *) mem_realloc((char *) tab,sizeof(double) * count,0);
  } 

  *number = count;
  return(tab);
}

/*****************************************************************************/
/*!
**  Generate a migration process
**
** \param[in]  tmin  minimum value
** \param[in]  tmax  maximum value
** \param[in]  scale scale of the exponential
**
** \param[out] tab    working array
** \param[out] number count of values generated
**
*****************************************************************************/
static void st_migration(double  tmin,
                         double  tmax,
                         double  scale,
                         double *tab,
                         int    *number)
{
  int count,i;
  double delta;
  
  /* Initializations */

  delta = tmax - tmin;
  if (scale < delta * EPS)
  {
    count = (int) ceil(delta / EPS);    
    for (i=0; i<count; i++) tab[i] = law_gaussian();
  }
  else
  {
    count  = 1;
    tab[0] = tmin + scale * log(law_uniform(0.,1.));
    tab[1] = tmin - scale * log(law_uniform(0.,1.));
    while (tab[count] <= tmax)
    {
      count++;
      tab[count] = tab[count-1] - scale * log(law_uniform(0.,1.));
    }
  }
  
  *number = count + 1;
  return;
}

/*****************************************************************************/
/*!
**  Generate a dilution process
**
** \param[in]  tmin minimum value
** \param[in]  tmax maximum value
** \param[in]  mesh mesh of the random walk
**
** \param[out] tab    working array
** \param[out] start  initial time
** \param[out] number count of values generated
**
*****************************************************************************/
static void st_dilution(double  tmin,
                        double  tmax,
                        double  mesh,
                        double *tab,
                        double *start,
                        int    *number)

{
  double tdeb;
  int    count;

  /* Initializations */

  count = 0;
  tdeb  = tmin - mesh * law_uniform(0.,1.);

  while (tdeb + count * mesh <= tmax)
  {
    tab[count] = (law_uniform(0.,1.) < 0.5) ? -1. : 1.;
    count++;
  }
  *start  = tdeb;
  *number = count;
  return;
}

/*****************************************************************************/
/*!
**  Generate a migration process. Perform the core allocation
**
** \param[in]  tmin minimum value
** \param[in]  tmax maximum value
** \param[in]  mesh mesh of the random walk
**
** \param[out] start  initial time
** \param[out] number count of values generated
**
** \remark  This procedure allocates memory that should be freed
**
*****************************************************************************/
static double *st_dilution_alloc(double  tmin,
                                 double  tmax,
                                 double  mesh,
                                 double *start,
                                 int    *number)
{
  double *tab,tdeb;
  int     count,n_quant;

  /* Initializations */

  n_quant = 1;
  count   = 0;
  tab = (double *) mem_alloc(sizeof(double) * n_quant * QUANT,0);
  if (tab == (double *) NULL) goto label_end;

  tdeb = tmin - mesh * law_uniform(0.,1.);

  while (tdeb + count * mesh <= tmax)
  {
    if (count >= n_quant * QUANT)
    {
      n_quant++;
      tab = (double *)
        mem_realloc((char *) tab,sizeof(double) * QUANT * n_quant,0);
      if (tab == (double *) NULL) goto label_end;
    }
    tab[count] = (law_uniform(0.,1.) < 0.5) ? -1. : 1.;
    count++;
  }

  /* Resize the array */

  tab = (double *) mem_realloc((char *) tab,sizeof(double) * count,0);
  *start  = tdeb;

label_end:
  *number = count;
  return(tab);
}

/*****************************************************************************/
/*!
**  Returns the rank of the point t0 in the Poisson point process
**
** \param[in]  def_rank Rank of the Poisson point
** \param[in]  t0 starting time
** \param[in]  t  Poisson point process
** \param[in]  nt number of Poisson point process
**
*****************************************************************************/
static int st_rank_in_poisson(int def_rank,
                              double  t0,
                              double *t,
                              int     nt)

{
  int it,itp,itn;

  /* First, try with the default interval then the next one and finally
     the previous one */

  if ( t0 >= t[def_rank] && t0 < t[def_rank+1] ) 
    return(def_rank);
  else if ( def_rank < (nt-2) &&  t0 >= t[def_rank+1] && t0 < t[def_rank+2] )
    return(def_rank+1);
  else if ( def_rank > 0 &&  t0 >= t[def_rank-1] && t0 < t[def_rank] )
    return(def_rank-1);

  /* The default value is not good ==> dichotomy */

  itp = 0;
  itn = nt - 1;
  while (itn - itp > 1)
  {
    it = (itn + itp) / 2;
    if (t0 >= t[it])
      itp = it;
    else
      itn = it;
  }
  return(itp);
}

/*****************************************************************************/
/*!
**  Returns the rank of the point t0 in a regular pavement
**
** \param[in]  t0    starting time
** \param[in]  tdeb  origin on the line
** \param[in]  scale scaling factor
** \param[in]  nt    number of Poisson point process
**
*****************************************************************************/
static int st_rank_regular(double  t0,
                           double  tdeb,
                           double  scale,
                           int     nt)

{
  int nt0;

  nt0 = (int) ((t0 - tdeb) / scale);
  if (nt0 < 0 || nt0 >= nt) messageAbort("Error in st_rank_regular");

  return(nt0);
}

/****************************************************************************/
/*!
**  Particular case of the stable model. It must be turned into:
**  - Exponential : when param is too close to 1
**  - Gaussian    : when param is too close to 2
**
**  Particular case of the K-Bessel model. It must be turned into:
**  - Exponential : when param is too close to 0.5
**
** \return  The modified type
**
*****************************************************************************/
static int st_particular_case(int type,
                              double param)
{
  static double eps = 1.e-7;
  
  switch (type)
  {
    case COV_STABLE:
      if (ABS(param - 1.) < eps) return(COV_EXPONENTIAL);
      if (ABS(param - 2.) < eps) return(COV_GAUSSIAN);
      return(COV_STABLE);
      break;
      
    case COV_BESSEL_K:
      if (ABS(param - 0.5) < eps) return(COV_EXPONENTIAL);
      break;
  }
  
  return(type);
}

/****************************************************************************/
/*!
**  Initialize the array of seeds for the generation of a simulation
**  using the Turning Bands method
**
** \return  Error return code : 1 for problem; 0 otherwise
**
** \param[in]  dbin     Input Db structure (NULL is absent)
** \param[in]  dbout    Output Db structure
** \param[in]  model    Model structure
**
** \param[out]  situba  Situba structure
**
*****************************************************************************/
static int st_initialize(Db     *dbin,
                         Db     *dbout,
                         Model  *model,
                         Situba *situba)
{
  double *t,*v0,*v1,*v2;
  double  scale,tdeb,tmin,tmax,omega,phi,correc,correc0,theta1,param;
  int     is,ib,ivar,ibs,ncova,nt,isimu,nvar,type;
  int     error,nbtuba,nbsimu,mem_seed;
    
  /* Initializations */
    
  st_density(situba);
  error  = 1;
  ncova  = model->getCovaNumber();
  nvar   = model->getVariableNumber();
  nbtuba = situba->nbtuba;
  nbsimu = situba->nbsimu;
  theta1 = 1. / situba->theta;
  t = v0 = v1 = v2 = (double *) NULL;
  nt     = 0;

  /* Loop on the turning bands */
  
  mem_seed = law_get_random_seed();
  for (ivar=0; ivar<nvar; ivar++)
    for (isimu=ibs=0; isimu<nbsimu; isimu++)
      for (is=0; is<ncova; is++)
        for (ib=0; ib<nbtuba; ib++,ibs++)
        {
          tmin    = situba->codir[ibs]->tmin;
          tmax    = situba->codir[ibs]->tmax;
          scale   = situba->codir[ibs]->scale;
          type    = model->getCovaType(is);
          param   = model->getParam(is);
          nt      = 0;
          correc0 = 0;
          type    = st_particular_case(type,param);
          SEEDS(ivar,is,ib,isimu) = law_get_random_seed();
	    
          switch (type)
          {
            case COV_NUGGET:
              // Next line is simply to let the random number cycle
              (void) law_gaussian();
              break;
	      
            case COV_EXPONENTIAL:
              scale *= 2.;
              t = st_migration_alloc(tmin,tmax,scale,&nt);
              (void) law_uniform(0.,1.);
              break;
		
            case COV_SPHERICAL:
              t = st_dilution_alloc(tmin,tmax,scale,&tdeb,&nt);
              break;
		
            case COV_CUBIC:
              t = st_dilution_alloc(tmin,tmax,scale,&tdeb,&nt);
              break;
		
            case COV_GAUSSIAN:
            case COV_SINCARD:
              st_spectral(type,scale,2.,&omega,&phi);
              break;
		
            case COV_BESSEL_J:
              st_spectral(type,scale,param,&omega,&phi);
              break;
		
            case COV_BESSEL_K:
              if (param >0.5)
                st_spectral(type,scale,param,&omega,&phi);
              else
              { 
                scale = st_compute_scale_Kb(param,scale) * 2;
                t = st_migration_alloc(tmin,tmax,scale,&nt);
                (void) law_uniform(0.,1.);
              }	
              break;
		
            case COV_STABLE:
              if (param > 1) 
                st_spectral(type,scale,param,&omega,&phi);
              else
              { 
                scale *=2;
                scale  = st_compute_scale(param,scale);
                t = st_migration_alloc(tmin,tmax,scale,&nt);
                (void) law_uniform(0.,1.);
              }
              break;
		
            case COV_POWER:
              st_set_power_1D(ib,scale,param,&omega,&phi,&correc,&correc0);
              break;
		
            case COV_SPLINE_GC:
              st_set_spline_1D(ib,scale,1,&omega,&phi,&correc,&correc0); 
              break;
	      
            case COV_LINEAR:
            case COV_ORDER1_GC:
            case COV_ORDER3_GC:
            case COV_ORDER5_GC:
              t = st_migration_alloc(tmin,tmax,theta1,&nt);
              st_irf_process_alloc(nt,type,&v0,&v1,&v2);
              break;
	      
            default:
              messerr("The structure (%s) cannot be simulated",model->getCovName(type).c_str());
              messerr("using the Turning Bands algorithm");
              goto label_end;
          }
          if (nt > situba->max_alloc) situba->max_alloc = nt;
          t  = (double *) mem_free((char *) t);
          v0 = (double *) mem_free((char *) v0);
          v1 = (double *) mem_free((char *) v1);
          v2 = (double *) mem_free((char *) v2);
        }
  law_set_random_seed(mem_seed);
  error = 0;
  
label_end:
  t  = (double *) mem_free((char *) t);
  v0 = (double *) mem_free((char *) v0);
  v1 = (double *) mem_free((char *) v1);
  v2 = (double *) mem_free((char *) v2);
  return(error);
}

/****************************************************************************/
/*!
**  Deallocates the Situba structure
**
** \return  Pointer to the freed Situba structure
**
** \param[in]  model   Model structure
** \param[in]  situba  Situba structure to be deallocated
**
*****************************************************************************/
static Situba *st_dealloc(Model  *model,
                          Situba *situba)

{
  int ibs;

  if (situba == (Situba *) NULL) return(situba);

  /* Deallocate the structures for the seeds */

  situba->seeds = (int *) mem_free((char *) situba->seeds);

  /* Deallocate the structures for the directions */

  for (ibs=0; ibs<situba->nbands; ibs++)
    situba->codir[ibs] = (Direction *) mem_free((char *) situba->codir[ibs]);
  situba->codir = (Direction **) mem_free((char *) situba->codir);

  /* Deallocate the Situba structure */

  situba = (Situba *) mem_free((char *) situba);

  return(situba);
}

/*****************************************************************************/
/*!
**  Allocates the arrays for the simulations using Turning Bands
**
** \return  Pointer to the Situba structure or NULL
**
** \param[in]  model   Model structure
** \param[in]  nbsimu  number of simulations
** \param[in]  nbtuba  number of turning bands
**
*****************************************************************************/
static Situba *st_alloc(Model *model,
                        int    nbsimu,
                        int    nbtuba)
{
  Situba *situba;
  int     i,ibs,nvar,ncova,error,size;

  error = 1;
  nvar  = model->getVariableNumber();
  ncova = model->getCovaNumber();

  /* Allocate the Situba structure */

  situba = (Situba *) mem_alloc(sizeof(Situba),0);
  if (situba == (Situba *) NULL) goto label_end;
  situba->nbands    = nbsimu * nbtuba * ncova;
  situba->nbsimu    = nbsimu;
  situba->nbtuba    = nbtuba;
  situba->max_alloc = 0;
  situba->nb_points_simu = 0;
  situba->theta     = 0.;
  situba->field     = 0.;
  situba->seeds     = (int *) NULL;
  situba->codir     = (Direction **) NULL;

  /* Allocate the structures for the seeds */

  size = nvar * ncova * nbtuba * nbsimu;
  situba->seeds = (int *) mem_alloc(sizeof(int) * size,0);
  if (situba->seeds == (int *) NULL) goto label_end;
  for (i=0; i<size; i++) situba->seeds[i] = 0;

  /* Allocate the structures for the directions */

  size = situba->nbands;
  situba->codir = (Direction **) mem_alloc(sizeof(Direction *) * size,0);
  if (situba->codir == (Direction **) NULL) goto label_end;
  for (ibs=0; ibs<size; ibs++) situba->codir[ibs] = (Direction *) NULL;

  for (ibs=0; ibs<size; ibs++)
  {
    situba->codir[ibs] = (Direction *) mem_alloc(sizeof(Direction),0);
    if (situba->codir[ibs] == (Direction *) NULL) goto label_end;
  }

  /* Set the error return flag */

  error = 0;

label_end:
  if (error) situba = st_dealloc(model,situba);
  return(situba);
}

/*****************************************************************************/
/*!
**  Perform the rotation of a set of normalized direction
**  coefficients
**
** \param[in]  situba Situba structure
** \param[in]  a      Rotation direction
** \param[in]  theta  Rotation angle
**
*****************************************************************************/
static void st_rotate_directions(Situba *situba,
                                 double  a[3],
                                 double  theta)
{
  double ct,st,dir[3];
  int    i,ibs;

  /* Initializations */

  ct = cos(theta);
  st = sin(theta);

  /* Loop on the direction coefficients */

  for (ibs=0; ibs<situba->nbands; ibs++)
  {
    for (i=0; i<3; i++) dir[i] = situba->codir[ibs]->ang[i];
    ut_rotation_direction(ct,st,a,dir);
    for (i=0; i<3; i++) situba->codir[ibs]->ang[i] = dir[i];
  }

  return;
}

/*****************************************************************************/
/*!
**  Calculates the projection of a point on a turning band
**
** \return  Projection value
**
** \param[in]  db      Db structure
** \param[in]  situba  Situba structure
** \param[in]  ibs     rank of the turning band
** \param[in]  iech    rank of the sample
**
*****************************************************************************/
static double st_project_point(Db     *db,
                               Situba *situba,
                               int     ibs,
                               int     iech)
{
  double t;
  int    idim;

  t = 0.;
  for (idim=0; idim<db->getNDim(); idim++)
    t += db->getCoordinate(iech,idim) * situba->codir[ibs]->ang[idim];

  return(t);
}

/*****************************************************************************/
/*!
**  Calculates the projection of a grid node on a turning band
**
** \return  Projection value
**
** \param[in]  db      Db structure
** \param[in]  situba  Situba structure
** \param[in]  ibs     rank of the turning band
** \param[in]  ix      grid index along X
** \param[in]  iy      grid index along Y
** \param[in]  iz      grid index along Z
**
*****************************************************************************/
static double st_project_grid(Db     *db,
                              Situba *situba,
                              int     ibs,
                              int     ix,
                              int     iy,
                              int     iz)
{
  double t,xyz[3];
  int    idim,indg[3];

  indg[0] = ix;
  indg[1] = iy;
  indg[2] = iz;
  grid_to_point(db,indg,(double *) NULL,xyz);

  t = 0.;
  for (idim=0; idim<db->getNDim(); idim++)
    t += xyz[idim] * situba->codir[ibs]->ang[idim];
  return(t);
}

/****************************************************************************/
/*!
**  Generate directions according to Van Der Corput algorithm.
**  The count of directions returned is the product of nbtuba by the
**  number of basic structures
**
** \param[in]  dbout   Output Db structure
** \param[in]  model   Model     structure
** \param[in]  situba  Situba    structure
**
*****************************************************************************/
static void st_gendir(Db     *dbout,
                      Model  *model,
                      Situba *situba)
{
  CovAniso* cova;
  double    axyz[3],x[2],d,r,theta,sqr,scale,val,t00,rot,range;
  int       i,j,n,ib,is,ibs,id,p,ncova,nbtuba,ndim,nbsimu,isimu;

  /* Initializations */

  ndim   = model->getDimensionNumber();
  ncova  = model->getCovaNumber();
  nbtuba = situba->nbtuba;
  nbsimu = situba->nbsimu;

  /* Loop on the directions */

  for (ibs=0; ibs<situba->nbands; ibs++)
  {

    /* Decomposition according to basis 2 then 3 */
    
    for (id=0; id<2; id++)
    {
      n     = 1 + ibs;
      x[id] = 0;
      d = p = id + 2;
      while (n > 0)
      {
        x[id] += (n % p) / d;
        d *= p;
        n /= p;
      }
    }
    
    /* Direction coefficients */
    
    sqr = sqrt (1. - x[1] * x[1]);
    situba->codir[ibs]->ang[0] = cos(2.* GV_PI * x[0]) * sqr;
    situba->codir[ibs]->ang[1] = sin(2.* GV_PI * x[0]) * sqr;
    situba->codir[ibs]->ang[2] = x[1];
    situba->codir[ibs]->tmin   =  1.e30;
    situba->codir[ibs]->tmax   = -1.e30;
    situba->codir[ibs]->scale  =  1.;
  }

  /* Random rotation of the directions */

  r = 0.;
  for (i=0; i<3; i++)
  {
    axyz[i] = law_gaussian();
    r += axyz[i] * axyz[i];
  }
  r = sqrt(r);
  for (i=0; i<3; i++) axyz[i] /= r;
  theta = 2. * GV_PI * law_uniform(0.,1.);
  st_rotate_directions(situba,axyz,theta);

  /* Take the anisotropy into account */

  for (isimu=ibs=0; isimu<nbsimu; isimu++)
    for (is=0; is<ncova; is++)
      for (ib=0; ib<nbtuba; ib++,ibs++)
      {
        cova = model->getCova(is);
        // If the covariance has no Range (i.e. Nugget Effect), the rest is non-sense.
        // Nevertheless this code is maintained to ensure in order not to disorganize
        // the possible drawing of random numbers.
        if (! cova->hasRange()) continue;
        if (cova->getFlagAniso())
        {
          VectorDouble ranges = cova->getScales();
          scale = 0.;
          for (i=0; i<3; i++)
          {
            val = 0.;
            if (cova->getFlagRotation())
              for (j=0; j<3; j++)
              {
                rot = (i == j) ? 1. : 0.;
                if (i < ndim && j < ndim)
                  rot = cova->getAnisoRotMat(i,j);
                range = 0.;
                if (j < ndim)
                  range = ranges[j];
                if (range > 0.)
                  val += (situba->codir[ibs]->ang[j] * rot / range);
              }
            else
            {
              range = 0.;
              if (i < ndim)
                range = ranges[i];
              if (range > 0.)
                val += (situba->codir[ibs]->ang[i] / range);
            }
            scale  += val * val;
            axyz[i] = val;
          }
          situba->codir[ibs]->scale = 1./ sqrt(scale);
          for (i=0; i<3; i++)
            situba->codir[ibs]->ang[i] = axyz[i] * situba->codir[ibs]->scale;
        }
        else
        {
          situba->codir[ibs]->scale = cova->getTheoretical();
        }

        if (is_grid(dbout))
        {
          situba->codir[ibs]->t00 = t00 = 
            st_project_grid(dbout,situba,ibs,0,0,0);
          situba->codir[ibs]->dxp = 
            st_project_grid(dbout,situba,ibs,1,0,0) - t00;
          situba->codir[ibs]->dyp = 
            st_project_grid(dbout,situba,ibs,0,1,0) - t00;
          situba->codir[ibs]->dzp = 
            st_project_grid(dbout,situba,ibs,0,0,1) - t00;
          if (cova->getType() == COV_SPHERICAL ||
              cova->getType() == COV_CUBIC )
          {
            situba->codir[ibs]->t00 /= situba->codir[ibs]->scale;
            situba->codir[ibs]->dxp /= situba->codir[ibs]->scale;
            situba->codir[ibs]->dyp /= situba->codir[ibs]->scale;
            situba->codir[ibs]->dzp /= situba->codir[ibs]->scale;
          }
        }
      }
  return;
}

/****************************************************************************/
/*!
**  Calculates the data extension for a set of turning bands
**
** \param[in]  db      Db structure
** \param[in]  model   Model structure
** \param[in]  situba  Situba structure
**
*****************************************************************************/
static void st_minmax(Db     *db,
                      Model  *model,
                      Situba *situba)
{
  double tt,delta;
  int    ibs,iech,nx,ny,nz,ix,iy,iz;

  /* Initializations */

  if (db == (Db *) NULL) return;

  if (is_grid(db))
  {
    nx = (db->getNDim() >= 1) ? db->getNX(0) : 1;
    ny = (db->getNDim() >= 2) ? db->getNX(1) : 1;
    nz = (db->getNDim() >= 3) ? db->getNX(2) : 1;

    /* Case when the data obeys to a grid organization */
    /* This test is programmed for 3-D (maximum) grid  */
    /* as the Turning Bands method is limited to 3-D   */

    for (ibs=0; ibs<situba->nbands; ibs++)
    {
      for (iz=0; iz<2; iz++)
        for (iy=0; iy<2; iy++)
          for (ix=0; ix<2; ix++)
          {
            tt = st_project_grid(db,situba,ibs,ix*(nx-1),iy*(ny-1),iz*(nz-1));
            if (tt < situba->codir[ibs]->tmin) situba->codir[ibs]->tmin = tt;
            if (tt > situba->codir[ibs]->tmax) situba->codir[ibs]->tmax = tt;
            delta = situba->codir[ibs]->tmax - situba->codir[ibs]->tmin;
            if (situba->field < delta) situba->field = delta;
          }
    }
  }
  else
  {

    /* Case of an isolated set of data */

    for (iech=0; iech<db->getSampleNumber(); iech++)
    {
      for (ibs=0; ibs<situba->nbands; ibs++)
      {
        if (! db->isActive(iech)) continue;
        tt = st_project_point(db,situba,ibs,iech);
        if (tt < situba->codir[ibs]->tmin) situba->codir[ibs]->tmin = tt;
        if (tt > situba->codir[ibs]->tmax) situba->codir[ibs]->tmax = tt;
        delta = situba->codir[ibs]->tmax - situba->codir[ibs]->tmin;
        if (situba->field < delta) situba->field = delta;
      }
    }
  }
  situba->nb_points_simu += db->getSampleNumber();

  return;
}

/*****************************************************************************/
/*!
**  Convert the non conditional simulations at the data points
**  into simulation error
**
** \param[in]  dbin       Input Db structure
** \param[in]  situba     Situba structure
** \param[in]  nvar       Number of variables
** \param[in]  icase      Case for PGS or GRF
** \param[in]  flag_pgs   1 if called from PGS
** \param[in]  flag_gibbs 1 if called from Gibbs
**
*****************************************************************************/
static void st_difference(Db     *dbin,
                          Situba *situba,
                          int     nvar,
                          int     icase,
                          int     flag_pgs,
                          int     flag_gibbs)
{
  double zvar,simval,simunc;
  int    iech,ivar,isimu,nbsimu;
  char string[100];

  /* Optional general title */

  nbsimu = situba->nbsimu;
  if (debug_query("simulate"))
  {
    mestitle(1,"Difference between Data and NC Simulation");
    tab_prints(NULL,1,GD_J_RIGHT,"Sample");
  }

  /* Transform the non conditional simulation into simulation error */

  if (! flag_pgs)
  {
    /********************************/
    /* Standard case (multivariate) */
    /********************************/
    
    /* Optional Header */

    if (debug_query("simulate"))
    {
      for (ivar=0; ivar<nvar; ivar++)
      {
        (void) sprintf(string,"Data%d",ivar+1);
        tab_prints(NULL,1,GD_J_RIGHT,string);
        
        for (isimu=0; isimu<nbsimu; isimu++)
        {
          (void) sprintf(string,"Simu%d",isimu+1);
          tab_prints(NULL,1,GD_J_RIGHT,string);
        }
      }
      message("\n");
    }
    
    /* Processing */
    
    for (iech=0; iech<dbin->getSampleNumber(); iech++)
    {
      if (! dbin->isActive(iech)) continue;
      if (debug_query("simulate")) tab_printi(NULL,1,GD_J_RIGHT,iech+1);
      for (ivar=0; ivar<nvar; ivar++)
      {
        zvar = TEST;
        if (! flag_gibbs)
        {
          zvar = dbin->getVariable(iech,ivar);
          if (debug_query("simulate")) tab_printg(NULL,1,GD_J_RIGHT,zvar);
        }
        for (isimu=0; isimu<nbsimu; isimu++)
        {
          if (flag_gibbs)
          {
            zvar = dbin->getSimvar(LOC_GAUSFAC,iech,isimu,ivar,0,nbsimu,nvar);
            if (debug_query("simulate")) tab_printg(NULL,1,GD_J_RIGHT,zvar);
          }
          simval = dbin->getSimvar(LOC_SIMU, iech, isimu, ivar, icase,
              nbsimu, nvar);
          if (FLAG_DGM)
          {
            simval = R_COEFF * simval + sqrt(1. - R_COEFF*R_COEFF) * law_gaussian();
          }
          if (debug_query("simulate") && !FFFF(zvar))
          {
            tab_printg(NULL, 1, GD_J_RIGHT, simval);
          }
          simunc = (!FFFF(zvar) && !FFFF(simval)) ?  simval - zvar : TEST;
          dbin->setSimvar(LOC_SIMU, iech, isimu, ivar, icase, nbsimu, nvar,
                          simunc);
        }
      }
      if (debug_query("simulate")) message("\n");
    }
  }
  else
  {
    
    /*********************************************************/
    /* Case of PGS: Data varies per simulation (monovariate) */
    /*********************************************************/

    /* Optional Header */
    
    if (debug_query("simulate"))
    {
      for (isimu=0; isimu<nbsimu; isimu++)
      {
        (void) sprintf(string,"Data%d",isimu+1);
        tab_prints(NULL,1,GD_J_RIGHT,string);
        
        (void) sprintf(string,"Simulation%d",isimu+1);
        tab_prints(NULL,1,GD_J_RIGHT,string);
      }
      message("\n");
    }
    
    /* Processing */
    
    for (iech=0; iech<dbin->getSampleNumber(); iech++)
    {
      if (! dbin->isActive(iech)) continue;
      if (debug_query("simulate")) tab_printi(NULL,1,GD_J_RIGHT,iech+1);
      for (isimu=0; isimu<nbsimu; isimu++)
      {
        zvar = dbin->getSimvar(LOC_GAUSFAC,iech,isimu,0,icase,nbsimu,1);
        if (debug_query("simulate"))
        {
          tab_printg(NULL,1,GD_J_RIGHT,zvar);
          if (! FFFF(zvar)) 
          {
            simval = dbin->getSimvar(LOC_SIMU,iech,isimu,0,icase,
                                nbsimu,1);
            tab_printg(NULL,1,GD_J_RIGHT,simval);
          }
        }
        if (! FFFF(zvar)) 
          dbin->updSimvar(LOC_SIMU,iech,isimu,0,icase,nbsimu,1,0,-zvar);
      }
      if (debug_query("simulate")) message("\n");
    }
  }

  return;
}

/*****************************************************************************/
/*!
**  Elaborates the title for the mes_process statement
**
** \param[in]  string  Resulting title
** \param[in]  title   Generic title
**
*****************************************************************************/
static void st_title_process(char *string,
                             const char *title)
{
  (void) strcpy(string,title);
  if (MES_NPGS > 1)
    (void) sprintf(&string[strlen(string)],
                   " (PGS %d/%d)",MES_IPGS+1,MES_NPGS);
  (void) sprintf(&string[strlen(string)],
                 " (GRF %d/%d)",MES_IGRF+1,MES_NGRF);
  return;
}

/*****************************************************************************/
/*!
**  Perform non-conditional simulations on a grid using the
**  Turning Bands method
**
** \return  Error return code :
** \return    0 no problem
** \return    1 a structure cannot be simulated
** \return    2 no structure to be simulated 1 core problem
**
** \param[in]  db         Db structure
** \param[in]  model      Model structure
** \param[in]  situba     Situba structure
** \param[in]  aic        Array 'aic'
** \param[in]  icase      Rank of PGS or GRF
** \param[in]  shift      Shift before writing the simulation result
**
*****************************************************************************/
static int st_simulate_grid(Db     *db,
                            Model  *model,
                            Situba *situba,
                            double *aic,
                            int     icase,
                            int     shift)
{
  double  * t,*v0,*v1,*v2,*tab,correc,correc0,vexp;
  double    tmin,tmax,tdeb,omega,phi,scale,param,t0,dt0,norme;
  double    t00,t0y,t0z,dxp,dyp,dzp,dt,theta1;
  double    c1,s1,c0x,s0x,c0y,s0y,c0z,s0z,cxp,sxp,cyp,syp,czp,szp;
  int       iech,ivar,jvar,is,ind,ib,ibs,ix,iy,iz,nvar,nt,nt0,nech;
  int       type,error,nbtuba,nbsimu,nx,ny,nz,isimu,ncova,istat;
  char      string[100];
  static double vexp1 = 0.1;
  static double vexp2 = 0.1967708298;

  /* Initializations */

  nbsimu = situba->nbsimu;
  nbtuba = situba->nbtuba;
  theta1 = 1. / situba->theta;
  nvar   = model->getVariableNumber();
  ncova  = model->getCovaNumber();
  nx     = (db->getNDim() >= 1) ? db->getNX(0) : 1;
  ny     = (db->getNDim() >= 2) ? db->getNX(1) : 1;
  nz     = (db->getNDim() >= 3) ? db->getNX(2) : 1;
  nech   = nx * ny * nz;
  norme  = sqrt(1. / nbtuba);
  nt     = iech = 0;
  vexp   = 0.;
  t = v0 = v1 = v2 = tab = (double *) NULL;
 
  /* Core allocation */

  error = 1;
  if (situba->max_alloc > 0)
  {
    t    = (double *) mem_alloc(sizeof(double) * situba->max_alloc,0);
    if (t  == (double *) NULL) goto label_end;
    v0   = (double *) mem_alloc(sizeof(double) * situba->max_alloc,0);
    if (v0 == (double *) NULL) goto label_end;
    v1   = (double *) mem_alloc(sizeof(double) * situba->max_alloc,0);
    if (v1 == (double *) NULL) goto label_end;
    v2   = (double *) mem_alloc(sizeof(double) * situba->max_alloc,0);
    if (v2 == (double *) NULL) goto label_end;
  }
  tab  = (double *) mem_alloc(sizeof(double) * nech,0);
  if (tab == (double *) NULL) goto label_end;
  
  /*****************************/
  /* Performing the simulation */
  /*****************************/
  
  st_title_process(string,"Non-conditional Simulation on Grid");
  for (ivar=istat=0; ivar<nvar; ivar++)
    for (isimu=ibs=0; isimu<nbsimu; isimu++)
      for (is=0; is<ncova; is++)
        for (ib=0; ib<nbtuba; ib++,ibs++,istat++)
        {
          mes_process(string,nbsimu*nvar*ncova*nbtuba,istat);
          tmin       = situba->codir[ibs]->tmin;
          tmax       = situba->codir[ibs]->tmax;
          dxp        = situba->codir[ibs]->dxp;
          dyp        = situba->codir[ibs]->dyp;
          dzp        = situba->codir[ibs]->dzp;
          t00        = situba->codir[ibs]->t00;
          scale      = situba->codir[ibs]->scale;
          type       = model->getCovaType(is);
          param      = model->getParam(is);
          correc     = 1.;
          correc0    = 0.;
          type       = st_particular_case(type,param);
          law_set_random_seed(SEEDS(ivar,is,ib,isimu));
          
          switch (type)
          {
            case COV_NUGGET:
              break;

            case COV_STABLE:
              if (param > 1)
              {
                correc = sqrt(2.);
                st_spectral(type,scale,param,&omega,&phi);
		
                cxp = cos(omega * dxp);
                sxp = sin(omega * dxp);
                cyp = cos(omega * dyp);
                syp = sin(omega * dyp);
                czp = cos(omega * dzp);
                szp = sin(omega * dzp);
		
                c0z = cos(omega * t00 + phi);
                s0z = sin(omega * t00 + phi);
                for (iz=ind=0; iz<nz; iz++)
                {
                  c0y = c0z;
                  s0y = s0z;
                  c1  = c0z * czp - s0z * szp;
                  s1  = s0z * czp + c0z * szp;
                  c0z = c1;
                  s0z = s1;
                  for (iy=0; iy<ny; iy++)
                  {
                    c0x = c0y;
                    s0x = s0y;
                    c1  = c0y * cyp - s0y * syp;
                    s1  = s0y * cyp + c0y * syp;
                    c0y = c1;
                    s0y = s1;
                    for (ix=0; ix<nx; ix++,ind++)
                    {
                      if (db->isActive(ind)) tab[ind] = c0x - correc0;
                      c1  = c0x * cxp - s0x * sxp;
                      s1  = s0x * cxp + c0x * sxp;
                      c0x = c1;
                      s0x = s1;
                    }
                  }
                }
              }
              else
              { 
                scale*=2;
                scale = st_compute_scale(param,scale);
                st_migration(tmin,tmax,scale,t,&nt);
                vexp = 1. - vexp1 + vexp2 * law_uniform(0.,1.);
                t0z = t00;
		
                nt0 = 0;
                for (iz=ind=0; iz<nz; iz++)
                {
                  t0y  = t0z;
                  t0z += dzp;
                  for (iy=0; iy<ny; iy++)
                  {
                    t0   = t0y;
                    t0y += dyp;
                    for (ix=0; ix<nx; ix++,ind++)
                    {
                      if (db->isActive(ind))
                      {
                        nt0 = st_rank_in_poisson(nt0,t0,t,nt);
                        tab[ind] = (2. * t0 > t[nt0+1]+t[nt0]) ? -vexp : vexp;
                      }
                      t0 += dxp;
                    }
                  }
                }
              }
              break;
	      
            case COV_BESSEL_K:
              if (param > 0.5)
              {
                correc = sqrt(2.);
                st_spectral(type,scale,param,&omega,&phi);
		
                cxp = cos(omega * dxp);
                sxp = sin(omega * dxp);
                cyp = cos(omega * dyp);
                syp = sin(omega * dyp);
                czp = cos(omega * dzp);
                szp = sin(omega * dzp);
		
                c0z = cos(omega * t00 + phi);
                s0z = sin(omega * t00 + phi);
                for (iz=ind=0; iz<nz; iz++)
                {
                  c0y = c0z;
                  s0y = s0z;
                  c1  = c0z * czp - s0z * szp;
                  s1  = s0z * czp + c0z * szp;
                  c0z = c1;
                  s0z = s1;
                  for (iy=0; iy<ny; iy++)
                  {
                    c0x = c0y;
                    s0x = s0y;
                    c1  = c0y * cyp - s0y * syp;
                    s1  = s0y * cyp + c0y * syp;
                    c0y = c1;
                    s0y = s1;
                    for (ix=0; ix<nx; ix++,ind++)
                    {
                      if (db->isActive(ind)) tab[ind] = c0x - correc0;
                      c1  = c0x * cxp - s0x * sxp;
                      s1  = s0x * cxp + c0x * sxp;
                      c0x = c1;
                      s0x = s1;
                    }
                  }
                }
              }
              else
              { 
                scale = st_compute_scale_Kb(param,scale) * 2;
                st_migration(tmin,tmax,scale,t,&nt);
                vexp = 1. - vexp1 + vexp2 * law_uniform(0.,1.);
                t0z = t00;
		
                nt0 = 0;
                for (iz=ind=0; iz<nz; iz++)
                {
                  t0y  = t0z;
                  t0z += dzp;
                  for (iy=0; iy<ny; iy++)
                  {
                    t0   = t0y;
                    t0y += dyp;
                    for (ix=0; ix<nx; ix++,ind++)
                    {
                      if (db->isActive(ind))
                      {
                        nt0 = st_rank_in_poisson(nt0,t0,t,nt);
                        tab[ind] = (2. * t0 > t[nt0+1]+t[nt0]) ? -vexp : vexp;
                      }
                      t0 += dxp;
                    }
                  }
                }
              }
              break;

            case COV_EXPONENTIAL:
              scale *= 2.;
              st_migration(tmin,tmax,scale,t,&nt);
              vexp = 1. - vexp1 + vexp2 * law_uniform(0.,1.);
	      
              t0z = t00;
              nt0 = 0;
              for (iz=ind=0; iz<nz; iz++)
              {
                t0y  = t0z;
                t0z += dzp;
                for (iy=0; iy<ny; iy++)
                {
                  t0   = t0y;
                  t0y += dyp;
                  for (ix=0; ix<nx; ix++,ind++)
                  {
                    if (db->isActive(ind))
                    {
                      nt0 = st_rank_in_poisson(nt0,t0,t,nt);
                      tab[ind] = (2. * t0 > t[nt0+1]+t[nt0]) ? -vexp : vexp;
                    }
                    t0 += dxp;
                  }
                }
              }
              break;
	      
            case COV_SPHERICAL:
              correc = sqrt(3.);
              st_dilution(tmin,tmax,scale,t,&tdeb,&nt);
              tdeb /= scale;
              t0z   = t00;
              for (iz=ind=0; iz<nz; iz++)
              {
                t0y  = t0z;
                t0z += dzp;
                for (iy=0; iy<ny; iy++)
                {
                  t0   = t0y;
                  t0y += dyp;
                  for (ix=0; ix<nx; ix++,ind++)
                  {
                    if (db->isActive(ind))
                    {
                      dt  = (t0 - tdeb);
                      nt0 = (int) dt;
                      dt0 = dt -  nt0;
                      tab[ind] = t[nt0] * (2. * dt0  - 1.);
                    }
                    t0 += dxp;
                  }
                }
              }
              break;

            case COV_CUBIC:
              correc = sqrt(840.);
              st_dilution(tmin,tmax,scale,t,&tdeb,&nt);
              tdeb /= scale;
              t0z   = t00;
              for (iz=ind=0; iz<nz; iz++)
              {
                t0y  = t0z;
                t0z += dzp;
                for (iy=0; iy<ny; iy++)
                {
                  t0   = t0y;
                  t0y += dyp;
                  for (ix=0; ix<nx; ix++,ind++)
                  {
                    if (db->isActive(ind))
                    {
                      dt  = (t0 - tdeb);
                      nt0 = (int) dt;
                      dt0 = dt -  nt0;
                      tab[ind] = t[nt0] * dt0 * (dt0 - 0.5) * (dt0 - 1.);
                    }
                    t0 += dxp;
                  }
                }
              }
              break;
	      
            case COV_GAUSSIAN:
            case COV_SINCARD:
            case COV_POWER:
            case COV_SPLINE_GC:
            case COV_BESSEL_J:
              correc = sqrt(2.);
              switch (type)
              {
                case COV_POWER:
                  st_set_power_1D(ib,scale,param,&omega,&phi,&correc,&correc0);
                  break;
		
                case  COV_SPLINE_GC:
                  st_set_spline_1D(ib,scale,1,&omega,&phi,&correc,&correc0); 
                  break;
		
                case COV_GAUSSIAN:
                case COV_SINCARD:
                  st_spectral(type,scale,2.,&omega,&phi);
                  break;
		
                case COV_BESSEL_J:
                  st_spectral(type,scale,param,&omega,&phi);
                  break;
              }

              cxp = cos(omega * dxp);
              sxp = sin(omega * dxp);
              cyp = cos(omega * dyp);
              syp = sin(omega * dyp);
              czp = cos(omega * dzp);
              szp = sin(omega * dzp);
	      
              c0z = cos(omega * t00 + phi);
              s0z = sin(omega * t00 + phi);
              for (iz=ind=0; iz<nz; iz++)
              {
                c0y = c0z;
                s0y = s0z;
                c1  = c0z * czp - s0z * szp;
                s1  = s0z * czp + c0z * szp;
                c0z = c1;
                s0z = s1;
                for (iy=0; iy<ny; iy++)
                {
                  c0x = c0y;
                  s0x = s0y;
                  c1  = c0y * cyp - s0y * syp;
                  s1  = s0y * cyp + c0y * syp;
                  c0y = c1;
                  s0y = s1;
                  for (ix=0; ix<nx; ix++,ind++)
                  {
                    if (db->isActive(ind)) tab[ind] = c0x - correc0;
                    c1  = c0x * cxp - s0x * sxp;
                    s1  = s0x * cxp + c0x * sxp;
                    c0x = c1;
                    s0x = s1;
                  }
                }
              }
              break;
	      
            case COV_LINEAR:
            case COV_ORDER1_GC:
            case COV_ORDER3_GC:
            case COV_ORDER5_GC:
              st_migration(tmin,tmax,theta1,t,&nt);
              st_irf_process(nt,type,t,v0,v1,v2);
              correc = st_irf_correc(type,theta1,scale);

              t0z = t00;
              nt0 = 0;
              for (iz=ind=0; iz<nz; iz++)
              {
                t0y  = t0z;
                t0z += dzp;
                for (iy=0; iy<ny; iy++)
                {
                  t0   = t0y;
                  t0y += dyp;
                  for (ix=0; ix<nx; ix++,ind++)
                  {
                    if (db->isActive(ind))
                    {
                      nt0 = st_rank_in_poisson(nt0,t0,t,nt);
                      tab[ind] = st_irf_process_sample(type,nt0,t0,t,v0,v1,v2);
                    }
                    t0 += dxp;
                  }
                }
              }
              break;
              
            default:
              break;
          }

          if (type != COV_NUGGET)
            for (iech=0; iech<nech; iech++)
              if (db->isActive(iech))
                for (jvar=0; jvar<nvar; jvar++)
                  db->updSimvar(LOC_SIMU, iech, shift + isimu, jvar, icase,
                                nbsimu, nvar, 0,
                                tab[iech] * correc * AIC(is, jvar, ivar));
        }

  /* Normation */

  for (isimu=0; isimu<nbsimu; isimu++)
    for (iech=0; iech<nech; iech++)
      for (jvar=0; jvar<nvar; jvar++)
        if (db->isActive(iech))
          db->updSimvar(LOC_SIMU, iech, shift + isimu, jvar, icase, nbsimu,
                        nvar, 1, norme);

  error = 0;

label_end:
  t    = (double *) mem_free((char *) t);
  v0   = (double *) mem_free((char *) v0);
  v1   = (double *) mem_free((char *) v1);
  v2   = (double *) mem_free((char *) v2);
  tab  = (double *) mem_free((char *) tab);

  return(error);
}

/*****************************************************************************/
/*!
**  Add the contribution of the nugget effect to the non-conditional 
**  simulations
**
** \param[in]  db         Db structure
** \param[in]  file_type  File type (DATA or RESULT)
** \param[in]  model      Model structure
** \param[in]  situba     Situba structure
** \param[in]  aic        Array 'aic'
** \param[in]  icase      Rank of PGS or GRF
**
*****************************************************************************/
static void st_simulate_nugget(Db     *db,
                               int     file_type,
                               Model  *model,
                               Situba *situba,
                               double *aic,
                               int     icase)
{
  double  nugget;
  int     iech,is,ivar,jvar,ncova,nech,isimu,nvar,type,nbtuba,nbsimu,flag_used;

  /* Initializations */

  nech   = db->getSampleNumber();
  ncova  = model->getCovaNumber();
  nvar   = model->getVariableNumber();
  nbtuba = situba->nbtuba;
  nbsimu = situba->nbsimu;

  /* Do nothing if there is no nugget effect in the model */

  flag_used = 0;
  for (is=0; is<ncova && flag_used == 0; is++)
  {
    if (model->getCovaType(is) == COV_NUGGET) flag_used = 1;
  }
  if (! flag_used) return;

  /* Performing the simulation */

  for (isimu=0; isimu<nbsimu; isimu++)
    for (ivar=0; ivar<nvar; ivar++)
      for (is=0; is<ncova; is++)
      {
        type = model->getCovaType(is);
        
        if (type != COV_NUGGET) continue;
        law_set_random_seed(SEEDS(ivar,is,0,isimu));
        
        for (iech=0; iech<nech; iech++)
        {
          if (! db->isActive(iech)) continue;
          nugget = law_gaussian();
          for (jvar=0; jvar<nvar; jvar++)
            db->updSimvar(LOC_SIMU, iech, isimu, jvar, icase, nbsimu, nvar, 0,
                          nugget * AIC(is, jvar, ivar));
        }
      }
  
  return;
}

/*****************************************************************************/
/*!
**  Perform non-conditional simulations on a set of points using
**  Turning Bands method.
**
** \return  Error return code :
** \return    0 no problem
** \return    1 a structure cannot be simulated
** \return    2 no structure to be simulated 1 core problem
**
** \param[in]  db         Db structure
** \param[in]  file_type  File type (DATA or RESULT)
** \param[in]  model      Model structure
** \param[in]  situba     Situba structure
** \param[in]  aic        Array 'aic'
** \param[in]  icase      Rank of PGS or GRF
** \param[in]  shift      Shift before writing the simulation result
**
*****************************************************************************/
static int st_simulate_point(Db     *db,
                             int     file_type,
                             Model  *model,
                             Situba *situba,
                             double *aic,
                             int     icase,
                             int     shift)
{
  double *t,*v0,*v1,*v2,*tab,correc,correc0,vexp;
  double  tmin,tmax,tdeb,omega,phi,param,scale,t0,dt0,norme,r,theta1;
  int     iech,is,ib,ivar,jvar,ibs,ncova,nt,nt0,nech,isimu,nvar,type;
  int     error,nbtuba,nbsimu,istat;
  char    string[100];
  static  double vexp1 = 0.1;
  static  double vexp2 = 0.1967708298;

  /* Initializations */

  nech   = db->getSampleNumber();
  ncova  = model->getCovaNumber();
  nvar   = model->getVariableNumber();
  nbtuba = situba->nbtuba;
  nbsimu = situba->nbsimu;
  theta1 = 1. / situba->theta;
  norme  = sqrt(1. / nbtuba);
  nt     = 0;
  vexp   = 0.;
  t = v0 = v1 = v2 = tab = (double *) NULL;

  /* Core allocation */

  error = 1;
  if (situba->max_alloc > 0)
  {
    t    = (double *) mem_alloc(sizeof(double) * situba->max_alloc,0);
    if (t  == (double *) NULL) goto label_end;
    v0   = (double *) mem_alloc(sizeof(double) * situba->max_alloc,0);
    if (v0 == (double *) NULL) goto label_end;
    v1   = (double *) mem_alloc(sizeof(double) * situba->max_alloc,0);
    if (v1 == (double *) NULL) goto label_end;
    v2   = (double *) mem_alloc(sizeof(double) * situba->max_alloc,0);
    if (v2 == (double *) NULL) goto label_end;
  }
  tab  = (double *) mem_alloc(sizeof(double) * nech,0);
  if (tab == (double *) NULL) goto label_end;

  /*****************************/
  /* Performing the simulation */
  /*****************************/

  st_title_process(string,"Non-conditional Simulation on Points");
  for (ivar=istat=0; ivar<nvar; ivar++)
    for (isimu=ibs=0; isimu<nbsimu; isimu++)
      for (is=0; is<ncova; is++)
        for (ib=0; ib<nbtuba; ib++,ibs++,istat++)
        {
          mes_process(string,nbsimu*nvar*ncova*nbtuba,istat);
          tmin    = situba->codir[ibs]->tmin;
          tmax    = situba->codir[ibs]->tmax;
          scale   = situba->codir[ibs]->scale;
          param   = model->getParam(is);
          type    = model->getCovaType(is);
          correc  = 1.;
          correc0 = 0.;
          type    = st_particular_case(type,param);
          law_set_random_seed(SEEDS(ivar,is,ib,isimu));

          switch (type)
          {
            case COV_NUGGET:
              break;
	      
            case COV_STABLE:
              if (param > 1)
              {
                correc = sqrt(2.);
                st_spectral(type,scale,param,&omega,&phi);
                for (iech=0; iech<nech; iech++)
                {
                  if (! db->isActive(iech)) continue;
                  t0  = st_project_point(db,situba,ibs,iech);
                  tab[iech] = cos(omega * t0 + phi);
                }
              }
              else
              { 
                scale*=2;
                scale = st_compute_scale(param,scale);
                st_migration(tmin,tmax,scale,t,&nt);
                vexp = 1. - vexp1 + vexp2 * law_uniform(0.,1.);
		
                for (iech=nt0=0; iech<nech; iech++)
                {
                  if (! db->isActive(iech)) continue;
                  t0  = st_project_point(db,situba,ibs,iech);
                  nt0 = st_rank_in_poisson(nt0,t0,t,nt);
                  tab[iech] = (2. * t0 > t[nt0+1]+t[nt0]) ? -vexp : vexp;
                }
              }
              break;
	      
            case COV_EXPONENTIAL:
              scale *= 2.;
              st_migration(tmin,tmax,scale,t,&nt);
              vexp = 1. - vexp1 + vexp2 * law_uniform(0.,1.);
              
              for (iech=nt0=0; iech<nech; iech++)
              {
                if (! db->isActive(iech)) continue;
                t0  = st_project_point(db,situba,ibs,iech);
                nt0 = st_rank_in_poisson(nt0,t0,t,nt);
                tab[iech] = (2. * t0 > t[nt0+1]+t[nt0]) ? -vexp : vexp;
              }
              break;

            case COV_SPHERICAL:
              correc = sqrt(3.);
              st_dilution(tmin,tmax,scale,t,&tdeb,&nt);
              for (iech=0; iech<nech; iech++)
              {
                if (! db->isActive(iech)) continue;
                t0  = st_project_point(db,situba,ibs,iech);
                nt0 = st_rank_regular(t0,tdeb,scale,nt);
                dt0 = (t0 - tdeb) - scale * nt0;
                r   = dt0 / scale;
                tab[iech] = t[nt0] * (2. * r - 1.);
              }
              break;

            case COV_CUBIC:
              correc = sqrt(840.);
              st_dilution(tmin,tmax,scale,t,&tdeb,&nt);
              for (iech=0; iech<nech; iech++)
              {
                if (! db->isActive(iech)) continue;
                t0  = st_project_point(db,situba,ibs,iech);
                nt0 = st_rank_regular(t0,tdeb,scale,nt);
                dt0 = (t0 - tdeb) - scale * nt0;
                r   = dt0 / scale;
                tab[iech] = t[nt0] * r * (r - 0.5) * (r - 1.);
              }
              break;
	      
            case COV_POWER:
              st_set_power_1D(ib,scale,param,&omega,&phi,&correc,&correc0);
              for (iech=0; iech<nech; iech++)
              {
                if (! db->isActive(iech)) continue;
                t0  = st_project_point(db,situba,ibs,iech);
                tab[iech] = cos(omega * t0 + phi)-correc0;
              }
              break;
	    
            case COV_SPLINE_GC:
              st_set_spline_1D(ib,scale,1,&omega,&phi,&correc,&correc0); 
              for (iech=0; iech<nech; iech++)
              {
                if (! db->isActive(iech)) continue;
                t0  = st_project_point(db,situba,ibs,iech);
                tab[iech] = cos(omega * t0 + phi)-correc0;
              }
              break; 
	    
            case COV_GAUSSIAN:
            case COV_SINCARD:
            case COV_BESSEL_J:
            case COV_BESSEL_K:
              correc = sqrt(2.);
              switch (type)
              {
                case COV_GAUSSIAN:
                case COV_SINCARD:
                  st_spectral(type,scale,2.,&omega,&phi);
                  break;
	      
                case COV_BESSEL_J:
                case COV_BESSEL_K:
                  st_spectral(type,scale,param,&omega,&phi);
                  break;
              }
	    
              for (iech=0; iech<nech; iech++)
              {
                if (! db->isActive(iech)) continue;
                t0  = st_project_point(db,situba,ibs,iech);
                tab[iech] = cos(omega * t0 + phi);
              }
              break;
	    
            case COV_LINEAR:
            case COV_ORDER1_GC:
            case COV_ORDER3_GC:
            case COV_ORDER5_GC:
              st_migration(tmin,tmax,theta1,t,&nt);
              st_irf_process(nt,type,t,v0,v1,v2);
              correc = st_irf_correc(type,theta1,scale);
	    
              for (iech=nt0=0; iech<nech; iech++)
              {
                if (! db->isActive(iech)) continue;
                t0  = st_project_point(db,situba,ibs,iech);
                nt0 = st_rank_in_poisson(nt0,t0,t,nt);
                tab[iech] = st_irf_process_sample(type,nt0,t0,t,v0,v1,v2);
              }
              break;
	    
            default:
              break;
          }
	  
          if (type != COV_NUGGET)
            for (iech=0; iech<nech; iech++)
              if (db->isActive(iech))
                for (jvar=0; jvar<nvar; jvar++)
                  db->updSimvar(LOC_SIMU, iech, shift + isimu, jvar, icase,
                                nbsimu, nvar, 0,
                                tab[iech] * correc * AIC(is, jvar, ivar));
        }
  
  /* Normation */
  
  for (isimu=0; isimu<nbsimu; isimu++)
    for (iech=0; iech<nech; iech++)
      for (jvar=0; jvar<nvar; jvar++)
        if (db->isActive(iech))
          db->updSimvar(LOC_SIMU, iech, shift + isimu, jvar, icase, nbsimu,
                        nvar, 1, norme);
  
  /* Set the error return code */
  
  error = 0;
  
label_end:
  t    = (double *) mem_free((char *) t);
  v0   = (double *) mem_free((char *) v0);
  v1   = (double *) mem_free((char *) v1);
  v2   = (double *) mem_free((char *) v2);
  tab  = (double *) mem_free((char *) tab);

  return(error);
}

/****************************************************************************/
/*!
**  Update the conditional simulations when the target coincides
**  with a data point
**
** \param[in]  dbin      Input Db structure
** \param[in]  dbout     Output Db structure
** \param[in]  model     Model structure
** \param[in]  situba    Situba structure
** \param[in]  icase     Case for PGS or GRF
** \param[in]  flag_pgs  1 if called from PGS
**
** \remarks This migration is not performed in the case where data point
** \remarks coincide with the target artificially. This is the case
** \remarks for the Discrete Gaussian Model (DGM) where data have been
** \remarks migrated to the cell center to mimic a point randomized 
** \remarks within a cell
**
*****************************************************************************/
static void st_update_data2target(Db     *dbin,
                                  Db     *dbout,
                                  Model  *model,
                                  Situba *situba,
                                  int     icase,
                                  int     flag_pgs)
{
  double *coor1,*coor2,dist,eps,eps2,valdat,radius,delta;
  int    *indg,ip,idim,ip_close,isimu,ivar,nvar,ndim,nbsimu,ik;

  /* Initialization */

  if (dbin->getSampleNumber() <= 0) return;
  if (FLAG_DGM) return;
  nvar   = model->getVariableNumber();
  ndim   = dbin->getNDim();
  nbsimu = situba->nbsimu;

  coor1 = db_vector_alloc(dbin);
  coor2 = db_vector_alloc(dbout);
  indg  = db_indg_alloc(dbout);

  /* Calculate the field extension */

  (void) db_extension_diag(dbin,&radius);
  eps  = radius * 1.e-6;
  eps2 = eps * eps;

  /* Dispatch according to the file type */

  if (is_grid(dbout))
  {

    /*********************************************/
    /* Case where the output file is a grid file */
    /*********************************************/

    for (ip=0; ip<dbin->getSampleNumber(); ip++)
    {
      if (! dbin->isActive(ip)) continue;
      db_sample_load(dbin,LOC_X,ip,coor2);
      if (point_to_grid(dbout,coor2,1,indg)) continue;
      ik = db_index_grid_to_sample(dbout,indg);
      if (! dbout->isActive(ik)) continue;
      grid_to_point(dbout,indg,(double *) NULL,coor1);

      /* Get the distance to the target point */
      
      dist = 0;
      for (idim=0; idim<ndim; idim++)
      {
        delta = coor1[idim] - coor2[idim];
        dist += delta * delta;
      }
      if (dist > eps2) continue;
      
      /* We have found a close data point: perform the assignment */
      
      for (isimu=0; isimu<nbsimu; isimu++)
        for (ivar=0; ivar<nvar; ivar++)
        {
          if (! flag_pgs)
            valdat = dbin->getVariable(ip,ivar);
          else
            valdat = dbin->getSimvar(LOC_GAUSFAC,ip,isimu,0,icase,nbsimu,1);
          if (FFFF(valdat)) continue;
          dbout->setSimvar(LOC_SIMU,ik,isimu,ivar,icase,nbsimu,nvar,valdat);
        }
    }
  }
  else
  {
    
    /**********************************************/
    /* Case where the output file is a point file */
    /**********************************************/
    
    for (ik=0; ik<dbout->getSampleNumber(); ik++)
    {
      if (! dbout->isActive(ik)) continue;
      db_sample_load(dbout,LOC_X,ik,coor1);

      /* Look for the closest data point */

      ip_close = -1;
      for (ip=0; ip<dbin->getSampleNumber() && ip_close<0; ip++)
      {
        if (! dbin->isActive(ip)) continue;
        db_sample_load(dbin,LOC_X,ip,coor2);

        /* Get the distance to the target point */

        dist = 0;
        for (idim=0; idim<ndim; idim++)
        {
          delta = coor1[idim] - coor2[idim];
          dist += delta * delta;
        }
        if (dist <= eps2) ip_close = ip;
      }

      if (ip_close < 0) continue;

      /* We have found a close data point: perform the assignment */

      for (isimu=0; isimu<nbsimu; isimu++)
        for (ivar=0; ivar<nvar; ivar++)
        {
          if (! flag_pgs)
            valdat = dbin->getVariable(ip_close,ivar);
          else
            valdat = dbin->getSimvar(LOC_GAUSFAC,ip_close,isimu,0,icase,
                                nbsimu,1);
          if (FFFF(valdat)) continue;
          dbout->setSimvar(LOC_SIMU,ik,isimu,ivar,icase,nbsimu,nvar,valdat);
        }
    }
  }

  coor1 = db_vector_free(coor1);
  coor2 = db_vector_free(coor2);
  indg  = db_indg_free(indg);
  return;
}

/****************************************************************************/
/*!
**  Correct for the mean in the case of non-conditional simualtions
**
** \param[in]  dbout     Output Db structure
** \param[in]  model     Model structure
** \param[in]  icase     Rank of PGS or GRF
** \param[in]  nbsimu    Number of simulations
**
*****************************************************************************/
static void st_mean_correct(Db    *dbout,
                            Model *model,
                            int    icase,
                            int    nbsimu)
{
  int isimu,ivar,ecr,iech;

  /* Loop on the simulations */

  for (isimu=ecr=0; isimu<nbsimu; isimu++)
  {

    /* Loop on the variables */

    for (ivar=0; ivar<model->getVariableNumber(); ivar++,ecr++)
    {

      /* Loop on the samples */

      for (iech=0; iech<dbout->getSampleNumber(); iech++)
      {
        if (! dbout->isActive(iech)) continue;
        dbout->updSimvar(LOC_SIMU, iech, isimu, ivar, icase, nbsimu,
                         model->getVariableNumber(), 0,
                         model->getContext().getMean(ivar));
      }  
    }
  }
}

/*****************************************************************************/
/*!
**  Perform non-conditional simulations on a set of gradient points using
**  Turning Bands method.
**
** \return  Error return code :
** \return    0 no problem
** \return    1 a structure cannot be simulated
** \return    2 no structure to be simulated 1 core problem
**
** \param[in]  dbgrd      Gradient Db structure
** \param[in]  model      Model structure
** \param[in]  situba     Situba structure
** \param[in]  aic        Array 'aic'
** \param[in]  delta      Value of the increment
**
** \remarks The simulated gradients are stored as follows:
** \remarks idim * nbsimu + isimu (for simulation at first point)
** \remarks idim * nbsimu + isimu + ndim * nbsimu (for simulation at 2nd point)
** \remarks At the end, the simulated gradient is stored at first point
**
*****************************************************************************/
static int st_simulate_gradient(Db     *dbgrd,
                                Model  *model,
                                Situba *situba,
                                double *aic,
                                double  delta)
{
  int icase,nbsimu,jsimu,ndim,error;
  double value1,value2;

  /* Initializations */

  error  = 1;
  icase  = 0;
  nbsimu = situba->nbsimu;
  ndim   = dbgrd->getNDim();

  for (int idim=0; idim<ndim; idim++)
  {

    /* Simulation at the initial location */

    for (int isimu=0; isimu<nbsimu; isimu++)
    {
      jsimu = isimu + idim * nbsimu;
      if (st_simulate_point(dbgrd,DATA,model,situba,aic,icase,jsimu)) 
        goto label_end;
    }

    /* Shift the information */

    for (int iech=0; iech<dbgrd->getSampleNumber(); iech++)
      if (dbgrd->isActive(iech))
        dbgrd->setCoordinate(iech,idim,dbgrd->getCoordinate(iech,idim) + delta);

    /* Simulation at the shift location */

    for (int isimu=0; isimu<nbsimu; isimu++)
    {
      jsimu = isimu + idim * nbsimu + ndim * nbsimu;
      if (st_simulate_point(dbgrd,DATA,model,situba,aic,icase,jsimu)) 
        goto label_end;
    }

    /* Un-Shift the information */

    for (int iech=0; iech<dbgrd->getSampleNumber(); iech++)
      if (dbgrd->isActive(iech))
        dbgrd->setCoordinate(iech,idim,dbgrd->getCoordinate(iech,idim) - delta);

    /* Scaling */

    for (int isimu=0; isimu<nbsimu; isimu++)
      for (int iech=0; iech<dbgrd->getSampleNumber(); iech++)
      {
        if (! dbgrd->isActive(iech)) continue;
        jsimu = isimu + idim * nbsimu + ndim * nbsimu;
        value2 = dbgrd->getSimvar(LOC_SIMU,iech,jsimu,0,icase,
                            2*ndim*nbsimu,1);
        jsimu = isimu + idim * nbsimu;
        value1 = dbgrd->getSimvar(LOC_SIMU,iech,jsimu,0,icase,
                            2*ndim*nbsimu,1);
        dbgrd->setSimvar(LOC_SIMU,iech,jsimu,0,icase,
                   2*ndim*nbsimu,1,(value2 - value1) / delta);
      }
  }

  /* Set the error return code */

  error = 0;

label_end:
  return(error);
}

/*****************************************************************************/
/*!
**  Perform non-conditional simulations on a set of tangent points using
**  Turning Bands method.
**
** \return  Error return code :
** \return    0 no problem
** \return    1 a structure cannot be simulated
** \return    2 no structure to be simulated 1 core problem
**
** \param[in]  dbtgt      Tangent Db structure
** \param[in]  model      Model structure
** \param[in]  situba     Situba structure
** \param[in]  aic        Array 'aic'
** \param[in]  delta      Value of the increment
**
** \remarks Warning: To perform the simulation of the tangent, we must
** \remarks simulated the gradients first. So we need to dimension the
** \remarks simulation outcome variables as for the gradients
**
*****************************************************************************/
static int st_simulate_tangent(Db     *dbtgt,
                               Model  *model,
                               Situba *situba,
                               double *aic,
                               double  delta)
{
  int error,nvar,nbsimu,icase;
  double value;

  /* Initializations */

  error  = 1;
  icase  = 0;
  nvar   = model->getVariableNumber();
  nbsimu = situba->nbsimu;

  /* Perform the simulation of the gradients at tangent points */

  if (st_simulate_gradient(dbtgt,model,situba,aic,delta)) goto label_end;

  /* Calculate the simulated tangent */

  for (int isimu=0; isimu<nbsimu; isimu++)
    for (int iech=0; iech<dbtgt->getSampleNumber(); iech++)
    {
      if (! dbtgt->isActive(iech)) continue;
      
      value = 0.;
      for (int idim=0; idim<dbtgt->getNDim(); idim++)
        value += dbtgt->getTangent(iech,idim) *
          dbtgt->getSimvar(LOC_SIMU,iech,isimu,0,icase,nbsimu,nvar);
      dbtgt->setSimvar(LOC_SIMU,iech,isimu,0,icase,nbsimu,nvar,value);
    }

  /* Set the error return code */

  error = 0;

label_end:
  return(error);
}

/****************************************************************************/
/*!
**  Check/Show the data (gaussian) against the closest grid node
**
** \param[in]  dbin       Input Db structure
** \param[in]  dbout      Output Db grid structure
** \param[in]  model      Model structure
** \param[in]  nbsimu     Number of simulations
** \param[in]  flag_pgs   1 if called from PGS
** \param[in]  flag_gibbs 1 if called from Gibbs
**
** \remark Attributes LOC_SIMU and LOC_GAUSFAC (for PGS) are mandatory
** \remark Tests have only been produced for icase=0
**
*****************************************************************************/
static void st_check_gaussian_data2grid(Db     *dbin,
                                        Db     *dbout,
                                        Model  *model,
                                        int     nbsimu,
                                        int     flag_pgs,
                                        int     flag_gibbs)
{
  int     iech,jech,isimu,number;
  double  valdat,valres,eps;
  VectorDouble coor;

  /* Initializations */

  if (dbin == (Db *) NULL || dbout == (Db *) NULL) return;
  number = 0;
  check_mandatory_attribute("st_check_gaussian_data2grid",dbout,LOC_SIMU);
  mestitle(1,"Checking Gaussian of data against closest grid node");

  /* Core allocation */

  coor.resize(dbin->getNDim());

  /* Loop on the data */

  for (iech=0; iech<dbin->getSampleNumber(); iech++)
  {
    if (! dbin->isActive(iech)) continue;


    // Find the index of the closest grid node and derive tolerance
    jech = index_point_to_grid(dbin,iech,0,dbout,coor.data());
    if (jech < 0) continue;
    eps = model_calcul_stdev(model,dbin,iech,dbout,jech,0,2.);
    if (eps < 1.e-6) eps = 1.e-6;
    
    for (isimu=0; isimu<nbsimu; isimu++)
    {
      if (! flag_pgs)
      {
        if (flag_gibbs)
          valdat = dbin->getSimvar(LOC_GAUSFAC,iech,isimu,0,0,nbsimu,1);
        else
          valdat = dbin->getVariable(iech,0);
      }
      else
      {
        valdat = dbin->getSimvar(LOC_GAUSFAC,iech,0,0,0,nbsimu,1);
      }

      valres = dbout->getSimvar(LOC_SIMU,jech,isimu,0,0,nbsimu,1);
      if (ABS(valdat - valres) < eps) continue;
      number++;
      
      /* The data facies is different from the grid facies */
      
      message("Inconsistency for Simulation (%d) between :\n",isimu+1);
      message("- Value (%lf) at Data (#%d) ",valdat,iech+1);
      message("at (");
      for (int idim=0; idim<dbin->getNDim(); idim++)
        message(" %lf",dbin->getCoordinate(iech,idim));
      message(")\n");

      message("- Value (%lf) at Grid (#%d) ",valres,jech+1);
      message("at (");
      for (int idim=0; idim<dbout->getNDim(); idim++)
        message(" %lf",dbout->getCoordinate(jech,idim));
      message(")\n");

      message("- Tolerance = %lf\n",eps);
    }
  }
  
  if (number <= 0) message("No problem found\n");
  return;
}

/****************************************************************************/
/*!
**  Perform the Simulation Process using the Turning Bands Method
**
** \return  Error return code
**
** \param[in]  dbin       Input Db structure
** \param[in]  dbout      Output Db structure
** \param[in]  model      Model structure
** \param[in]  neigh      Neigh structure
** \param[in]  situba     Situba structure
** \param[in]  dmean      Array giving the prior means for the drift terms
** \param[in]  dcov       Array containing the prior covariance matrix
**                        for the drift terms
** \param[in]  nbsimu     Number of simulations
** \param[in]  icase      Case for PGS or -1
** \param[in]  flag_pgs   1 if called from PGS
** \param[in]  flag_gibbs 1 if called from Gibbs
** \param[in]  flag_check 1 to check the proximity in Gaussian scale
**
*****************************************************************************/
static int st_simtub_process(Db     *dbin,
                             Db     *dbout,
                             Model  *model,
                             Neigh  *neigh,
                             Situba *situba,
                             double *dmean,
                             double *dcov,
                             int     nbsimu,
                             int     icase,
                             int     flag_pgs,
                             int     flag_gibbs,
                             int     flag_check)
{
  int     flag_cond,error,ncova,nvar;
  double *aic,*valpro,*vecpro;
  char    string[100];

  /* Initializations */

  error = 1;
  ncova = model->getCovaNumber();
  nvar  = model->getVariableNumber();
  aic   = valpro = vecpro = (double *) NULL;
  flag_cond  = (dbin != (Db *) NULL) && (neigh != (Neigh *) NULL);
  st_gendir(dbout,model,situba);
  st_minmax(dbout,model,situba);
  st_minmax(dbin ,model,situba);
  if (st_initialize(dbin,dbout,model,situba)) goto label_end;

  /* Calculate the 'aic' array */

  valpro = (double *) mem_alloc(sizeof(double) * nvar,1);
  vecpro = (double *) mem_alloc(sizeof(double) * nvar  * nvar,1);
  aic    = (double *) mem_alloc(sizeof(double) * nvar  * nvar  * ncova,1);
  if (model_update_coreg(model,aic,valpro,vecpro))
    messageAbort("model_update_coreg");

  /* Non conditional simulations on the data points */

  if (flag_cond)
  {
    if (st_simulate_point(dbin,DATA,model,situba,aic,icase,0)) 
      goto label_end;

    /* Calculate the simulated error */
    
    st_difference(dbin,situba,nvar,icase,flag_pgs,flag_gibbs);
  }

  /* Non conditional simulations on the grid */
  
  if (is_grid(dbout))
  {
    if (st_simulate_grid(dbout,model,situba,aic,icase,0)) 
      goto label_end;
  }
  else
  {
    if (st_simulate_point(dbout,RESULT,model,situba,aic,icase,0))
      goto label_end;
  }

  /* Add the contribution of nugget effect (optional) */

  st_simulate_nugget(dbout,RESULT,model,situba,aic,icase);

  /* Conditional simulations */

  if (flag_cond)
  {
    st_title_process(string,"Conditioning by Kriging");
    if (krigsim(string,dbin,dbout,model,neigh,dmean,dcov,icase,
                nbsimu,FLAG_DGM,R_COEFF)) goto label_end;
  }
  else
  {

    /* In non-conditional case, correct for the mean */

    st_mean_correct(dbout,model,icase,nbsimu);
  }

  /* Copy value from data to coinciding grid node */

  if (flag_cond)
    st_update_data2target(dbin,dbout,model,situba,icase,flag_pgs);

  /* Check consistency between data and resulting simulations (optional) */

  if (flag_check)
    st_check_gaussian_data2grid(dbin,dbout,model,nbsimu,flag_pgs,flag_gibbs);
    
  /* Set the error return code */

  error = 0;

label_end:
  aic    = (double *) mem_free((char *) aic);
  valpro = (double *) mem_free((char *) valpro);
  vecpro = (double *) mem_free((char *) vecpro);
  return(error);
}

/****************************************************************************/
/*!
**  Perform the (non-conditional) Simulation(s) using the Turning Bands Method
**
** \return  Error return code
**
** \param[in]  dbiso     Isovalues Db structure
** \param[in]  dbgrd     Gradient Db structure
** \param[in]  dbtgt     Tangent Db structure
** \param[in]  dbout     Output Db structure
** \param[in]  model     Model structure
** \param[in]  nbsimu    Number of simulations
** \param[in]  nbtuba    Number of turning bands
** \param[in]  delta     Value of the increment
**
*****************************************************************************/
GEOSLIB_API int simtub_potential(Db     *dbiso,
                                 Db     *dbgrd,
                                 Db     *dbtgt,
                                 Db     *dbout,
                                 Model  *model,
                                 int     nbsimu,
                                 int     nbtuba,
                                 double  delta)
{
  Situba *situba;
  int     error,ncova,nvar,icase;
  double *aic,*valpro,*vecpro;

  /* Initializations */

  error = 1;
  icase = 0;
  ncova = model->getCovaNumber();
  nvar  = model->getVariableNumber();
  aic   = valpro = vecpro = (double *) NULL;
  situba = (Situba *) NULL;

  /* Processing the Turning Bands algorithm */

  situba = st_alloc(model,nbsimu,nbtuba);
  if (situba == (Situba *) NULL) goto label_end;

  st_gendir(dbout,model,situba);
  st_minmax(dbout,model,situba);
  st_minmax(dbiso,model,situba);
  st_minmax(dbgrd,model,situba);
  st_minmax(dbtgt,model,situba);
  if (st_initialize(dbiso,dbout,model,situba)) goto label_end;

  /* Calculate the 'aic' array */

  valpro = (double *) mem_alloc(sizeof(double) * nvar,1);
  vecpro = (double *) mem_alloc(sizeof(double) * nvar  * nvar,1);
  aic    = (double *) mem_alloc(sizeof(double) * nvar  * nvar  * ncova,1);
  if (model_update_coreg(model,aic,valpro,vecpro))
    messageAbort("model_update_coreg");

  /* Non conditional simulations on the data points */

  if (dbiso != (Db *) NULL)
  {
    if (st_simulate_point(dbiso,DATA,model,situba,aic,icase,0)) goto label_end;
  }

  /* Non conditional simulations on the gradient points */

  if (dbgrd != (Db *) NULL) 
  {
    if (st_simulate_gradient(dbgrd,model,situba,aic,delta)) goto label_end;
  }

  /* Non conditional simulations on the tangent points */

  if (dbtgt != (Db *) NULL) 
  {
    if (st_simulate_tangent(dbtgt,model,situba,aic,delta)) goto label_end;
  }

  /* Non conditional simulations on the grid */

  if (is_grid(dbout))
  {
    if (st_simulate_grid(dbout,model,situba,aic,icase,0)) goto label_end;
  }
  else
  {
    if (st_simulate_point(dbout,RESULT,model,situba,aic,icase,0)) 
      goto label_end;
  }

  /* Add the contribution of nugget effect (optional) */

  st_simulate_nugget(dbout,RESULT,model,situba,aic,icase);

  /* Set the error return code */

  error = 0;

label_end:
  situba = st_dealloc(model,situba);
  aic    = (double *) mem_free((char *) aic);
  valpro = (double *) mem_free((char *) valpro);
  vecpro = (double *) mem_free((char *) vecpro);
  return(error);
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
** \param[in]  neigh      Neigh structure (optional)
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
GEOSLIB_API int simtub(Db *dbin,
                       Db *dbout,
                       Model *model,
                       Neigh *neigh,
                       int nbsimu,
                       int seed,
                       int nbtuba,
                       int flag_check,
                       NamingConvention namconv)
{
  Situba *situba;
  int     flag_cond,nvar,error,iext,inostat,iptr_in,iptr_out;

  /* Initializations */

  error = 1;
  nvar  = model->getVariableNumber();
  iptr_in = iptr_out = -1;
  flag_cond = (dbin != (Db *) NULL && neigh != (Neigh *) NULL);
  situba = (Situba *) NULL;
  law_set_random_seed(seed);
  if (st_check_simtub_environment(dbin,dbout,model,neigh)) goto label_end;
  if (manage_external_info(1,LOC_F,dbin,dbout,&iext)) goto label_end;
  if (manage_external_info(1,LOC_NOSTAT,dbin,dbout,&inostat)) goto label_end;

  /* Define the environment variables for printout */

  st_simulation_environment();
  MES_NPGS = MES_NGRF = 1;

  /* Add the attributes for storing the results */

  if (flag_cond)
  {
    if (db_locator_attribute_add(dbin, LOC_SIMU, nvar * nbsimu, 0, 0.,
                                 &iptr_in)) goto label_end;
  }
  if (db_locator_attribute_add(dbout, LOC_SIMU, nvar * nbsimu, 0, 0.,
                               &iptr_out)) goto label_end;

  /* Processing the Turning Bands algorithm */

  situba = st_alloc(model,nbsimu,nbtuba);
  if (situba == (Situba *) NULL) goto label_end;
  if (st_simtub_process(dbin,dbout,model,neigh,situba,
                        (double *) NULL,(double *) NULL,
                        nbsimu,0,0,0,flag_check)) goto label_end;

  /* Free the temporary variables */

  if (flag_cond) dbin->deleteFieldByLocator(LOC_SIMU);

  /* Set the error return flag */

  error = 0;
  namconv.setNamesAndLocators(dbin,LOC_Z,model->getVariableNumber(),
                              dbout,iptr_out,String(),nbsimu);

label_end:
  (void) manage_external_info(-1,LOC_F,dbin,dbout,&iext);
  (void) manage_external_info(-1,LOC_NOSTAT,dbin,dbout,&inostat);
  situba = st_dealloc(model,situba);
  return(error);
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
** \param[in]  neigh      Neigh structure (optional)
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
GEOSLIB_API int simdgm(Db    *dbin,
                       Db    *dbout,
                       Model *model,
                       Neigh *neigh,
                       double rval,
                       int    seed,
                       int    nbsimu,
                       int    nbtuba,
                       int    flag_check)
{
  Situba *situba;
  int     flag_cond,nvar,error,iext,inostat,iptr;

  /* Initializations */

  error = 1;
  nvar  = model->getVariableNumber();
  iptr = -1;
  flag_cond = (dbin != (Db *) NULL && neigh != (Neigh *) NULL);
  situba = (Situba *) NULL;
  law_set_random_seed(seed);
  if (st_check_simtub_environment(dbin,dbout,model,neigh)) goto label_end;
  if (manage_external_info(1,LOC_F,dbin,dbout,&iext)) goto label_end;
  if (manage_external_info(1,LOC_NOSTAT,dbin,dbout,
                           &inostat)) goto label_end;

  /* Define the environment variables for printout */

  st_simulation_environment();
  MES_NPGS = MES_NGRF = 1;
  FLAG_DGM = 1;
  R_COEFF = rval;

  /* Add the attributes for storing the results */

  if (flag_cond)
  {
    if (db_locator_attribute_add(dbin,LOC_SIMU,nvar*nbsimu,0,0.,
                                 &iptr)) goto label_end;
  }
  if (db_locator_attribute_add(dbout,LOC_SIMU,nvar*nbsimu,0,0.,
                               &iptr)) goto label_end;

  /* Processing the Turning Bands algorithm */

  situba = st_alloc(model,nbsimu,nbtuba);
  if (situba == (Situba *) NULL) goto label_end;
  if (st_simtub_process(dbin,dbout,model,neigh,situba,
                        (double *) NULL,(double *) NULL,
                        nbsimu,0,0,0,flag_check)) goto label_end;

  /* Free the temporary variables */

  if (flag_cond) dbin->deleteFieldByLocator(LOC_SIMU);

  /* Set the error return flag */

  error = 0;

label_end:
  (void) manage_external_info(-1,LOC_F,dbin,dbout,&iext);
  (void) manage_external_info(-1,LOC_NOSTAT,dbin,dbout,&inostat);
  situba = st_dealloc(model,situba);
  return(error);
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
** \param[in]  neigh      Neigh structure (optional)
** \param[in]  dmean      Array giving the prior means for the drift terms
** \param[in]  dcov       Array containing the prior covariance matrix
**                        for the drift terms
** \param[in]  seed       Seed for random number generator
** \param[in]  nbsimu     Number of simulations
** \param[in]  nbtuba     Number of turning bands
** \param[in]  flag_check 1 to check the proximity in Gaussian scale
**
** \remark  The arguments 'dbout' and 'neigh' are optional: they must
** \remark  be defined only for conditional simulations
**
*****************************************************************************/
GEOSLIB_API int simbayes(Db     *dbin,
                         Db     *dbout,
                         Model  *model,
                         Neigh  *neigh,
                         double *dmean,
                         double *dcov,
                         int     seed,
                         int     nbsimu,
                         int     nbtuba,
                         int     flag_check)
{
  Situba *situba;
  int     flag_cond,nvar,error,iptr;

  /* Initializations */

  error = 1;
  nvar  = model->getVariableNumber();
  iptr  = -1;
  flag_cond = (dbin != (Db *) NULL && neigh != (Neigh *) NULL);
  situba = (Situba *) NULL;
  law_set_random_seed(seed);
  if (st_check_simtub_environment(dbin,dbout,model,neigh)) goto label_end;

  /* Define the environment variables for printout */

  st_simulation_environment();
  MES_NPGS = MES_NGRF = 1;

  /* Add the attributes for storing the results */

  if (flag_cond)
  {
    if (db_locator_attribute_add(dbin,LOC_SIMU,nvar*nbsimu,0,0,
                                 &iptr)) goto label_end;
  }
  if (db_locator_attribute_add(dbout,LOC_SIMU,nvar*nbsimu,0,0,
                               &iptr)) goto label_end;

  /* Processing the Turning Bands algorithm */

  situba = st_alloc(model,nbsimu,nbtuba);
  if (situba == (Situba *) NULL) goto label_end;
  if (st_simtub_process(dbin,dbout,model,neigh,situba,dmean,dcov,
                        nbsimu,0,0,0,flag_check)) goto label_end;

  /* Free the temporary variables */

  if (flag_cond) dbin->deleteFieldByLocator(LOC_SIMU);

  /* Set the error return flag */

  error = 0;

label_end:
  situba = st_dealloc(model,situba);
  return(error);
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
GEOSLIB_API int get_rank_from_propdef(PropDef *propdef,
                                      int    ipgs,
                                      int    igrf)
{
  if (ipgs <= 0 || propdef == (PropDef *) NULL)
    return(igrf);
  else
    return(propdef->ngrf[0] + igrf);
}

/****************************************************************************/
/*!
**  Suppresses the added samples
**
** \param[in]  db      Db structure
** \param[in]  nech    initial number of samples
**
*****************************************************************************/
static void st_suppress_added_samples(Db *db,
                                      int nech)
{
  int iech;

  if (nech <= 0) return;
  for (iech=db->getSampleNumber()-1; iech>=nech; iech--)
    db->deleteSample(iech);
  return;
}

/****************************************************************************/
/*!
**  Check/Show the data against facies at the closest grid node
**
** \param[in]  propdef    PropDef structure
** \param[in]  dbin       Input Db structure
** \param[in]  dbout      Output Db grid structure
** \param[in]  rule       Lithotype Rule definition
** \param[in]  flag_used  Tell if a GRF is used or not
** \param[in]  flag_stat  1 for stationary; 0 otherwise
** \param[in]  flag_check 1 check the consistency between data and grid
** \param[in]  flag_show  1 show the data on grid
** \param[in]  ipgs       Rank of the PGS
** \param[in]  nechin     Initial number of data
** \param[in]  nvar       Number of variables
** \param[in]  nfacies    Number of facies
** \param[in]  nbsimu     Number of simulations
**
** \remark Attributes LOC_FACIES are mandatory
** \remark Attributes LOC_GAUSFAC are mandatory
**
*****************************************************************************/
static void st_check_facies_data2grid(PropDef  *propdef,
                                      Db     *dbin,
                                      Db     *dbout,
                                      const Rule   *rule,
                                      int    *flag_used,
                                      int     flag_stat,
                                      int     flag_check,
                                      int     flag_show,
                                      int     ipgs,
                                      int     nechin,
                                      int     nvar,
                                      int     nfacies,
                                      int     nbsimu)
{
  int     iech,jech,isimu,facdat,facres,number;
  double *coor;

  /* Initializations */

  check_mandatory_attribute("st_check_facies_data2grid",dbout,LOC_FACIES);
  number = 0;
  coor   = (double *) NULL;
  if (flag_check) 
    mestitle(1,"Checking facies of data against closest grid node (PGS=%d)",
             ipgs+1);

  /* Core allocation */

  coor = db_sample_alloc(dbin,LOC_X);
  if (coor == (double *) NULL) goto label_end;

  /* Loop on the data */

  for (iech=0; iech<nechin; iech++)
  {
    if (! dbin->isActive(iech)) continue;
    facdat = (int) dbin->getVariable(iech,0);
    if (facdat < 1 || facdat > nfacies) continue;
    jech = index_point_to_grid(dbin,iech,0,dbout,coor);
    if (jech < 0) continue;
    
    for (isimu=0; isimu<nbsimu; isimu++)
    {
      facres = (int) dbout->getSimvar(LOC_FACIES, jech, isimu, 0, ipgs, nbsimu,
                                      1);
      if (flag_show)
      {
        if (facdat == facres)
          dbout->setSimvar(LOC_FACIES,jech,isimu,0,ipgs,nbsimu,1,-facdat);
        else
          dbout->setSimvar(LOC_FACIES,jech,isimu,0,ipgs,nbsimu,1,0.);
      }
      
      if (facdat == facres) continue;
      number++;
      
      /* The data facies is different from the grid facies */
      
      if (flag_check) 
      {
        message("Inconsistency for Simulation (%d) between :\n",isimu+1);
        message("- Facies (%d) at Data (#%d)\n",facdat,iech+1);
        message("- Facies (%d) at Grid (#%d)\n",facres,jech+1);
      }
    }
  }
  
label_end:
  if (flag_check && number <= 0) message("No problem found\n");
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
** \param[in]  neigh       Neighborhood structure
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
** \param[in]  gibbs_eps   Relative immobile criterion
** \param[in]  delta       Spatial increment (used for generating increments)
** \param[in]  namconv     Naming convention
**
** \remark  When conditional, the unique variable in the input Db structure
** \remark  should correspond to the facies index (starting from 1)
** \remark  The argument 'dbin' is optional: it must be defined only for
** \remark  conditional simulations
**
*****************************************************************************/
GEOSLIB_API int simpgs(Db *dbin,
                       Db *dbout,
                       RuleProp* ruleprop,
                       Model *model1,
                       Model *model2,
                       Neigh *neigh,
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
                       double gibbs_eps,
                       double delta,
                       NamingConvention namconv)
{
  int     iptr,igrf,icase,nfacies,flag_used[2];
  int     iptr_RP,iptr_RF,iptr_DF,iptr_DN,iptr_RN;
  Situba *situba;
  Model  *models[2];
  PropDef  *propdef;
  std::vector<Model *> modvec;

  /* Initializations */

  int error     = 1;
  int nvar      = 1;
  int nechin    = 0;
  int ngrf      = 0;
  situba    = (Situba *) NULL;
  propdef   = (PropDef *) NULL;
  models[0] = model1;
  models[1] = model2;
  bool flag_cond = (dbin != (Db *) NULL);
  iptr_RP   = iptr_RF = iptr_DF = nfacies = 0;
  iptr      = -1;
  bool    verbose = false;
  law_set_random_seed(seed);

  if (ruleprop == nullptr)
  {
    messerr("RuleProp must be defined");
    return 1;
  }
  int flag_stat = ruleprop->isFlagStat();
  const Rule* rule = ruleprop->getRule();
  const VectorDouble& propcst = ruleprop->getPropCst();
  const Db* dbprop = ruleprop->getDbprop();

  ngrf = rule->getGRFNumber();
  if (rule->particularities(dbout,dbprop,model1,1,flag_stat))
    goto label_end;
  if (st_check_simtub_environment(dbin,dbout,model1,neigh)) goto label_end;

  /**********************/
  /* Preliminary checks */
  /**********************/

  /* Input Db */
  if (flag_cond)
  {
    nechin = dbin->getSampleNumber();
    if (! dbin->isVariableNumberComparedTo(1)) goto label_end;
  }

  /* Output Db */
  if (flag_modif && flag_gaus)
  {
    messerr("Calculating the facies proportions is incompatible with storing the Gaussian values");
    goto label_end;
  }

  /* Model */
  for (int igrf=0; igrf<2; igrf++)
  {
    flag_used[igrf] = rule->isYUsed(igrf);
    if (! flag_used[igrf]) continue;
    if (models[igrf] == (Model *) NULL)
    {
      messerr("The Underlying GRF #%d is needed",igrf+1);
      messerr("No corresponding Model is provided");
      goto label_end;
    }
    if (models[igrf]->getVariableNumber() != 1)
    {
      messerr("The number of variables in the model #%d (%d) should be 1",
              igrf+1,model1->getVariableNumber());
      goto label_end;
    }
    if (model_stabilize(models[igrf],1,percent)) goto label_end;
    if (model_normalize(models[igrf],1)) goto label_end;
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

  MES_NPGS = 1;
  MES_IPGS = 0;
  MES_NGRF = ngrf;

  /**********************/
  /* Add the attributes */
  /**********************/

  nfacies = rule->getFaciesNumber();

  /* Storage of the facies proportions */
  if (flag_modif)
  {
    if (db_locator_attribute_add(dbout,LOC_P,nfacies,0,0.,&iptr_RP)) 
      goto label_end;
  }

  /* Storage of the facies simulations in the Output Db */
  if (db_locator_attribute_add(dbout,LOC_FACIES,nbsimu,0,0.,&iptr_RF)) 
    goto label_end;

  /* Storage of the facies simulations in the input file */
  if (flag_cond)
  {
    if (db_locator_attribute_add(dbin,LOC_FACIES,nbsimu,0,0.,&iptr_DF)) 
      goto label_end;
  }

  if (flag_cond)
  {
    /* Gaussian transform of the facies input data */
    if (db_locator_attribute_add(dbin,LOC_GAUSFAC,ngrf*nbsimu,
                                 0,0.,&iptr)) goto label_end;

    /* Non-conditional simulations at data points */
    if (db_locator_attribute_add(dbin,LOC_SIMU,ngrf*nbsimu,
                                 0,0.,&iptr_DN)) goto label_end;
  }

  /* (Non-) Conditional simulations at target points */
  if (db_locator_attribute_add(dbout,LOC_SIMU,ngrf*nbsimu,
                               0,0.,&iptr_RN)) goto label_end;

  if (flag_cond)
  {
    /* Lower bound at input data points */
    if (db_locator_attribute_add(dbin,LOC_L,ngrf,0,0.,&iptr)) goto label_end;

    /* Upper bound at input data points */
    if (db_locator_attribute_add(dbin,LOC_U,ngrf,0,0.,&iptr)) goto label_end;
  }

  propdef = proportion_manage(1,1,flag_stat,ngrf,0,nfacies,0,dbin,dbprop,
                              propcst,propdef);
  if (propdef == (PropDef *) NULL) goto label_end;
  simu_define_func_update(simu_func_categorical_update);
  simu_define_func_scale (simu_func_categorical_scale);
  ModCat.propdef = propdef;
  ModCat.rule    = rule;
  ModCat.ipgs    = 0;
  ModCat.flag_used[0] = flag_used[0];
  ModCat.flag_used[1] = flag_used[1];

  /****************************************/
  /* Convert facies into gaussian at data */
  /****************************************/

  proportion_rule_process(propdef,PROCESS_COPY);

  if (flag_cond)
  {
    int npgs = 1;
    int ipgs = 0;

    // Create the Gibbs sampler (multi-mono case)
    AGibbs* gibbs = GibbsFactory::createGibbs(dbin, modvec, rule->getRho(), false);

    /* Initialize the Gibbs calculations */

    gibbs->init(npgs, ngrf, gibbs_nburn, gibbs_niter, 0, true, true,
                rule->getRho());
      
    /* Allocate the covariance matrix inverted */
  
    if (gibbs->covmatAlloc(verbose)) goto label_end;

    /* Allocate the Gaussian vector */

    VectorVectorDouble y = gibbs->allocY();

    /* Loop on the simulations */

    for (int isimu=0; isimu<nbsimu; isimu++)
    {
      for (int igrf = 0; igrf < ngrf; igrf++)
        if (rule->evaluateBounds(propdef, dbin, dbout, isimu, igrf, ipgs,
                                 nbsimu)) goto label_end;
      
      /* Initialization for the Gibbs sampler */

      if (gibbs->calculInitialize(y, isimu, ipgs, verbose)) goto label_end;
      if (verbose) gibbs->print(true,y,isimu,ipgs);

      /* Iterations of the Gibbs sampler */

      for (int iter = 0; iter < gibbs->getNiter(); iter++)
        gibbs->update(y, isimu, ipgs,  iter);
        
      /* Check the validity of the Gibbs results (optional) */
        
      if (flag_check) gibbs->checkGibbs(y,isimu,ipgs);
      if (verbose) gibbs->print(false, y, isimu, ipgs);

      // Store the results

      gibbs->storeResult(y, isimu, ipgs);
    }
  }

  /***************************************************/
  /* Perform the conditional simulation for each GRF */
  /***************************************************/

  for (MES_IGRF=0; MES_IGRF<2; MES_IGRF++)
  {
    if (! flag_used[MES_IGRF]) continue;
    icase  = get_rank_from_propdef(propdef,0,MES_IGRF);
    situba = st_alloc(models[MES_IGRF],nbsimu,nbtuba);
    if (situba == (Situba *) NULL) goto label_end;
    if (st_simtub_process(dbin,dbout,models[MES_IGRF],neigh,situba,
                          (double *) NULL,(double *) NULL,
                          nbsimu,icase,1,0,flag_check)) goto label_end;
    situba = st_dealloc(models[MES_IGRF],situba);
  }
  
  /* Convert gaussian to facies at target point */
  
  if (! flag_gaus)
  {
    for (int isimu=0; isimu<nbsimu; isimu++)
      simu_func_categorical_transf(dbout,0,isimu,nbsimu);
  }
  
  /* Update facies proportions at target points */
  
  if (flag_modif)
  {
    for (int isimu=0; isimu<nbsimu; isimu++)
      simu_func_categorical_update(dbout,0,isimu,nbsimu);
    simu_func_categorical_scale(dbout,0,nbsimu);
  }
  
  /* Check/show facies at data against facies at the closest grid node */
  
  if (flag_cond && ! flag_gaus && (flag_check || flag_show))
    st_check_facies_data2grid(propdef,dbin,dbout,rule,flag_used,
                              flag_stat,flag_check,flag_show,0,
                              nechin,nvar,nfacies,nbsimu);

  /********************************/
  /* Free the temporary variables */
  /********************************/

  if (dbout != nullptr)
  {
    if (!st_keep(flag_gaus, flag_modif, RESULT, PROP) && iptr_RP)
      dbout->deleteFieldByLocator(LOC_P);
    else
      namconv.setNamesAndLocators(NULL, LOC_Z, -1, dbout, iptr_RP, "Props",
                                  nbsimu, false);

    if (!st_keep(flag_gaus, flag_modif, RESULT, GAUS))
      dbout->deleteFieldByLocator(LOC_SIMU);
    else
      namconv.setNamesAndLocators(NULL, LOC_Z, -1, dbout, iptr_RN, "Gaus",
                                  ngrf * nbsimu, false);

    if (!st_keep(flag_gaus, flag_modif, RESULT, FACIES) && iptr_RF)
      dbout->deleteFieldByLocator(LOC_FACIES);
    else
      namconv.setNamesAndLocators(NULL, LOC_Z, -1, dbout, iptr_RF, String(),
                                  nbsimu);
  }

  if (dbin != nullptr)
  {
    if (!st_keep(flag_gaus, flag_modif, DATA, GAUS))
      dbin->deleteFieldByLocator(LOC_SIMU);
    else
      namconv.setNamesAndLocators(NULL, LOC_Z, -1, dbin, iptr_DN, "Gaus",
                                  ngrf * nbsimu, false);

    if (!st_keep(flag_gaus, flag_modif, DATA, FACIES) && iptr_DF)
      dbin->deleteFieldByLocator(LOC_FACIES);
    else
      namconv.setNamesAndLocators(NULL, LOC_Z, -1, dbin, iptr_DF, String(),
                                  nbsimu, false);

    dbin->deleteFieldByLocator(LOC_GAUSFAC);
    dbin->deleteFieldByLocator(LOC_L);
    dbin->deleteFieldByLocator(LOC_U);
  }

  /* Set the error return flag */

  error = 0;

label_end:
  propdef = proportion_manage(-1,1,flag_stat,ngrf,0,nfacies,0,dbin,dbprop,
                              propcst,propdef);
  st_suppress_added_samples(dbin,nechin);
  for (igrf=0; igrf<2; igrf++)
    situba = st_dealloc(models[igrf],situba);
  return(error);
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
** \param[in]  neigh       Neighborhood structure
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
** \param[in]  gibbs_eps   Relative stability criterion
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
GEOSLIB_API int simbipgs(Db       *dbin,
                         Db       *dbout,
                         RuleProp *ruleprop,
                         Model    *model11,
                         Model    *model12,
                         Model    *model21,
                         Model    *model22,
                         Neigh    *neigh,
                         int       nbsimu,
                         int       seed,
                         int       flag_gaus,
                         int       flag_modif,
                         int       flag_check,
                         int       flag_show,
                         int       nbtuba,
                         int       gibbs_nburn,
                         int       gibbs_niter,
                         double    percent,
                         double    gibbs_eps,
                         NamingConvention namconv)
{
  int     iptr,igrf,iatt_z[2];
  int     nvar,ipgs,npgs,flag_cond,error,isimu,icase;
  int     nfac[2],nfactot,flag_used[2][2],nechin,ngrf[2],ngrftot;
  int     iptr_RP,iptr_RF,iptr_DF,iptr_RN,iptr_DN;
  bool    verbose;
  Rule   *rules[2];
  Model  *models[2][2];
  std::vector<Model *> modvec[2];
  Situba *situba;
  PropDef  *propdef;

  /* Initializations */

  error     = 1;
  nvar      = 1;
  npgs      = 2;
  nechin    = 0;
  verbose   = false;
  situba    = (Situba *) NULL;
  propdef   = (PropDef  *) NULL;

  if (ruleprop == nullptr)
  {
    messerr("RuleProp must be defined");
    return 1;
  }
  if (ruleprop->getRule(0)->getModeRule() != RULE_STD)
  {
    messerr("SimuBiPgs is restricted to Standard Lithotype Rule");
    return 1;
  }
  int flag_stat = ruleprop->isFlagStat();
  Rule rule1(*ruleprop->getRule(0));
  Rule rule2(*ruleprop->getRule(1));
  const VectorDouble& propcst = ruleprop->getPropCst();
  const Db* dbprop = ruleprop->getDbprop();

  nfac[0]   = rule1.getFaciesNumber();
  nfac[1]   = rule2.getFaciesNumber();
  rules[0]  = &rule1;
  rules[1]  = &rule2;
  models[0][0] = model11;
  models[0][1] = model12;
  models[1][0] = model21;
  models[1][1] = model22;
  nfactot   = nfac[0] + nfac[1];
  flag_cond = (dbin != (Db *) NULL);
  iptr_RP   = iptr_RF = iptr_DF = iptr_RN = iptr_DN = 0;
  iptr      = -1;
  for (ipgs=0; ipgs<2; ipgs++)
  {
    ngrf[ipgs] = 0;
    iatt_z[ipgs] = -1;
  }
  law_set_random_seed(seed);
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
    if (! dbin->isVariableNumberComparedTo(2)) goto label_end;
    iatt_z[0] = db_attribute_identify(dbin,LOC_Z,0);
    iatt_z[1] = db_attribute_identify(dbin,LOC_Z,1);
  }

  /* Output Db */
  if (flag_modif && flag_gaus)
  {
    messerr("Calculating the facies proportions is incompatible with storing the Gaussian values");
    goto label_end;
  }

  /* Model */

  ngrftot = 0;
  for (ipgs=0; ipgs<npgs; ipgs++)
  {
    ngrf[ipgs] = rules[ipgs]->getGRFNumber();
    ngrftot += ngrf[ipgs];

    /* Check the validity of the model */

    for (igrf=0; igrf<2; igrf++)
    {
      flag_used[ipgs][igrf] = rules[ipgs]->isYUsed(igrf);
      if (! flag_used[ipgs][igrf]) continue;
      if (models[ipgs][igrf] == (Model *) NULL)
      {
        messerr("Variable #%d needs the underlying GRF #%d",ipgs+1,igrf+1);
        messerr("No corresponding Model is provided");
        goto label_end;
      }
      if (models[ipgs][igrf]->getVariableNumber() != 1)
      {
        messerr("The number of variables in Model #%d (%d) for Variable %d should be 1",
                igrf+1,ipgs+1,models[ipgs][igrf]->getVariableNumber());
        goto label_end;
      }
      if (model_stabilize(models[ipgs][igrf],1,percent)) goto label_end;
      if (model_normalize(models[ipgs][igrf],1)) goto label_end;

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

  for (ipgs=0; ipgs<npgs; ipgs++)
  {
    // Check the Rules (only RULE_STD case is authorized)
    if (rules[ipgs]->getModeRule() != RULE_STD)
    {
      messerr("In the Bi-PGS application, only Standard Rule is authorized");
    }
  }

  /* Final checks */

  for (ipgs=0; ipgs<2; ipgs++)
  {
    if (flag_cond)
    {
      dbin->clearLocators(LOC_Z);
      dbin->setLocatorByAttribute(iatt_z[ipgs],LOC_Z);
    }
    if (st_check_simtub_environment(dbin,dbout,models[ipgs][0],neigh))
      goto label_end;
  }

  /* Core allocation */

  propdef = proportion_manage(1,1,flag_stat,ngrf[0],ngrf[1],nfac[0],nfac[1],
                              dbin,dbprop,propcst,propdef);
  if (propdef == (PropDef *) NULL) goto label_end;
  simu_define_func_update(simu_func_categorical_update);
  simu_define_func_scale (simu_func_categorical_scale);
  ModCat.propdef = propdef;

  /**********************/
  /* Add the attributes */
  /**********************/

  /* Storage of the proportions */
  if (flag_modif)
  {
    if (db_locator_attribute_add(dbout,LOC_P,nfactot,0,0.,&iptr_RP)) 
      goto label_end;
  }

  /* Storage of the facies simulations in the output file */
  if (db_locator_attribute_add(dbout,LOC_FACIES,npgs*nbsimu,
                               0,0.,&iptr_RF)) goto label_end;

  /* Storage of the facies simulations in the input file */
  if (flag_cond)
  {
    if (db_locator_attribute_add(dbin,LOC_FACIES,npgs*nbsimu,
                                 0,0.,&iptr_DF)) goto label_end;
  }
  
  /* Gaussian transform of the facies input data */
  if (flag_cond)
  {
    if (db_locator_attribute_add(dbin,LOC_GAUSFAC,ngrftot*nbsimu,
                                 0,0.,&iptr)) goto label_end;
  }
  
  /* Non-conditional gaussian simulations at data points */
  if (flag_cond)
  {
    if (db_locator_attribute_add(dbin,LOC_SIMU,ngrftot * nbsimu,
                                     0,0.,&iptr_DN)) goto label_end;
  }
  
  /* Non-conditional gaussian simulations at target points */
  if (db_locator_attribute_add(dbout,LOC_SIMU,ngrftot * nbsimu,
                               0,0.,&iptr_RN)) goto label_end;

  if (flag_cond)
  {
    /* Lower bound at input data points */
    if (db_locator_attribute_add(dbin,LOC_L,ngrftot,0,0.,&iptr)) goto label_end;

    /* Upper bound at input data points */
    if (db_locator_attribute_add(dbin,LOC_U,ngrftot,0,0.,&iptr)) goto label_end;
  }

  /* Define the environment variables for printout */

  MES_NPGS = npgs;

  /************************/
  /* Main loop on the PGS */
  /************************/

  for (ipgs=0; ipgs<npgs; ipgs++)
  {
    if (flag_cond)
    {
      dbin->clearLocators(LOC_Z);
      dbin->setLocatorByAttribute(iatt_z[ipgs],LOC_Z);
    }

    if (ipgs == 0)
      proportion_rule_process(propdef,PROCESS_MARGINAL);
    else
      proportion_rule_process(propdef,PROCESS_CONDITIONAL);

    ModCat.rule  = rules[ipgs];
    ModCat.ipgs  = ipgs;
    ModCat.flag_used[0] = flag_used[ipgs][0];
    ModCat.flag_used[1] = flag_used[ipgs][1];

    /****************************************/
    /* Convert facies into gaussian at data */
    /****************************************/

    if (flag_cond)
    {

      // Create the Gibbs sampler

      AGibbs* gibbs = GibbsFactory::createGibbs(dbin, modvec[ipgs], rules[ipgs]->getRho(), false);

      /* Initialize the Gibbs calculations */

      gibbs->init(npgs, ngrf[ipgs], gibbs_nburn, gibbs_niter,
                  0, true, true, rules[ipgs]->getRho());

      /* Allocate the covariance matrix inverted */

      if (gibbs->covmatAlloc(0)) goto label_end;

      // Core allocation

      VectorVectorDouble y = gibbs->allocY();

      /* Loop on the simulations */

      for (isimu=0; isimu<nbsimu; isimu++)
      {
	  
        /* Update the proportions */

        for (int igrf = 0; igrf < ngrf[ipgs]; igrf++)
          if (rules[ipgs]->evaluateBounds(propdef, dbin, dbout, isimu, igrf,
                                          ipgs, nbsimu)) goto label_end;

        /* Initialization for the Gibbs sampler */

        if (gibbs->calculInitialize(y, isimu, ipgs, verbose)) goto label_end;
        if (verbose) gibbs->print(true,y,isimu,ipgs);

        /* Iterations of the Gibbs sampler */

        for (int iter = 0; iter < gibbs->getNiter(); iter++)
          gibbs->update(y, isimu, ipgs, iter);

        /* Check the validity of the Gibbs results (optional) */

        if (flag_check) gibbs->checkGibbs(y,isimu,ipgs);
        if (verbose) gibbs->print(false,y,isimu,ipgs);

        // Store the results

        gibbs->storeResult(y, isimu, ipgs);
      }
      
      /* Convert gaussian to facies on data point */
      
      for (int isimu=0; isimu<nbsimu; isimu++)
      {
        if (rules[ipgs]->gaus2facData(propdef,dbin,dbout,flag_used[ipgs],
                                      ipgs,isimu,nbsimu)) goto label_end;
      }
    }
    
    /***************************************************/
    /* Perform the conditional simulation for each GRF */
    /***************************************************/

    /* Define the environment variables for printout */
    
    MES_IPGS = ipgs;
    MES_NGRF = ngrf[ipgs];
    
    for (MES_IGRF=0; MES_IGRF<2; MES_IGRF++)
    {
      if (! flag_used[ipgs][MES_IGRF]) continue;
      icase  = get_rank_from_propdef(propdef,ipgs,MES_IGRF);
      situba = st_alloc(models[ipgs][MES_IGRF],nbsimu,nbtuba);
      if (situba == (Situba *) NULL) goto label_end;
      if (st_simtub_process(dbin,dbout,models[ipgs][MES_IGRF],neigh,situba,
                            (double *) NULL,(double *) NULL,
                            nbsimu,icase,1,0,flag_check)) goto label_end;
      situba = st_dealloc(models[ipgs][MES_IGRF],situba);
    }
    
    /* Convert gaussian to facies at target point */
    
    if (! flag_gaus)
      for (int isimu=0; isimu<nbsimu; isimu++)
        simu_func_categorical_transf(dbout,0,isimu,nbsimu);
    
    /* Update facies proportions at target points */
    
    if (flag_modif)
    {
      for (int isimu=0; isimu<nbsimu; isimu++)
        simu_func_categorical_update(dbout,0,isimu,nbsimu);
      simu_func_categorical_scale(dbout,0,nbsimu);
    }
    
    /* Check/show facies at data against facies at the closest grid node */
    
    if (flag_cond && ! flag_gaus && (flag_check || flag_show))
      st_check_facies_data2grid(propdef,dbin,dbout,
                                rules[ipgs],flag_used[ipgs],
                                flag_stat,flag_check,flag_show,ipgs,
                                nechin,nvar,nfac[ipgs],nbsimu);
  }

  /********************************/
  /* Free the temporary variables */
  /********************************/

  if (dbout != nullptr)
  {
    if (!st_keep(flag_gaus, flag_modif, RESULT, PROP) && iptr_RP)
      dbout->deleteFieldByLocator(LOC_P);
    else
      namconv.setNamesAndLocators(NULL, LOC_Z, -1, dbout, iptr_RP, "Props",
                                  nfactot, false);

    if (!st_keep(flag_gaus, flag_modif, RESULT, GAUS) && iptr_RN)
      dbout->deleteFieldByLocator(LOC_SIMU);
    else
      namconv.setNamesAndLocators(NULL, LOC_Z, -1, dbout, iptr_RN, "Gaus",
                                  ngrftot * nbsimu, false);

    if (!st_keep(flag_gaus, flag_modif, RESULT, FACIES) && iptr_RF)
      dbout->deleteFieldByLocator(LOC_FACIES);
    else
      namconv.setNamesAndLocators(NULL, LOC_Z, -1, dbout, iptr_RF, String(),
                                  npgs * nbsimu);
  }

  if (dbin != nullptr)
  {
    if (!st_keep(flag_gaus, flag_modif, DATA, GAUS) && iptr_DN)
      dbin->deleteFieldByLocator(LOC_SIMU);
    else
      namconv.setNamesAndLocators(NULL, LOC_Z, -1, dbin, iptr_DN, "Gaus",
                                  ngrftot * nbsimu, false);

    if (!st_keep(flag_gaus, flag_modif, DATA, FACIES) && iptr_DF)
      dbin->deleteFieldByLocator(LOC_FACIES);
    else
      namconv.setNamesAndLocators(NULL, LOC_Z, -1, dbin, iptr_DF, String(),
                                  npgs * nbsimu, false);

    dbin->deleteFieldByLocator(LOC_GAUSFAC);
    dbin->deleteFieldByLocator(LOC_L);
    dbin->deleteFieldByLocator(LOC_U);
  }

  /* Set the error return flag */

  error = 0;

label_end:
  st_suppress_added_samples(dbin,nechin);
  for (ipgs=0; ipgs<npgs; ipgs++)
    for (igrf=0; igrf<2; igrf++)
      situba = st_dealloc(models[ipgs][igrf],situba);
  propdef = proportion_manage(-1,1,flag_stat,ngrf[0],ngrf[1],nfac[0],nfac[1],
                              dbin,dbprop,propcst,propdef);
  return(error);
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
GEOSLIB_API int db_simulations_to_ce(Db    *db,
                                     ENUM_LOCS locatorType,
                                     int    nbsimu,
                                     int    nvar,
                                     int   *iptr_ce_arg,
                                     int   *iptr_cstd_arg)
{
  int    error,iptr_ce,iptr_cstd,iptr_nb,nech;
  double value,count,mean,var;

  // Initializations

  error = 1;
  iptr_ce = iptr_cstd = iptr_nb = -1;
  if (db == (Db *) NULL) goto label_end;
  nech = db->getSampleNumber();
  if (nbsimu <= 0 || nvar <= 0 || nech <= 0) return(1);

  // Allocate the new attributes:

  iptr_ce = db->addFields(nvar,0.);
  if (iptr_ce < 0) goto label_end;
  iptr_cstd = db->addFields(nvar,0.);
  if (iptr_cstd < 0) goto label_end;
  iptr_nb = db->addFields(nvar,0.);
  if (iptr_nb < 0) goto label_end;

  // Loop on the simulations

  for (int isimu=0; isimu<nbsimu; isimu++)
  {
    // Loop on the samples

    for (int iech=0; iech<nech; iech++)
    {
      if (! db->isActive(iech)) continue;

      // Loop on the variables

      for (int ivar=0; ivar<nvar; ivar++)
      {
        // Arguments 'simu' and 'nvar' are interchanged to keep correct order
        value = db->getSimvar(locatorType,iech,ivar,isimu,0,nvar,nbsimu);
        if (FFFF(value)) continue;
        db->updArray(iech,iptr_ce  +ivar,0,value);
        db->updArray(iech,iptr_cstd+ivar,0,value * value);
        db->updArray(iech,iptr_nb  +ivar,0,1.);
      }
    }
  }
  
  // Scale the conditional expectation and variance

  for (int iech=0; iech<nech; iech++)
  {
    if (! db->isActive(iech)) continue;

    for (int ivar=0; ivar<nvar; ivar++)
    {
      count = db->getArray(iech,iptr_nb+ivar);
      if (count <= 0)
      {
        db->setArray(iech,iptr_ce  +ivar,TEST);
        db->setArray(iech,iptr_cstd+ivar,TEST);
      }
      else
      {
        mean = db->getArray(iech,iptr_ce  +ivar) / count;
        db->setArray(iech,iptr_ce+ivar,mean);
        var  = db->getArray(iech,iptr_cstd+ivar) / count - mean * mean;
        var  = (var > 0.) ? sqrt(var) : 0.;
        db->setArray(iech,iptr_cstd+ivar,var);
      }
    }
  }
  
  // Set the error return code

  error = 0;
  
label_end:
  (void) db_attribute_del_mult(db,iptr_nb,nvar);
  iptr_nb = -1;
  if (error)
  {
    (void) db_attribute_del_mult(db,iptr_ce  ,nvar);
    (void) db_attribute_del_mult(db,iptr_cstd,nvar);
    *iptr_ce_arg   = -1;
    *iptr_cstd_arg = -1;
  }
  else
  {
    *iptr_ce_arg   = iptr_ce;
    *iptr_cstd_arg = iptr_cstd;
  }
  return(error);
}

/****************************************************************************/
/*!
**  Perform the Gibbs sampler
**
** \return  Error return code
**
** \param[in]  dbin        Db structure
** \param[in]  model       Model structure
** \param[in]  neigh       Neigh structure (optional)
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
** \param[in]  gibbs_eps   Relative convergence criterion
** \param[in]  flag_ce     1 if the conditional expectation
**                         should be returned instead of simulations
** \param[in]  flag_cstd   1 if the conditional standard deviation
**                         should be returned instead of simulations
** \param[in]  verbose     Verbose flag
** \param[in]  namconv     Naming convention
**
*****************************************************************************/
GEOSLIB_API int gibbs_sampler(Db     *dbin,
                              Model  *model,
                              Neigh  *neigh,
                              int     nbsimu,
                              int     seed,
                              int     gibbs_nburn,
                              int     gibbs_niter,
                              bool    flag_norm,
                              bool    flag_multi_mono,
                              bool    flag_propagation,
                              int     gibbs_optstats,
                              double  percent,
                              double  gibbs_eps,
                              bool    flag_ce,
                              bool    flag_cstd,
                              bool    verbose,
                              NamingConvention namconv)
{
  int      error,iptr,isimu,npgs,nvar,iptr_ce,iptr_cstd;
  PropDef *propdef;

  /* Initializations */

  error   = 1;
  npgs    = 1;
  iptr_ce = iptr_cstd = -1;
  propdef = (PropDef *) NULL;

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

  if (model == (Model *) NULL)
  {
    messerr("No Model is provided");
    goto label_end;
  }
  nvar = model->getVariableNumber();
  if (! flag_propagation)
  {
    if (model_stabilize(model,1,percent)) goto label_end;
  }
  if (flag_norm)
  {
    if (model_normalize(model,1)) goto label_end;
  }

  /*******************/
  /* Core allocation */
  /*******************/

  law_set_random_seed(seed);

  propdef = proportion_manage(1,0,1,1,0,nvar,0,dbin,NULL,VectorDouble(),propdef);
  if (propdef == (PropDef *) NULL) goto label_end;

  /**********************/
  /* Add the attributes */
  /**********************/

  if (db_locator_attribute_add(dbin,LOC_GAUSFAC,nbsimu*nvar,0,0.,&iptr))
    goto label_end;

  /*****************/
  /* Gibbs sampler */
  /*****************/
  
  {
    AGibbs *gibbs;
    if (flag_multi_mono)
      gibbs = GibbsFactory::createGibbs(dbin, model, neigh);
    else
    {
      std::vector<Model *> modvec;
      modvec.push_back(model);
      gibbs = GibbsFactory::createGibbs(dbin, modvec, 0., flag_propagation);
    }
    if (gibbs == nullptr) goto label_end;
    gibbs->setOptionStats(gibbs_optstats);

    /* Initialize the Gibbs calculations */

    gibbs->init(npgs, nvar, gibbs_nburn, gibbs_niter,
                0, false, true, 0.);

    /* Allocate the covariance matrix inverted */

    if (gibbs->covmatAlloc(verbose)) goto label_end;

    // Allocate the Gaussian vector

    VectorVectorDouble y = gibbs->allocY();

    /* Loop on the simulations */

    for (isimu = 0; isimu < nbsimu; isimu++)
    {
      message("Processing Simulation %d/%d\n", isimu + 1, nbsimu);
      int ipgs = 0;

      // Initialize the iterations

      if (gibbs->calculInitialize(y, isimu, ipgs, verbose)) goto label_end;
      if (verbose) gibbs->print(true,y,isimu,ipgs);

      /* Iterations of the Gibbs sampler */

      for (int iter = 0; iter < gibbs->getNiter(); iter++)
        gibbs->update(y, isimu, ipgs, iter);
      if (verbose) gibbs->print(false,y,isimu,ipgs);

      // Store the results

      gibbs->storeResult(y, isimu, ipgs);
    }
  }

  /* Convert the simulations */

  if (flag_ce || flag_cstd)
  {
    if (db_simulations_to_ce(dbin,LOC_GAUSFAC,nbsimu,nvar,
                             &iptr_ce,&iptr_cstd)) goto label_end;

    // We release the attributes dedicated to simulations on Dbout

    if (! flag_ce)  
    {
      (void) db_attribute_del_mult(dbin,iptr_ce  ,nvar);
      iptr_ce = -1;
    }
    if (! flag_cstd)
    {
      (void) db_attribute_del_mult(dbin,iptr_cstd,nvar);
      iptr_cstd = -1;
    } 
    dbin->deleteFieldByLocator(LOC_GAUSFAC);
 }

  /* Set the error return flag */

  error = 0;
  namconv.setNamesAndLocators(dbin,LOC_UNKNOWN,nvar,
                              dbin,iptr,String(),nbsimu);

label_end:
  propdef = proportion_manage(-1,0,1,1,0,model->getVariableNumber(),0,dbin,NULL,
                              VectorDouble(),propdef);
  return(error);
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
** \param[in]  neigh        Neigh structure (optional)
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
GEOSLIB_API int simtub_constraints(Db       *dbin,
                                   Db       *dbout,
                                   Model    *model,
                                   Neigh    *neigh,
                                   int       seed,
                                   int       nbtuba,
                                   int       nbsimu_min,
                                   int       nbsimu_quant,
                                   int       niter_max,
                                   VectorInt& cols,
                                   int     (*func_valid)(int     flag_grid,
                                                         int     nDim,
                                                         int     nech,
                                                         int    *nx,
                                                         double *dx,
                                                         double *x0,
                                                         double  nonval,
                                                         double  percent,
                                                         VectorDouble& tab))
{
  int    *nx,iatt,retval,nbtest;
  int     error,nbsimu,nvalid,isimu,ndim,iter,nech,flag_grid,i;
  double *dx,*x0,percent;
  VectorDouble tab;

  /* Initializations */

  error  = 1;
  law_set_random_seed(seed);
  nx = (int *) NULL;
  dx = x0 = (double *) NULL;
  cols.clear();

  /* Preliminary check */

  flag_grid = is_grid(dbout);
  ndim = dbout->getNDim();
  nech = dbout->getSampleNumber();
  tab.resize(dbout->getSampleNumber());
  if (flag_grid)
  {
    nx   = (int    *) mem_alloc(sizeof(int)    * ndim,0);
    if (nx  == (int    *) NULL) goto label_end;
    dx   = (double *) mem_alloc(sizeof(double) * ndim,0);
    if (dx  == (double *) NULL) goto label_end;
    x0   = (double *) mem_alloc(sizeof(double) * ndim,0);
    if (x0  == (double *) NULL) goto label_end;

    for (i=0; i<ndim; i++)
    {
      nx[i] = dbout->getNX(i);
      dx[i] = dbout->getDX(i);
      x0[i] = dbout->getX0(i);
    }
  }

  /* Implicit loop on the simulations */

  iatt   = dbout->getFieldNumber();
  nvalid = iter = nbtest = 0;
  nbsimu = nbsimu_min + nbsimu_quant;
  while (nvalid < nbsimu_min && iter < niter_max)
  {

    /* Performing the simulations */

    iter++;
    nbtest += nbsimu;
    if (simtub(dbin,dbout,model,neigh,nbsimu,0,nbtuba,0)) goto label_end;

    /* Check if the simulated outcomes are valid */

    for (isimu=0; isimu<nbsimu; isimu++,iatt++)
    {

      /* Load the target simulation into the interface buffer */

      if (db_vector_get_att_sel(dbout,iatt,tab.data())) goto label_end;

      /* Check if the simulation is valid */

      percent = 100. * nvalid / nbsimu_min;
      retval  = func_valid(flag_grid,ndim,nech,nx,dx,x0,TEST,percent,tab);
      if (retval == 0)
      {

        /* Delete the current simulation */

        dbout->deleteFieldByAttribute(iatt);

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

          dbout->setFieldByAttribute(tab,iatt);
        }
        cols.push_back(iatt);
        nvalid++;
      }
    }

    /* Optional printout */

    if (debug_query("converge"))
      message("Iteration #%2d - Simulations %3d tested, %2d valid\n",
              iter,nbtest,nvalid);

    /* Define the number of simulations for the next batch */

    nbsimu = nbsimu_quant;
  }

  /* Set the error return code */

  error = 0;

label_end:
  nx  = (int    *) mem_free((char *) nx);
  dx  = (double *) mem_free((char *) dx);
  x0  = (double *) mem_free((char *) x0);
  return(error);
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
static int st_maxstable_mask(Db     *dbout,
                             double  seuil,
                             double  scale,
                             int     iptrv,
                             int     iptrs)
{
  int    iech,number;
  double valsim;

  for (iech=number=0; iech<dbout->getSampleNumber(); iech++)
  {
    if (! dbout->isActive(iech)) continue;
    valsim = dbout->getArray(iech,iptrv);
    if (valsim > seuil / scale) 
      dbout->setArray(iech,iptrs,0.);
    else
      number++;
  }
  return(number);
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
static void st_maxstable_combine(Db     *dbout,
                                 double  scale,
                                 int     iter0,
                                 int     iptrg,
                                 int     iptrv,
                                 int     iptrr,
                                 int    *last)
{
  int    iech;
  double valsim,valold;

  for (iech=0; iech<dbout->getSampleNumber(); iech++)
  {
    if (! dbout->isActive(iech)) continue;
    valold = dbout->getArray(iech,iptrv);
    valsim = dbout->getArray(iech,iptrg) / scale;
    if (valsim > valold)
    {
      dbout->setArray(iech,iptrv,valsim);
      dbout->setArray(iech,iptrr,iter0);
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
GEOSLIB_API int simmaxstable(Db    *dbout,
                             Model *model,
                             double ratio,
                             int    seed,
                             int    nbtuba,
                             int    flag_simu,
                             int    flag_rank,
                             int    verbose)
{
  Situba *situba;
  double  tpois,seuil;
  int     error,iptrg,iptrv,iptrr,iptrs,niter,nleft,icov,last;
  static double seuil_ref = 5.;

  /* Initializations */

  error  = 1;
  iptrv  = iptrg = iptrs = iptrr = -1;
  situba = (Situba *) NULL;
  law_set_random_seed(seed);
  if (st_check_simtub_environment(NULL,dbout,model,NULL)) goto label_end;
  seuil = get_keypone("MaxStableThresh",seuil_ref);

  /* Preliminary checks */

  if (model->getVariableNumber() != 1)
  {
    messerr("This feature is limited to the monovariate case");
    goto label_end;
  }
  if (! flag_simu && ! flag_rank)
  {
    messerr("You must choose 'flag_simu' or  'flag_rank' or both");
    goto label_end;
  }

  /* Define the environment variables for printout */

  st_simulation_environment();
  MES_NPGS = MES_NGRF = 1;

  /* Add the attributes for storing the results */

  iptrv = dbout->addFields(1,0.);
  if (iptrv < 0) goto label_end;
  iptrr = dbout->addFields(1,0.);
  if (iptrr < 0) goto label_end;
  if (db_locator_attribute_add(dbout,LOC_SEL,1,0,0.,&iptrs)) 
    goto label_end;
  if (db_locator_attribute_add(dbout,LOC_SIMU,1,0,0.,&iptrg)) 
    goto label_end;

  /* Implicit loop on the simulations */

  if (verbose)
  {
    message("Total number of cells = %d\n",dbout->getSampleNumber());
    message("Maximum simulation value = %lf\n",seuil);
  }

  tpois = 0.;
  niter = nleft = last = 0;
  while (1)
  {
    niter++;
    tpois -= log(law_uniform(0.,1.));

    /* Mask the nodes that cannot be accessed anymore */

    nleft = st_maxstable_mask(dbout,seuil,tpois,iptrv,iptrs);
    if (nleft <= 0) break;

    /* Processing the Turning Bands algorithm */
    
    situba = st_alloc(model,1,nbtuba);
    if (situba == (Situba *) NULL) goto label_end;
    if (st_simtub_process(NULL,dbout,model,NULL,situba,
                          (double *) NULL,(double *) NULL,
                          1,0,0,0,0)) goto label_end;
    situba = st_dealloc(model,situba);
    
    /* Combine the newly simulated outcome to the background */
    
    st_maxstable_combine(dbout,tpois,niter,iptrg,iptrv,iptrr,&last);

    if (verbose)
      message("Iteration #%2d - Scale = %6.3lf - Nb. cells left = %d\n",
              niter,tpois,nleft);

    /* Update the model for next iteration */

    for (icov=0; icov<model->getCovaNumber(); icov++)
      model->getCova(icov)->setRange(model->getCova(icov)->getRange() * ratio);
  }

  if (verbose)
  {
    message("Number of iterations = %d\n",niter);
    message("Rank of the last covering iteration = %d\n",last);
  }
  
  /* Set the error return flag */

  error = 0;

label_end:
  if (iptrs >= 0) dbout->deleteFieldByAttribute(iptrs);
  if (iptrg >= 0) dbout->deleteFieldByAttribute(iptrg);
  if (! flag_rank && iptrr >= 0) dbout->deleteFieldByAttribute(iptrr);
  if (! flag_simu && iptrv >= 0) dbout->deleteFieldByAttribute(iptrv);
  return(error);
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
static double st_quantile(Db     *dbout,
                          double  proba,
                          double *sort)
{
  int iech,nech,nval,rank;

  /* Initializations */

  nech = dbout->getSampleNumber();

  /* Load the non-masked simulated values */

  for (iech=nval=0; iech<nech; iech++)
  {
    if (! dbout->isActive(iech)) continue;
    sort[nval++] = dbout->getSimvar(LOC_SIMU,iech,0,0,0,1,1);
  }

  /* Sorting the array */

  ut_sort_double(0,nval,NULL,sort);

  /* Calculate the quantile */

  rank = (int) (proba * (double) nval);
  return(sort[rank]);
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
GEOSLIB_API int simRI(Db     *dbout,
                      Model  *model,
                      int     ncut,
                      double *zcut,
                      double *wcut,
                      int     seed,
                      int     nbtuba,
                      int     verbose)
{
  Situba *situba;
  double *pres,*pton,*sort,cumul,simval,proba,seuil;
  int     icut,error,iptrg,iptrs,nech,iech,count,total;

  /* Initializations */

  error  = 1;
  iptrg  = iptrs = -1;
  pres   = pton = sort = (double *) NULL;
  situba = (Situba *) NULL;
  nech   = dbout->getSampleNumber();
  law_set_random_seed(seed);
  if (st_check_simtub_environment(NULL,dbout,model,NULL)) goto label_end;

  /* Preliminary checks */

  if (model->getVariableNumber() != 1)
  {
    messerr("This feature is limited to the monovariate case");
    goto label_end;
  }

  /* Define the environment variables for printout */

  st_simulation_environment();
  MES_NPGS = MES_NGRF = 1;

  /* Add the attributes for storing the results */

  sort = (double *) mem_alloc(sizeof(double) * nech,0);
  if (sort == (double *) NULL) goto label_end;
  pton = (double *) mem_alloc(sizeof(double) * ncut,0);
  if (pton == (double *) NULL) goto label_end;
  pres = (double *) mem_alloc(sizeof(double) * (ncut-1),0);
  if (pres == (double *) NULL) goto label_end;
  if (db_locator_attribute_add(dbout,LOC_SEL,1,0,0.,&iptrs)) 
    goto label_end;
  if (db_locator_attribute_add(dbout,LOC_SIMU,1,0,0.,&iptrg)) 
    goto label_end;

  /* Preliminary calculations */

  cumul = 0.;
  for (icut=0; icut<ncut; icut++)
  {
    if (icut > 0 && zcut[icut] <= zcut[icut-1])
    {
      messerr("The cutoff values must be ordered increasingly");
      goto label_end;
    }
    if (wcut[icut] < 0) 
    {
      messerr("The weight of class (%d) cannot be negative",icut+1);
      goto label_end;
    }
    cumul += wcut[icut];
  }
  if (cumul <= 0.)
  {
    messerr("The sum of weights cannot be negative or null");
    goto label_end;
  }
  for (icut=0; icut<ncut; icut++) wcut[icut] /= cumul;
  pton[0] = 1.;
  for (icut=1; icut<ncut; icut++)   pton[icut] = pton[icut-1] - wcut[icut];
  for (icut=0; icut<ncut-1; icut++) pres[icut] = pton[icut+1] / pton[icut];

  /* Set the mask to the whole set of grid nodes */

  for (iech=0; iech<nech; iech++) dbout->setSelection(iech,1);

  /* Loop on the cutoff classes */

  total = 0;
  for (icut=0; icut<ncut; icut++)
  {

    /* Simulation in the non-masked part of the grid */
    
    situba = st_alloc(model,1,nbtuba);
    if (situba == (Situba *) NULL) goto label_end;
    if (st_simtub_process(NULL,dbout,model,NULL,situba,
                          (double *) NULL,(double *) NULL,
                          1,0,0,0,0)) goto label_end;
    situba = st_dealloc(model,situba);

    /* Look for the quantile */

    proba = 1. - pres[icut];
    seuil = (icut < ncut-1) ? st_quantile(dbout,proba,sort) : TEST;

    /* Update the current selection */

    for (iech=count=0; iech<nech; iech++)
    {
      if (! dbout->getSelection(iech)) continue;
      simval = dbout->getSimvar(LOC_SIMU,iech,0,0,0,1,1);
      if (! FFFF(seuil) && simval >= seuil) continue;
      dbout->setSimvar(LOC_SIMU,iech,0,0,0,1,1,(double) (icut+1));
      dbout->setSelection(iech,0);
      count++;
    }
    total += count;
    if (verbose)
      message("Level %3d - Proba=%lf - Affected=%7d - Total=%7d\n",
              icut+1,proba,count,total);
  }
  
  /* Perform the final coding */
  
  for (iech=0; iech<nech; iech++)
  {
    icut = (int) dbout->getSimvar(LOC_SIMU,iech,0,0,0,1,1);
    if (icut < 1 || icut > ncut) 
      dbout->setSimvar(LOC_SIMU,iech,0,0,0,1,1,TEST);
    else
      dbout->setSimvar(LOC_SIMU,iech,0,0,0,1,1,zcut[icut-1]);
  }
  
  /* Set the error return flag */

  error = 0;

label_end:
  sort = (double *) mem_free((char *) sort);
  pton = (double *) mem_free((char *) pton);
  pres = (double *) mem_free((char *) pres);
  if (iptrs >= 0) dbout->deleteFieldByAttribute(iptrs);
  return(error);
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
GEOSLIB_API int simpgs_spde(Db       *dbin,
                            Db       *dbout,
                            RuleProp *ruleprop,
                            Model    *model1,
                            Model    *model2,
                            const String& triswitch,
                            const VectorDouble& gext,
                            int     flag_gaus,
                            int     flag_modif,
                            int     flag_check,
                            int     flag_show,
                            int     nfacies,
                            int     seed,
                            int     nbsimu,
                            int     gibbs_nburn,
                            int     gibbs_niter,
                            int     ngibbs_int,
                            int     verbose,
                            double  percent)
{
  int     iptr,ngrf,igrf,nechin,error,nvar,flag_used[2],flag_cond;
  int     iptr_RF,iptr_RP;
  Model  *models[2];
  PropDef  *propdef;
  SPDE_Option s_option;

  /* Initializations */

  error     = 1;
  nvar      = 1;
  nechin    = 0;
  ngrf      = 0;
  propdef   = (PropDef *) NULL;
  models[0] = model1;
  models[1] = model2;
  iptr_RF   = iptr_RP = 0;
  iptr      = -1;
  flag_cond = (dbin != (Db *) NULL);
  law_set_random_seed(seed);

  if (ruleprop == nullptr)
  {
    messerr("RuleProp must be defined");
    return 1;
  }
  int flag_stat = ruleprop->isFlagStat();
  const Rule* rule = ruleprop->getRule();
  const VectorDouble& propcst = ruleprop->getPropCst();
  const Db* dbprop = ruleprop->getDbprop();

  if (rule->getModeRule() == RULE_SHADOW)
  {
    messerr("The 'Shadow' rule is not authorized");
    goto label_end;
  }
  if (rule->particularities(dbout,dbprop,model1,1,flag_stat)) goto label_end;

  /**********************/
  /* Preliminary checks */
  /**********************/

  /* Input Db */
  if (flag_cond)
  {
    nechin = dbin->getSampleNumber();
    if (! dbin->isVariableNumberComparedTo(1)) goto label_end;
  }

  /* Output Db */
  if (dbout == (Db *) NULL)
  {
    messerr("'dbout' is compulsory");
    goto label_end;
  }
  if (flag_modif && flag_gaus)
  {
    messerr("Calculating the facies proportions is incompatible with storing the Gaussian values");
    goto label_end;
  }

  /* Model */
  for (igrf=0; igrf<2; igrf++)
  {
    flag_used[igrf] = rule->isYUsed(igrf);
    if (! flag_used[igrf]) continue;
    ngrf++;
    if (models[igrf] == (Model *) NULL)
    {
      messerr("The Underlying GRF #%d is needed",igrf+1);
      messerr("No corresponding Model is provided");
      goto label_end;
    }
    if (models[igrf]->getVariableNumber() != 1)
    {
      messerr("The number of variables in the model #%d (%d) should be 1",
              igrf+1,model1->getVariableNumber());
      goto label_end;
    }
    if (model_stabilize(models[igrf],1,percent)) goto label_end;
    if (model_normalize(models[igrf],1)) goto label_end;
  }
  if (spde_check(dbin,dbout,model1,model2,verbose,gext,
                 1,1,1,0,0,1,flag_modif)) goto label_end;
  s_option = spde_option_alloc();
  spde_option_update(s_option,triswitch);

  /* Define the environment variables for printout */

  MES_NPGS = 1;
  MES_IPGS = 0;
  MES_NGRF = ngrf;

  /**********************/
  /* Add the attributes */
  /**********************/

  /* Storage of the facies proportions */
  if (flag_modif)
  {
    if (db_locator_attribute_add(dbout,LOC_P,nfacies,0,0.,&iptr_RP)) 
      goto label_end;
  }

  /* Storage of the simulations in the Output Db */
  if (db_locator_attribute_add(dbout,LOC_SIMU,nbsimu*ngrf,0,0.,&iptr_RF)) 
    goto label_end;

  if (flag_cond)
  {
    /* Lower bound at input data points */
    if (db_locator_attribute_add(dbin,LOC_L,ngrf,0,0.,&iptr)) 
      goto label_end;
    
    /* Upper bound at input data points */
    if (db_locator_attribute_add(dbin,LOC_U,ngrf,0,0.,&iptr)) 
      goto label_end;
  }

  /* Storage of the facies simulations in the Output Db */
  if (db_locator_attribute_add(dbout,LOC_FACIES,nbsimu,0,0.,&iptr_RF)) 
    goto label_end;

  propdef = proportion_manage(1,1,flag_stat,ngrf,0,nfacies,0,dbin,dbprop,
                              propcst,propdef);
  if (propdef == (PropDef *) NULL) goto label_end;
  if (! flag_gaus) simu_define_func_transf(simu_func_categorical_transf);
  simu_define_func_update(simu_func_categorical_update);
  simu_define_func_scale (simu_func_categorical_scale);
  ModCat.propdef = propdef;
  ModCat.rule    = rule;
  ModCat.ipgs    = 0;
  ModCat.flag_used[0] = flag_used[0];
  ModCat.flag_used[1] = flag_used[1];

  /****************************************/
  /* Convert facies into gaussian at data */
  /****************************************/

  proportion_rule_process(propdef,PROCESS_COPY);

  /* Initialize the Gibbs calculations */

  st_init_gibbs_params(rule->getRho());
    
  if (flag_cond)
    for (igrf=0; igrf<2; igrf++)
    {
      if (! flag_used[igrf]) continue;
      for (int isimu=0; isimu<nbsimu; isimu++)
      {
        if (rule->evaluateBounds(propdef,dbin,dbout,isimu,
                                 igrf,0,nbsimu)) goto label_end;
      }
    }
  
  /***************************************************/
  /* Perform the conditional simulation for each GRF */
  /***************************************************/

  if (spde_prepar(dbin,dbout,gext,s_option)) goto label_end;
  if (spde_process(dbin,dbout,s_option,nbsimu,gibbs_nburn,gibbs_niter,
                   ngibbs_int)) goto label_end;

  /* Check/show facies at data against facies at the closest grid node */
  
  if (flag_cond && ! flag_gaus && (flag_check || flag_show))
    st_check_facies_data2grid(propdef,dbin,dbout,rule,flag_used,
                              flag_stat,flag_check,flag_show,0,
                              nechin,nvar,nfacies,nbsimu);
  
  /********************************/
  /* Free the temporary variables */
  /********************************/

  if (! st_keep(flag_gaus,flag_modif,RESULT,PROP) && iptr_RP)
    dbout->deleteFieldByLocator(LOC_P);

  if (! st_keep(flag_gaus,flag_modif,RESULT,FACIES) && iptr_RF)
    dbout->deleteFieldByLocator(LOC_FACIES);

  if (! st_keep(flag_gaus,flag_modif,RESULT,GAUS))
    dbout->deleteFieldByLocator(LOC_SIMU);

  dbin->deleteFieldByLocator(LOC_L);
  dbin->deleteFieldByLocator(LOC_U);

  /* Set the error return flag */

  error = 0;

label_end:
  propdef = proportion_manage(-1,1,flag_stat,ngrf,0,nfacies,0,
                              dbin,dbprop,propcst,propdef);
  st_suppress_added_samples(dbin,nechin);
  return(error);
}

/****************************************************************************/
/*!
**  Check if the Model can be simulated using Turning Bands
**
** \return  1 if the Model is valid; 0 otherwise
**
** \param[in]  model    Model structure
**
*****************************************************************************/
GEOSLIB_API int simtub_workable(Model  *model)

{
  int workable,type;

  /* Initializations */

  workable = 1;

  /* Loop on the structures */
  
  for (int is=0; is<model->getCovaNumber(); is++)
  {
    type = model->getCovaType(is);
	    
    switch (type)
    {
      case COV_NUGGET:
      case COV_EXPONENTIAL:
      case COV_SPHERICAL:
      case COV_CUBIC:
      case COV_GAUSSIAN:
      case COV_SINCARD:
      case COV_BESSEL_J:
      case COV_BESSEL_K:
      case COV_STABLE:
      case COV_POWER:
      case COV_SPLINE_GC:
      case COV_LINEAR:
      case COV_ORDER1_GC:
      case COV_ORDER3_GC:
      case COV_ORDER5_GC:
        break;
	      
      default:
        workable = 0;
    }
  }
  return(workable);
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
** \param[in]  gibbs_eps   Relative immobile criterion
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
GEOSLIB_API int simcond(Db    *dbin,
                        Db    *dbout,
                        Model *model,
                        int    seed,
                        int    nbsimu,
                        int    nbtuba,
                        int    gibbs_nburn,
                        int    gibbs_niter,
                        double gibbs_eps,
                        int    flag_check,
                        int    flag_ce,
                        int    flag_cstd,
                        int    verbose)
{
  Neigh  *neigh;
  Situba *situba;
  PropDef  *propdef;
  int     nvar,error,iext,inostat,iptr,iptr_ce,iptr_cstd,ndim;

  /* Initializations */

  error   = 1;
  neigh   = (Neigh *) NULL;
  nvar    = model->getVariableNumber();
  ndim    = model->getDimensionNumber();
  iptr    = -1;
  situba  = (Situba *) NULL;
  propdef = (PropDef *) NULL;

  /* Preliminary checks */

  neigh = neigh_init_unique(ndim);
  law_set_random_seed(seed);
  if (st_check_simtub_environment(dbin,dbout,model,NULL)) goto label_end;
  if (manage_external_info(1,LOC_F,dbin,dbout,&iext)) goto label_end;
  if (manage_external_info(1,LOC_NOSTAT,dbin,dbout,
                           &inostat)) goto label_end;

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

  propdef = proportion_manage(1,0,1,1,0,nvar,0,dbin,NULL,VectorDouble(),propdef);
  if (propdef == (PropDef *) NULL) goto label_end;

  /* Define the environment variables for printout */

  st_simulation_environment();
  MES_NPGS = MES_NGRF = 1;

  /* Add the attributes for storing the results */

  if (db_locator_attribute_add(dbin,LOC_GAUSFAC,nbsimu,0,0.,
                               &iptr)) goto label_end;
  if (db_locator_attribute_add(dbin,LOC_SIMU,nvar*nbsimu,0,0.,
                               &iptr)) goto label_end;
  if (db_locator_attribute_add(dbout,LOC_SIMU,nvar*nbsimu,0,0.,
                               &iptr)) goto label_end;

  /*****************/
  /* Gibbs sampler */
  /*****************/

  {
    AGibbs* gibbs = GibbsFactory::createGibbs(dbin, model, nullptr);

    /* Initialize the Gibbs calculations */

    gibbs->init(1, 1, gibbs_nburn, gibbs_niter, 0, false, true, 0.);

    /* Allocate the covariance matrix inverted */

    if (gibbs->covmatAlloc(verbose)) goto label_end;

    // Allocate the Gaussian Vector

    VectorVectorDouble y = gibbs->allocY();

    /* Loop on the simulations */

    for (int isimu = 0; isimu < nbsimu; isimu++)
    {
      int ipgs = 0;

      /* Initialization for the Gibbs sampler */

      if (gibbs->calculInitialize(y,isimu, ipgs, verbose)) goto label_end;
      if (verbose) gibbs->print(true,y,isimu,ipgs);

      /* Iterations of the gibbs sampler */

      for (int iter = 0; iter < gibbs->getNiter(); iter++)
        gibbs->update(y,isimu, ipgs, iter);

      if (verbose) gibbs->print(false,y,isimu,ipgs);

      // Store the results

      gibbs->storeResult(y, isimu, ipgs);
    }
  }

  /* Processing the Turning Bands algorithm */

  situba = st_alloc(model,nbsimu,nbtuba);
  if (situba == (Situba *) NULL) goto label_end;
  if (st_simtub_process(dbin,dbout,model,neigh,situba,
                        (double *) NULL,(double *) NULL,
                        nbsimu,0,0,1,flag_check)) goto label_end;

  /* Free the temporary variables not used anymore */

  dbin->deleteFieldByLocator(LOC_GAUSFAC);
  dbin->deleteFieldByLocator(LOC_SIMU);

  /* Convert the simulations */

  if (flag_ce || flag_cstd)
  {
    if (db_simulations_to_ce(dbout,LOC_SIMU,nbsimu,nvar,
                             &iptr_ce,&iptr_cstd)) goto label_end;

    // We release the attributes dedicated to simulations on Dbout

    dbout->deleteFieldByLocator(LOC_SIMU);
    if (! flag_ce) 
    {
      (void) db_attribute_del_mult(dbout,iptr_ce  ,nvar);
      iptr_ce = -1;
    }
    if (! flag_cstd) 
    {
      (void) db_attribute_del_mult(dbout,iptr_cstd,nvar);
      iptr_cstd = -1;
    }
  }

  /* Set the error return flag */

  error = 0;

label_end:
  neigh = neigh_free(neigh);
  (void) manage_external_info(-1,LOC_F,dbin,dbout,&iext);
  (void) manage_external_info(-1,LOC_NOSTAT,dbin,dbout,&inostat);
  situba = st_dealloc(model,situba);
  return(error);
}
