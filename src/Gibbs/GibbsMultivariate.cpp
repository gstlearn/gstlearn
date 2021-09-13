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
#include "Gibbs/GibbsMultivariate.hpp"
#include "Model/Model.hpp"
#include "Db/Db.hpp"
#include "Basic/Law.hpp"
#include "Morpho/Morpho.hpp"
#include "geoslib_f.h"

#define MEAN(iech,ivar)          (mean[(ivar) * nech + (iech)])
#define COVMAT(i,j)              (_covmat[(i) * nactive + (j)])

GibbsMultivariate::GibbsMultivariate()
  : AGibbs()
  , _covmat()
{
}

GibbsMultivariate::GibbsMultivariate(const GibbsMultivariate &r)
  : AGibbs(r)
  , _covmat(r._covmat)
{
}

GibbsMultivariate& GibbsMultivariate::operator=(const GibbsMultivariate &r)
{
  if (this != &r)
  {
    AGibbs::operator=(r);
    _covmat = r._covmat;
  }
  return *this;
}

GibbsMultivariate::~GibbsMultivariate()
{
}

/****************************************************************************/
/*!
**  Establish the covariance matrix for Gibbs
**
** \return  Error returned code
**
** \param[in]  db          Db structure
** \param[in]  model       Model structure
** \param[in]  verbose     Verbose flag
**
*****************************************************************************/
int GibbsMultivariate::covmatAlloc(Db *db, Model *model, bool verbose)
{
  // Initialization

  int nvar    = model->getVariableNumber();
  int nactive = db->getActiveSampleNumber();
  int neq     = nvar * nactive;

  // Core allocation

  _covmat.resize(neq * nactive,0.);

  // Establish Covariance Matrix

  if (verbose) message("Establish Covariance matrix\n");

  /* Establish the covariance matrix and invert it */

  model_covmat_multivar(model,db,0,1,_covmat.data());

  // Invert Covariance Matrix

  if (verbose) message("Invert Covariance matrix\n");
  if (matrix_invert(_covmat.data(),neq,-1))
  {
    messerr("Error during the covariance matrix inversion");
    return 1;
  }
  return 0;
}

/****************************************************************************/
/*!
**  Initializes the Gibbs sampler for a set of inequalities
**
** \return  Error return code
**
** \param[in]  flag_category 1 for categorical; 0 for continuous
** \param[in]  flag_order    Order relationship
** \li                        1 if the ascending order must be honored
** \li                       -1 if the descending order must be honored
** \li                        0 if no order relationship must be honored
** \param[in]  db            Db structure
** \param[in]  model         Model structure
** \param[in]  isimu         Rank of the simulation
** \param[in]  igrf          Rank of the GRF
** \param[in]  ipgs          Rank of the GS
** \param[in]  verbose       Verbose flag
**
** \remark Only available for monovariate case
** \remark This method is not documented on purpose. It should remain private
** \remark Attributes LOC_GAUSFAC are mandatory
**
*****************************************************************************/
int GibbsMultivariate::calculInitialize(int flag_category,
                                        int flag_order,
                                        Db *db,
                                        Model *model,
                                        int isimu,
                                        int igrf,
                                        int ipgs,
                                        bool verbose)
{
  int    iech,nech,icov,nvar,ivar,icase;
  double simval,vmin,vmax,sk,pmin,pmax;

  /* Core allocation */

  _checkMandatoryAttribute("Gibbs_Multivariate_Initialization",db,LOC_GAUSFAC);
  nech  = db->getSampleNumber();
  nvar  = model->getVariableNumber();
  icase = getRank(ipgs,0);

  /* Print the title */

  if (debug_query("converge"))
    mestitle(1,"Initial Values for Gibbs Sampler (GS:%d - Simu:%d)",
             ipgs+1,isimu+1);

  /* Loop on the samples */

  for (ivar=0; ivar<nvar; ivar++)
  {

    /* Get the theoretical standard deviation */

    sk = 0.;
    for (icov=0; icov<model->getCovaNumber(); icov++)
      sk += model->getSill(icov,ivar,ivar);
    sk = sqrt(sk);

    for (iech=0; iech<nech; iech++)
    {
      if (! db->isActive(iech)) continue;
      if (_correctBoundsOrder(flag_category, flag_order, db, iech, ivar, icase,
                              nvar, &vmin, &vmax)) return (1);

      /* Compute the median value of the interval */

      pmin   = (FFFF(vmin)) ? 0. : law_cdf_gaussian(vmin);
      pmax   = (FFFF(vmax)) ? 1. : law_cdf_gaussian(vmax);
      simval = sk * law_invcdf_gaussian((pmin + pmax) / 2.);
      db->setSimvar(LOC_GAUSFAC,iech,isimu,ivar,icase,getNbsimu(),nvar,simval);
    }
  }

  /* Optional printout */

  if (verbose)
    _gibbsInitPrint("Gibbs Sampler Initial",db,1,getNbsimu(),isimu,icase);

  return(0);
}

/****************************************************************************/
/*!
**  Perform iteration of the Gibbs sampler
**
** \return  Error return code
**
** \param[in]  db          Db structure
** \param[in]  model       Model structure
** \param[in]  isimu       Rank of the simulation
** \param[in]  ipgs        Rank of the GS
** \param[in]  igrf        Rank of the bounds (starting from 0)
** \param[in]  verbose     Verbose flag
**
** \param[out] mean        Working array for convergence criterion
**                         (Dimension: nech [even if masked samples])
**
** \remark Only available for monovariate case
** \remark Attributes LOC_GAUSFAC are mandatory
**
*****************************************************************************/
int GibbsMultivariate::calculIteration(Db *db,
                                       Model *model,
                                       int isimu,
                                       int ipgs,
                                       int igrf,
                                       int verbose,
                                       double *mean)
{
  int error,iech,jech,nech,nbdiv,nfois,nactive,iter,ncumul,neq,icase;
  int ivar,jvar,iecr,jecr,nvar;
  double  vmin,vmax,delloc,old_mean,new_mean,refe,yk,sk;
  double *y,*covmat;

  /* Initializations */

  _checkMandatoryAttribute("Gibbs_Multivariate_Iteration",db,LOC_GAUSFAC);
  error   = 1;
  nech    = db->getSampleNumber();
  nvar    = model->getVariableNumber();
  nactive = db->getActiveSampleNumber();
  neq     = nvar * nactive;
  y       = covmat = (double *) NULL;
  delloc  = new_mean = 0.;
  icase   = getRank(ipgs,0);

  /* Core allocation */

  y     = (double *) mem_alloc(sizeof(double) * neq,0);
  if (y     == (double *) NULL) goto label_end;

  /* Print the title */

  if (debug_query("converge"))
    mestitle(1,"Iterative Conditional Expectation (GS:%d - Simu:%d)",
             ipgs+1,isimu+1);

  /* Load the vector in memory */

  iecr = 0;
  for (ivar=0; ivar<nvar; ivar++)
    for (iech=0; iech<nech; iech++)
    {
      MEAN(ivar,iech) = 0.;
      if (! db->isActive(iech)) continue;
      y[iecr] = db->getSimvar(LOC_GAUSFAC,iech,0,ivar,icase,1,nvar);
      iecr++;
    }

  /* Loop on the iterations */

  nfois = ncumul = 0;
  for (iter=0; iter<getNburn()+getNiter(); iter++)
  {
    nfois++;

    /* Loop on the samples */

    iecr = 0;
    for (ivar=0; ivar<nvar; ivar++)
      for (iech=0; iech<nech; iech++)
      {
        if (! db->isActive(iech)) continue;
        vmin = db->getLowerBound(iech,icase);
        vmax = db->getUpperBound(iech,icase);
        sk = 1. / COVMAT(iecr,iecr);

        /* Perform the estimation from the other informations */

        yk = 0.;
        jecr = 0;
        for (jvar=0; jvar<nvar; jvar++)
          for (jech=0; jech<nech; jech++)
          {
            if (! db->isActive(jech)) continue;
            if (iecr != jecr) yk -= y[jecr] * COVMAT(iecr,jecr);
            jecr++;
          }
        yk *= sk;

        /* Draw an authorized normal value */

        if (_boundsCheck(db,iech,TEST,&vmin,&vmax,ITEST,TEST,TEST,TEST))
        {
          messerr("Bounds for sample #%d and variable %d",iech+1,ivar+1);
          messerr("are inverted during iteration #%d",iter+1);
          goto label_end;
        }
        y[iecr] = yk + sk * law_gaussian_between_bounds(vmin,vmax);
        iecr++;
      }

    /* Update the convergence criterion */

    if (iter+1 >= getNburn())
    {
      ncumul++;
      nbdiv = iecr = 0;
      for (ivar=0; ivar<nvar; ivar++)
        for (iech=0; iech<nech; iech++)
        {
          if (! db->isActive(iech)) continue;
          old_mean         = (ncumul > 1) ? MEAN(ivar,iech) / (ncumul - 1) : 0.;
          MEAN(ivar,iech) += y[iecr];
          new_mean         = MEAN(ivar,iech) / (ncumul);
          refe             = ABS(new_mean + old_mean) / 2.;
          delloc           = ABS(new_mean - old_mean);
          delloc           = (refe != 0.) ? 100. * delloc / refe : 0.;
          if (delloc > getEps()) nbdiv++;

          /* Optional printout */

          if (debug_query("converge"))
          {
            vmin = db->getLowerBound(iech,icase);
            vmax = db->getUpperBound(iech,icase);
            _printInequalities(db,iecr,iech,ivar,nfois,ncumul>1,y[iecr],
                               vmin,vmax,new_mean,delloc);
          }
          iecr++;
        }
      if (nbdiv <= 0) goto label_store;
    }
  }

  /* Storage */

label_store:
  iecr = 0;
  for (ivar=0; ivar<nvar; ivar++)
    for (iech=0; iech<nech; iech++)
    {
      if (! db->isActive(iech)) continue;
      db->setSimvar(LOC_GAUSFAC,iech,0,ivar,icase,1,nvar,y[iecr]);
      iecr++;
    }

  /* Scale the mean array */

  if (ncumul > 0)
    for (ivar=0; ivar<nvar; ivar++)
      for (iech=0; iech<nech; iech++)
        MEAN(ivar,iech) /= (ncumul);

  /* Optional printout */

  if (verbose)
    _gibbsIterPrint("Gibbs Sampler Results",db,nvar,isimu,iter,icase);

  /* Set the error return code */

  error = 0;

  /* Core deallocation */

label_end:
  y      = (double *) mem_free((char *) y);
  covmat = (double *) mem_free((char *) covmat);
  return(error);
}
