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
#include "Gibbs/GibbsStandard.hpp"
#include "Model/Model.hpp"
#include "Db/Db.hpp"
#include "Basic/Law.hpp"
#include "Morpho/Morpho.hpp"
#include "geoslib_f.h"

#define COVMAT(i,j)              (_covmat[(i) * nactive + (j)])

GibbsStandard::GibbsStandard()
  : AGibbs()
  , _covmat()
{
}

GibbsStandard::GibbsStandard(const GibbsStandard &r)
  : AGibbs(r)
  , _covmat(r._covmat)
{
}

GibbsStandard& GibbsStandard::operator=(const GibbsStandard &r)
{
  if (this != &r)
  {
    AGibbs::operator=(r);
    _covmat = r._covmat;
  }
  return *this;
}

GibbsStandard::~GibbsStandard()
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
int GibbsStandard::covmatAlloc(Db *db, Model *model, bool verbose)
{
  // Initialization

  int nactive = db->getActiveSampleNumber();

  // Core allocation

  _covmat.resize(nactive * nactive,0.);

  // Establish Covariance Matrix

  if (verbose) message("Establish Covariance matrix\n");
  if (model->isNoStat())
    model_covmat_nostat(model,db,db,-1,-1,0,1,_covmat.data());
  else
    model_covmat(model,db,db,-1,-1,0,1,_covmat.data());

  // Invert Covariance Matrix

  if (verbose) message("Invert Covariance matrix\n");
  if (matrix_invert(_covmat.data(),nactive,-1))
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
int GibbsStandard::calculInitialize(int flag_category,
                                    int flag_order,
                                    Db *db,
                                    Model *model,
                                    int isimu,
                                    int igrf,
                                    int ipgs,
                                    bool verbose)
{
  int    iech,nech,error,icov,icase,icase0;
  double simval,vmin,vmax,val1,ratio,sk,pmin,pmax;

  /* Core allocation */

  _checkMandatoryAttribute("Gibbs_Monovariate_Initialization",db,LOC_GAUSFAC);
  error  = 1;
  nech   = db->getSampleNumber();
  icase  = getRank(ipgs,igrf);
  icase0 = getRank(ipgs,0);

  /* Print the title */

  if (debug_query("converge"))
    mestitle(1,"Initial Values for Gibbs Sampler (GRF:%d - Simu:%d)",
             igrf+1,isimu+1);

  /* Get the theoretical standard deviation */

  sk = 0.;
  for (icov=0; icov<model->getCovaNumber(); icov++)
    sk += model->getSill(icov,0,0);
  sk = sqrt(sk);

  /* Loop on the samples */

  val1  = 0.;
  ratio = (igrf == 0) ? 1. : getSqr();
  for (iech=0; iech<nech; iech++)
  {
    if (! db->isActive(iech)) continue;
    if (_correctBoundsOrder(flag_category, flag_order, db, iech, 0, icase, 1,
                            &vmin, &vmax)) goto label_end;

    /* Draw a value as the median value of the interval */

    if (igrf == 1)
      val1 = getRho() * db->getSimvar(LOC_GAUSFAC,iech,isimu,0,icase0,getNbsimu(),1);
    if (! FFFF(vmin)) vmin = (vmin - val1) / ratio;
    if (! FFFF(vmax)) vmax = (vmax - val1) / ratio;

    /* Compute the median value of the interval */

    pmin   = (FFFF(vmin)) ? 0. : law_cdf_gaussian(vmin);
    pmax   = (FFFF(vmax)) ? 1. : law_cdf_gaussian(vmax);
    simval = sk * law_invcdf_gaussian(law_uniform(pmin,pmax));
    db->setSimvar(LOC_GAUSFAC,iech,isimu,0,icase,getNbsimu(),1,simval);
  }

  /* Optional printout */

  if (verbose)
    _gibbsInitPrint("Gibbs Sampler Initial",db,1,getNbsimu(),isimu,icase);

  /* Set the error return code */

  error = 0;

label_end:
  return(error);
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
** \remark Only available for monovariate case
** \remark Attributes LOC_GAUSFAC are mandatory
**
*****************************************************************************/
int GibbsStandard::calculIteration(Db *db,
                                   Model *model,
                                   int isimu,
                                   int ipgs,
                                   int igrf,
                                   int verbose)
{
  int     error,iech,iiech,jech,jjech,nech,nbdiv,nfois,nactive,iter,ncumul;
  int    *flag_h,icase,icase0,itest;
  double  vmin,vmax,delloc,old_mean,new_mean,refe,yk,sk,yval,ratio;
  double *y,*yhard;
  VectorDouble mean;

  /* Initializations */

  _checkMandatoryAttribute("Gibbs_Monovariate_Iteration",db,LOC_GAUSFAC);
  error   = 1;
  nech    = db->getSampleNumber();
  nactive = db->getActiveSampleNumber();
  y       = yhard = (double *) NULL;
  flag_h  = (int *) NULL;
  ratio   = (igrf == 0) ? 1. : getSqr();
  delloc  = new_mean = 0.;
  icase   = getRank(ipgs,igrf);
  icase0  = getRank(ipgs,0);

  /* Core allocation */

  mean.resize(nech,0.);
  y      = (double *) mem_alloc(sizeof(double) * nactive,0);
  if (y      == (double *) NULL) goto label_end;
  yhard  = (double *) mem_alloc(sizeof(double) * nactive,0);
  if (yhard  == (double *) NULL) goto label_end;
  flag_h = (int    *) mem_alloc(sizeof(double) * nech,0);
  if (flag_h == (int    *) NULL) goto label_end;

  /* Print the title */

  if (debug_query("converge"))
    mestitle(1,"Iterative Conditional Expectation (GRF:%d - Simu:%d)",
             igrf+1,isimu+1);

  /* Load the vector in memory */

  if (verbose) message("Starting Gibbs...\n");
  for (iech=iiech=0; iech<nech; iech++)
  {
    if (! db->isActive(iech)) continue;
    y[iiech] = db->getSimvar(LOC_GAUSFAC,iech,isimu,0,icase,getNbsimu(),1);
    iiech++;
  }

  /* Pre-calculations of kriging from hard to soft data */

  if (verbose) message("Counting the number of Hard and Soft data\n");
  for (iech=0; iech<nech; iech++)
  {
    flag_h[iech] = 0;
    if (! db->isActive(iech)) continue;
    vmin = db->getLowerBound(iech,icase);
    vmax = db->getUpperBound(iech,icase);
    flag_h[iech] = (vmin >= vmax) ?  1 : -1;
  }

  /* Perform the estimation at soft point from the hard information only */

  if (verbose)
    message("Contribution of Hard Data to the estimation of Soft Data\n");
  for (iech=iiech=0; iech<nech; iech++)
  {
    if (flag_h[iech] == 0) continue;
    if (flag_h[iech]  < 0)
    {
      yk = 0.;
      for (jech=jjech=0; jech<nech; jech++)
      {
        if (flag_h[jech] == 0) continue;
        if (flag_h[jech]  > 0) yk -= y[jjech] * COVMAT(iiech,jjech);
        jjech++;
      }
      yhard[iiech] = yk;
    }
    else
    {
      yhard[iiech] = TEST;
    }
    iiech++;
  }

  /* Loop on the iterations */

  if (verbose) message("Loop on the iterations and the samples\n");
  nfois = ncumul = itest = 0;
  for (iter=0; iter<getNburn()+getNiter(); iter++)
  {
    nfois++;

    /* Loop on the samples */

    for (iech=iiech=0; iech<nech; iech++,itest++)
    {
      mes_process("Gibbs processing",(getNburn()+getNiter())*nech,itest);
      if (flag_h[iech] == 0) continue;
      if (flag_h[iech]  < 0)
      {
        vmin = db->getLowerBound(iech,icase);
        vmax = db->getUpperBound(iech,icase);

        /* Perform the estimation from the other informations */

        sk = 1. / COVMAT(iiech,iiech);
        yk = yhard[iiech];
        for (jech=jjech=0; jech<nech; jech++)
        {
          if (flag_h[jech] == 0) continue;
          if (flag_h[jech]  < 0 && iiech != jjech)
            yk -= y[jjech] * COVMAT(iiech,jjech);
          jjech++;
        }
        yk *= sk;

        /* Correct from the already simulated Gaussian (Multigaussian case) */

        yval = yk;
        if (igrf > 0 && getRho() > 0.)
          yval = yk * getSqr() + getRho() *
            db->getSimvar(LOC_GAUSFAC,iech,isimu,0,icase0,getNbsimu(),1);

        /* Update the definition interval */

        sk = sqrt(sk);
        if (! FFFF(vmin)) vmin  = (vmin - yval) / (sk * ratio);
        if (! FFFF(vmax)) vmax  = (vmax - yval) / (sk * ratio);

        /* Draw an authorized normal value */

        if (_boundsCheck(db,iech,TEST,&vmin,&vmax,ITEST,TEST,TEST,TEST))
          messageAbort("Bounds for sample #%d are inverted during iteration #%d",
                    iech+1,iter+1);
        if (vmax > vmin)
          y[iiech] = yk + sk * law_gaussian_between_bounds(vmin,vmax);
        else
          y[iiech] = yk + sk * vmin;
      }
      iiech++;
    }

    /* Update the convergence criterion */

    if (iter+1 >= getNburn())
    {
      nbdiv = 0;
      ncumul++;
      for (iech=iiech=0; iech<nech; iech++)
      {
        if (! db->isActive(iech)) continue;
        old_mean     = (ncumul > 1) ? mean[iech] / (ncumul - 1) : 0.;
        mean[iech]  += y[iiech];
        new_mean     = mean[iech] / (ncumul);
        refe         = ABS(new_mean + old_mean) / 2.;
        delloc       = ABS(new_mean - old_mean);
        delloc       = (refe != 0.) ? 100. * delloc / refe : 0.;
        if (delloc > getEps()) nbdiv++;

        /* Optional printout */

        if (debug_query("converge"))
        {
          vmin = db->getLowerBound(iech,icase);
          vmax = db->getUpperBound(iech,icase);
          _printInequalities(db, iiech, iech, 0, nfois, ncumul > 1, y[iiech],
                             vmin, vmax, new_mean, delloc);
        }
        iiech++;
      }
      if (nbdiv <= 0)
      {
        message("Convergence reached after %d iterations\n",iter);
        goto label_store;
      }
    }
  }

  /* Storage */

label_store:
  for (iech=iiech=0; iech<nech; iech++)
  {
    if (! db->isActive(iech)) continue;
    db->setSimvar(LOC_GAUSFAC,iech,isimu,0,icase,getNbsimu(),1,y[iiech]);
    iiech++;
  }

  /* Scale the mean array */

  if (ncumul > 0)
    for (iech=0; iech<nech; iech++) mean[iech] /= (ncumul);

  /* Optional printout */

  if (verbose)
    _gibbsIterPrint("Gibbs Sampler Results", db, 1, isimu, iter, icase);

  /* Set the error return code */

  error = 0;

  /* Core deallocation */

label_end:
  y      = (double *) mem_free((char *) y);
  yhard  = (double *) mem_free((char *) yhard);
  flag_h = (int    *) mem_free((char *) flag_h);
  return(error);
}

/****************************************************************************/
/*!
**  Perform iteration of the Gibbs sampler (Propagative algorithm)
**
** \return  Error return code
**
** \param[in]  db          Db structure
** \param[in]  model       Model structure
** \param[in]  isimu       Rank of the simulation
** \param[in]  verbose     Verbose flag
**
** \remark The coefficient 'r' of the Gibbs Propagative algorithm
** \remark can be defined using:
** \remark set_keypair("gibbsPropaR",newval). Default 0.
**
** \remark The relative tolerance 'eps' of the Gibbs Propagative algorithm
** \remark can be defined using:
** \remark set_keypair("gibbsEps",newval). Default 0.
** \remark Attributes LOC_GAUSFAC are mandatory
**
*****************************************************************************/
int GibbsStandard::calculatePropagation(Db *db,
                                        Model *model,
                                        int isimu,
                                        bool verbose)
{
  int     iech,iiech,jech,jjech,iter,nech,nbdiv,nfois,ncumul,itest;
  int     ndim,idim,error,nactive,flag_affect,npart,icase;
  double  delloc,old_mean,new_mean,refe,eps;
  double *y,delta,sigval,sigloc,sqr,r;
  VectorUChar img;
  VectorDouble d1, mean;
  VectorInt nx;
  CovCalcMode mode;

  /* Initializations */

  _checkMandatoryAttribute("Gibbs_Monovariate_Propagation",db,LOC_GAUSFAC);
  error   = 1;
  nfois   = 0;
  nech    = db->getSampleNumber();
  nactive = db->getActiveSampleNumber();
  y       = (double *) NULL;
  delloc  = new_mean = 0.;
  r       = get_keypone("gibbsPropaR",0.);
  eps     = get_keypone("gibbsEps",0.);
  sqr     = sqrt(1. - r * r);
  ndim    = model->getDimensionNumber();
  icase   = getRank(0,0);

  /* Core allocation */

  nx.resize(2);
  nx[0] = nactive;
  nx[1] = nactive;
  mean.resize(nech,0.);
  y     = (double *) mem_alloc(sizeof(double) * nech,0);
  if (y     == (double *) NULL) goto label_end;
  d1.resize(ndim);
  img = morpho_image_manage(nx);

  /* Print the title */

  if (debug_query("converge"))
    mestitle(1,"Iterative Conditional Expectation (Simu:%d)",isimu+1);

  /* Load the vector in memory */

  for (iech=0; iech<nech; iech++)
  {
    y[iech] = (! db->isActive(iech)) ?
      TEST : db->getSimvar(LOC_GAUSFAC,iech,isimu,0,icase,getNbsimu(),1);
  }

  /* Loop on the iterations */

  nfois = ncumul = itest = npart = 0;
  for (iter=0; iter<getNburn()+getNiter(); iter++)
  {
    nfois++;

    /* Loop on the samples */

    for (iech=iiech=0; iech<nech; iech++,itest++)
    {
      mes_process("Propagation",(getNburn()+getNiter())*nech,itest);
      if (! db->isActive(iech)) continue;

      /* Covariance vector between the current datum and the other samples */

      for (idim=0; idim<ndim; idim++) d1[idim] = 0.;
      if (model->isNoStat())
      {
        CovInternal covint(1,iech,1,iech,ndim,db,db);
        model_calcul_cov_nostat(model,mode,&covint,1,1.,d1,&sigval);
      }
      else
        model_calcul_cov(model,mode,1,1.,d1,&sigval);
      if (sigval <= 0) continue;
      delta = (r - 1.) * y[iech] + sqrt(sigval) * sqr * law_gaussian();

      /* Update the gaussian vector */

      for (jech=jjech=0; jech<nech; jech++)
      {
        if (! db->isActive(jech)) continue;

        if (iter > 0 && ! bitmap_get_value(nx,img,iiech,jjech,0)) continue;

        for (idim=0; idim<ndim; idim++)
          d1[idim] = db->getCoordinate(iech,idim) - db->getCoordinate(jech,idim);
        if (model->isNoStat())
        {
          CovInternal covint(1,iech,1,jech,ndim,db,db);
          model_calcul_cov_nostat(model,mode,&covint,1,1.,d1,&sigloc);
        }
        else
          model_calcul_cov(model,mode,1,1.,d1,&sigloc);

        flag_affect = (ABS(sigloc) > sigval * eps);
        if (iter <= 0)
        {
          bitmap_set_value(nx,img,iiech,jjech,0,flag_affect);
          npart += flag_affect;
        }
        if (flag_affect) y[jech] += delta * sigloc / sigval;
        jjech++;
      }
      iiech++;
    }

    /* Update the convergence criterion */

    if (iter+1 >= getNburn())
    {
      nbdiv = 0;
      ncumul++;
      for (iech=iiech=0; iech<nech; iech++)
      {
        if (! db->isActive(iech)) continue;
        old_mean    = (ncumul > 1) ? mean[iech] / (ncumul - 1) : 0.;
        mean[iech] += y[iech];
        new_mean    = mean[iech] / (ncumul);
        refe        = ABS(old_mean + new_mean) / 2.;
        delloc      = ABS(new_mean - old_mean);
        delloc      = (refe != 0.) ? 100. * delloc / refe : 0.;
        if (delloc > getEps()) nbdiv++;

        /* Optional printout */

        if (debug_query("converge"))
          _printInequalities(db,iiech,iech,0,nfois,ncumul>1,y[iech],
                             TEST,TEST,new_mean,delloc);
        iiech++;
      }
      if (nbdiv <= 0) goto label_store;
    }
  }

  /* Storage */

label_store:
  for (iech=0; iech<nech; iech++)
    db->setSimvar(LOC_GAUSFAC,iech,isimu,0,icase,getNbsimu(),1,y[iech]);

  /* Scale the mean array */

  if (nfois > 0)
    for (iech=0; iech<nech; iech++) mean[iech] /= (ncumul);

  /* Optional printout */

  if (verbose)
    _gibbsIterPrint("Gibbs Sampler Results",db,1,isimu,iter,icase);

  /* Set the error return code */

  error = 0;

  /* Core deallocation */

label_end:
  y = (double *) mem_free((char *) y);
  return(error);
}

