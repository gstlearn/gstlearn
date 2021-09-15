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
#include "Gibbs/GibbsMoving.hpp"
#include "Model/Model.hpp"
#include "Db/Db.hpp"
#include "Basic/Law.hpp"
#include "Morpho/Morpho.hpp"
#include "geoslib_f.h"

#define MEAN(iech,ivar)          (mean[(ivar) * nech + (iech)])
#define COVMAT(i,j)              (covmat[(i) * nsize + (j)])

GibbsMoving::GibbsMoving()
  : AGibbs()
  , _wgt()
{
}

GibbsMoving::GibbsMoving(const GibbsMoving &r)
  : AGibbs(r)
  , _wgt(r._wgt)
{
}

GibbsMoving& GibbsMoving::operator=(const GibbsMoving &r)
{
  if (this != &r)
  {
    AGibbs::operator=(r);
    _wgt = r._wgt;
  }
  return *this;
}

GibbsMoving::~GibbsMoving()
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
** \param[in]  neigh       Neigh structure
** \param[in]  verbose     Verbose flag
**
*****************************************************************************/
int GibbsMoving::covmatAlloc(Db* db, Model* model, Neigh* neigh, bool verbose)
{
  VectorInt ivars, iechs;
  double *covmat;

  // Initialization

  int error = 1;
  int nvar = model->getVariableNumber();
  int nactive = db->getActiveSampleNumber();
  covmat = (double *) NULL;

  // Core allocation

  _wgt.resize(nactive);

  // Loop on the active samples

  int rank = 0;
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    GibbsWeights ww = _wgt[rank++];

    // Initialization
    ww._ivars = VectorInt();
    ww._iechs = VectorInt();
    ww._pivot = VectorInt();
    ww._sigma = VectorDouble();
    ww._ll    = VectorVectorDouble();

    // Select the Neighborhood
    getGeneralNeigh(db, neigh, iech, ww._ivars, ww._iechs);
    int nsize = ww._iechs.size();

    // Establishing the (moving) Covariance matrix
    covmat = model_covmat_by_variable_and_ranks(model, db, ww._ivars, ww._iechs, 0, 0);
    if (covmat == (double *) NULL) goto label_end;

    // Inverting the (moving) Covariance matrix
    if (matrix_invert(covmat, nsize, 0)) goto label_end;

    // Establishing the Vector of data
    ww._pivot.resize(nvar);
    ww._ll.resize(nvar);

    // Loop on the target variable
    for (int ivar = 0; ivar < nvar; ivar++)
    {

      // Look for the rank of the target point within the list
      int found = -1;
      for (int i = 0; i < nsize; i++)
      {
        int kech = ww._iechs[i];
        int kvar = ww._ivars[i];
        if (iech == kech && ivar == kvar) found = i;
      }
      ww._pivot[ivar] = found;

      int ecr = 0;
      ww._ll[ivar].resize(nsize);
      for (int i = 0; i < nsize; i++)
        ww._ll[ivar][ecr++] = COVMAT(i,found);
      ww._sigma[ivar] = 1. / COVMAT(found,found);
    }
  }

  error = 0;

  label_end: covmat = (double *) mem_free((char * ) covmat);
  return error;
}

/****************************************************************************/
/*!
**  Initializes the Gibbs sampler for a set of inequalities
**
** \return  Error return code
**
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
int GibbsMoving::calculInitialize(Db *db,
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
    if (_correctBoundsOrder(getFlagCategory(), getFlagOrder(),
                            db, iech, 0, icase, 1,
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
int GibbsMoving::calculIteration(Db *db,
                                   Model *model,
                                   int isimu,
                                   int ipgs,
                                   int igrf,
                                   int verbose)
{
  int error,nbdiv,nfois,ncumul,icase,iecr,pivot,rank;
  double  vmin,vmax,delloc,old_mean,new_mean,refe,yk,sk;
  double *y,*covmat;
  VectorDouble mean;

  /* Initializations */

  _checkMandatoryAttribute("Gibbs_Multivariate_Iteration",db,LOC_GAUSFAC);
  error   = 1;
  int nech    = db->getSampleNumber();
  int nvar    = model->getVariableNumber();
  int nactive = db->getActiveSampleNumber();
  int neq     = nvar * nactive;
  y       = covmat = (double *) NULL;
  delloc  = new_mean = 0.;
  icase   = getRank(ipgs,0);

  /* Core allocation */

  y = (double *) mem_alloc(sizeof(double) * neq, 0);
  if (y == (double *) NULL) goto label_end;
  mean.resize(nvar * nech, 0);

  /* Print the title */

  if (debug_query("converge"))
    mestitle(1, "Iterative Conditional Expectation (GS:%d - Simu:%d)", ipgs + 1,
             isimu + 1);

  /* Load the vector in memory */

  iecr = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int iech = 0; iech < nech; iech++)
    {
      if (!db->isActive(iech)) continue;
      y[iecr] = db->getSimvar(LOC_GAUSFAC, iech, 0, ivar, icase, 1, nvar);
      iecr++;
    }

  /* Loop on the iterations */

  nfois = ncumul = 0;
  for (int iter = 0; iter < getNburn() + getNiter(); iter++)
  {
    nfois++;

    /* Loop on the samples */

    rank = 0;
    for (int iech = 0; iech < nech; iech++)
    {
      if (!db->isActive(iech)) continue;
      GibbsWeights ww = _wgt[rank++];
      int nsize = ww._iechs.size();

      for (int ivar = 0; ivar < nvar; ivar++)
      {
        vmin = db->getLowerBound(iech, ivar);
        vmax = db->getUpperBound(iech, ivar);
        sk = ww._sigma[ivar];
        pivot = ww._pivot[ivar];

        /* Perform the estimation from the other informations */

        yk = 0.;
        for (int i2 = 0; i2 < nsize; i2++)
        {
          if (i2 != pivot) yk -= y[i2] * ww._ll[ivar][i2];
        }
        yk *= sk;

        /* Draw an authorized normal value */

        if (_boundsCheck(db, iech, TEST, &vmin, &vmax, ITEST, TEST, TEST, TEST))
        {
          messerr("Bounds for sample #%d and variable %d", iech + 1, ivar + 1);
          messerr("are inverted during iteration #%d", iter + 1);
          goto label_end;
        }
        y[iecr] = yk + sk * law_gaussian_between_bounds(vmin, vmax);
        iecr++;
      }

      /* Update the convergence criterion */

      if (iter + 1 >= getNburn())
      {
        ncumul++;
        nbdiv = iecr = 0;
        for (int ivar = 0; ivar < nvar; ivar++)
          for (int iech = 0; iech < nech; iech++)
          {
            if (!db->isActive(iech)) continue;
            old_mean = (ncumul > 1) ? MEAN(ivar,iech) / (ncumul - 1) :
                                      0.;
            MEAN(ivar,iech) += y[iecr];
            new_mean = MEAN(ivar,iech) / (ncumul);
            refe = ABS(new_mean + old_mean) / 2.;
            delloc = ABS(new_mean - old_mean);
            delloc = (refe != 0.) ? 100. * delloc / refe :
                                    0.;
            if (delloc > getEps()) nbdiv++;

            /* Optional printout */

            if (debug_query("converge"))
            {
              vmin = db->getLowerBound(iech, ivar);
              vmax = db->getUpperBound(iech, ivar);
              _printInequalities(db, iecr, iech, ivar, nfois, ncumul > 1,
                                 y[iecr], vmin, vmax, new_mean, delloc);
            }
            iecr++;
          }
        if (nbdiv <= 0) goto label_store;
      }
    }

    /* Optional printout */

    if (verbose)
      _gibbsIterPrint("Gibbs Sampler Results", db, nvar, isimu, iter, icase);
  }

  /* Storage */

  label_store:
  iecr = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int iech = 0; iech < nech; iech++)
    {
      if (!db->isActive(iech)) continue;
      db->setSimvar(LOC_GAUSFAC, iech, 0, ivar, icase, 1, nvar, y[iecr]);
      iecr++;
    }

  /* Scale the mean array */

  if (ncumul > 0)
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int iech = 0; iech < nech; iech++)
        MEAN(ivar,iech) /= (ncumul);

  /* Set the error return code */

  error = 0;

  /* Core deallocation */

  label_end: y = (double *) mem_free((char * ) y);
  covmat = (double *) mem_free((char * ) covmat);
  return (error);
}
