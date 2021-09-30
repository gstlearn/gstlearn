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
#include "../../include/Gibbs/GibbsUPropMono.hpp"
#include "Model/Model.hpp"
#include "Db/Db.hpp"
#include "Basic/Law.hpp"
#include "Morpho/Morpho.hpp"
#include "geoslib_f.h"

GibbsUPropMono::GibbsUPropMono()
  : AGibbs()
  , _eps(EPSILON3)
{
}

GibbsUPropMono::GibbsUPropMono(Db* db, Model* model)
  : AGibbs(db, model)
  , _eps(EPSILON3)
{
}

GibbsUPropMono::GibbsUPropMono(const GibbsUPropMono &r)
  : AGibbs(r)
  , _eps(r._eps)
{
}

GibbsUPropMono& GibbsUPropMono::operator=(const GibbsUPropMono &r)
{
  if (this != &r)
  {
    AGibbs::operator=(r);
    _eps = r._eps;
  }
  return *this;
}

GibbsUPropMono::~GibbsUPropMono()
{
}

/****************************************************************************/
/*!
**  Establish the covariance matrix for Gibbs
**
** \return  Error returned code
**
** \param[in]  verbose     Verbose flag
**
*****************************************************************************/
int GibbsUPropMono::covmatAlloc(bool verbose)
{
  if (verbose) mestitle(1,"Gibbs using Unique Neighborhood in Propagative case");

  // Initialize the statistics (optional)

  statsInit();

  return 0;
}

/****************************************************************************/
/*!
**  Perform one update of the Gibbs sampler (Propagative algorithm)
**
** \param[in]  y           Gaussian vector
** \param[in]  isimu       Rank of the simulation
** \param[in]  ipgs        Rank of the GS (should be 0)
** \param[in]  iter        Rank of the iteration
**
*****************************************************************************/
void GibbsUPropMono::update(VectorVectorDouble& y,
                            int isimu,
                            int ipgs,
                            int iter)
{
  CovCalcMode mode;

  /* Initializations */

  Db* db = getDb();
  Model* model = getModel();
  int nactive = db->getActiveSampleNumber();
  int ndim    = model->getDimensionNumber();
  int icase   = getRank(ipgs,0);

  double sqr  = getSqr();
  double eps  = getEps();
  double r    = 0.1;

  /* Core allocation */

  VectorDouble d1(ndim);
  VectorBool img(nactive * nactive);

  /* Print the title */

  if (debug_query("converge"))
    mestitle(1,"Iterative Conditional Expectation (Simu:%d)",isimu+1);

  /* Loop on the samples */

  for (int iact = 0; iact < nactive; iact++)
  {
    int iech = getSampleRank(iact);

    /* Covariance vector between the current datum and the other samples */

    double sigval;
    for (int idim = 0; idim < ndim; idim++)
      d1[idim] = 0.;
    if (model->isNoStat())
    {
      CovInternal covint(1, iech, 1, iech, ndim, db, db);
      model_calcul_cov_nostat(model, mode, &covint, 1, 1., d1, &sigval);
    }
    else
      model_calcul_cov(model, mode, 1, 1., d1, &sigval);
    if (sigval <= 0) continue;
    sigval = sqrt(sigval);
    double delta = (r - 1.) * y[icase][iact] + sigval * sqr * law_gaussian();

    /* Update the gaussian vector */

    for (int jact = 0; jact < nactive; jact++)
    {
      if (iter > 0 && ! img[nactive * iact + jact]) continue;
      int jech = getSampleRank(jact);

      double sigloc;
      for (int idim = 0; idim < ndim; idim++)
        d1[idim] = db->getCoordinate(iech, idim) - db->getCoordinate(jech, idim);
      if (model->isNoStat())
      {
        CovInternal covint(1, iech, 1, jech, ndim, db, db);
        model_calcul_cov_nostat(model, mode, &covint, 1, 1., d1, &sigloc);
      }
      else
        model_calcul_cov(model, mode, 1, 1., d1, &sigloc);

      bool flag_affect = (ABS(sigloc) > sigval * eps);
      if (iter <= 0) img[nactive * iact + jact] = flag_affect;
      if (flag_affect) y[icase][jact] += delta * sigloc / sigval;
    }
  }

  // Update statistics (optional)

  updateStats(y, isimu, ipgs, iter);
}
