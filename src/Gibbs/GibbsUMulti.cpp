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
#include "../../include/Gibbs/GibbsUMulti.hpp"
#include "Model/Model.hpp"
#include "Db/Db.hpp"
#include "Basic/Law.hpp"
#include "Morpho/Morpho.hpp"
#include "geoslib_f.h"

#define MEAN(ivar,iech)          (mean[(ivar) * nech + (iech)])
#define COVMAT(i,j)              (_covmat[(i) * neq + (j)])

GibbsUMulti::GibbsUMulti()
  : AGibbs()
  , _covmat()
{
}

GibbsUMulti::GibbsUMulti(Db* db, Model* model)
  : AGibbs(db, model)
  , _covmat()
{
}

GibbsUMulti::GibbsUMulti(const GibbsUMulti &r)
  : AGibbs(r)
  , _covmat(r._covmat)
{
}

GibbsUMulti& GibbsUMulti::operator=(const GibbsUMulti &r)
{
  if (this != &r)
  {
    AGibbs::operator=(r);
    _covmat = r._covmat;
  }
  return *this;
}

GibbsUMulti::~GibbsUMulti()
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
int GibbsUMulti::covmatAlloc(bool verbose)
{
  // Initialization

  if (verbose) mestitle(1,"Gibbs using Unique Neighborhood");
  Db* db = getDb();
  Model* model = getModel();
  int nvar    = model->getVariableNumber();
  int nactive = db->getActiveSampleNumber();
  int neq     = nvar * nactive;

  // Core allocation

  _covmat.resize(neq * neq,0.);

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
**  Perform one update of the Gibbs sampler
**
** \param[in]  y           Gaussian vector
** \param[in]  isimu       Rank of the simulation
** \param[in]  ipgs        Rank of the GS
** \param[in]  iter        Rank of the iteration
**
*****************************************************************************/
void GibbsUMulti::update(VectorVectorDouble& y,
                         int isimu,
                         int ipgs,
                         int iter)
{
  Db* db = getDb();
  Model* model = getModel();
  int nvar    = getNvar();
  int nactive = db->getActiveSampleNumber();
  int neq     = nvar * nactive;

  /* Print the title */

  if (debug_query("converge"))
    mestitle(1,"Iterative Conditional Expectation (GS:%d - Simu:%d)",
             ipgs+1,isimu+1);

  /* Loop on the samples */

  for (int ivar = 0, iecr = 0; ivar < nvar; ivar++)
  {
    int icase   = getRank(ipgs,ivar);
    for (int iact = 0; iact < nactive; iact++, iecr++)
    {

      /* Perform the estimation from the other informations */

      double sk = 1. / COVMAT(iecr, iecr);
      double yk = 0.;
      for (int jvar = 0, jecr = 0; jvar < nvar; jvar++)
        for (int jact = 0; jact < nactive; jact++, jecr++)
        {
          if (iecr != jecr) yk -= y[icase][jact] * COVMAT(iecr, jecr);
        }
      yk *= sk;

      /* Draw an authorized normal value */

      y[icase][iact] = getSimulate(y, yk, sk, iact, ipgs, ivar, iter);
    }
  }
}
