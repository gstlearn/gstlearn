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
#include "../../include/Gibbs/GibbsUMultiMono.hpp"
#include "Model/Model.hpp"
#include "Db/Db.hpp"
#include "Basic/Law.hpp"
#include "Morpho/Morpho.hpp"
#include "geoslib_f.h"

#define COVMAT(ivar,i,j)              (_covmat[ivar][(i) * nactive + (j)])

GibbsUMultiMono::GibbsUMultiMono()
  : AGibbs()
  , _covmat()
{
}

GibbsUMultiMono::GibbsUMultiMono(Db* db, Model* model)
  : AGibbs(db, model)
  , _covmat()
{
}

GibbsUMultiMono::GibbsUMultiMono(const GibbsUMultiMono &r)
  : AGibbs(r)
  , _covmat(r._covmat)
{
}

GibbsUMultiMono& GibbsUMultiMono::operator=(const GibbsUMultiMono &r)
{
  if (this != &r)
  {
    AGibbs::operator=(r);
    _covmat = r._covmat;
  }
  return *this;
}

GibbsUMultiMono::~GibbsUMultiMono()
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
int GibbsUMultiMono::covmatAlloc(bool verbose)
{
  Db* db = getDb();
  Model* model = getModel();

  // Initialization

  if (verbose) mestitle(1,"Gibbs using Unique Neighborhood in MultiMono case");
  int nactive = db->getActiveSampleNumber();
  int nvar    = model->getVariableNumber();
  _covmat.resize(nvar);

  // Loop on the variables

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    _covmat[ivar].resize(nactive * nactive, 0.);

    // Establish Covariance Matrix

    if (verbose) message("Establish Covariance matrix (Var=%d)\n",ivar+1);
    model_covmat(model, db, db, ivar, ivar, 0, 1, _covmat[ivar].data());

    // Invert Covariance Matrix

    if (verbose) message("Invert Covariance matrix (Var=%d)\n",ivar+1);
    if (matrix_invert(_covmat[ivar].data(), nactive, -1))
    {
      messerr("Error during the covariance matrix inversion");
      return 1;
    }
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
void GibbsUMultiMono::update(VectorVectorDouble& y,
                             int isimu,
                             int ipgs,
                             int iter)
{
  Db* db = getDb();
  int nactive = db->getActiveSampleNumber();
  int nvar    = getNvar();

  /* Print the title */

  if (debug_query("converge"))
    mestitle(1,"Iterative Conditional Expectation (PGS=%d - Simu:%d - Iter=%d)",
             ipgs+1,isimu+1,iter+1);

  /* Loop on the variables */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int icase   = getRank(ipgs,ivar);

    /* Loop on the samples */

    for (int iact = 0; iact < nactive; iact++)
    {

      /* Perform the estimation from the other informations */

      double sk = 1. / COVMAT(ivar, iact, iact);
      double yk = 0.;
      for (int jact = 0; jact < nactive; jact++)
      {
        if (iact != jact) yk -= y[icase][jact] * COVMAT(ivar, iact, jact);
      }
      yk *= sk;

      /* Draw the simulated Gaussian */

      y[icase][iact] = getSimulate(y, yk, sk, iact, ipgs, ivar, iter);
    }
  }
}
