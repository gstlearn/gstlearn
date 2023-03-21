/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "geoslib_old_f.h"

#include "Gibbs/GibbsUMulti.hpp"
#include "Model/Model.hpp"
#include "Db/Db.hpp"
#include "Basic/Law.hpp"
#include "Basic/OptDbg.hpp"
#include "Morpho/Morpho.hpp"

#include <math.h>

#define COVMAT(i,j)              (_covmat[(i) * neq + (j)])

GibbsUMulti::GibbsUMulti()
  : GibbsMulti()
  , _covmat()
{
}

GibbsUMulti::GibbsUMulti(Db* db, Model* model)
  : GibbsMulti(db, model)
  , _covmat()
{
}

GibbsUMulti::GibbsUMulti(const GibbsUMulti &r)
  : GibbsMulti(r)
  , _covmat(r._covmat)
{
}

GibbsUMulti& GibbsUMulti::operator=(const GibbsUMulti &r)
{
  if (this != &r)
  {
    GibbsMulti::operator=(r);
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
** \param[in]  verbose      Verbose flag
** \param[in]  verboseTimer True to show elapse times
**
*****************************************************************************/
int GibbsUMulti::covmatAlloc(bool verbose, bool /*verboseTimer*/)
{
  // Initialization

  if (verbose) mestitle(1,"Gibbs using Unique Neighborhood");
  Db* db = getDb();
  Model* model = getModel();
  int nvar = model->getVariableNumber();
  int nact = getSampleRankNumber();
  int neq  = nvar * nact;

  // Core allocation

  _covmat.resize(neq * neq,0.);

  // Establish Covariance Matrix

  if (verbose) message("Establish Covariance matrix\n");

  /* Establish the covariance matrix and invert it */

  model_covmat(model,db,db,-1,-1,0,1,_covmat.data());

  // Invert Covariance Matrix

  if (verbose) message("Invert Covariance matrix\n");
  if (matrix_invert(_covmat.data(),neq,-1))
  {
    messerr("Error during the covariance matrix inversion");
    return 1;
  }

  // Initialize the statistics (optional)

  statsInit();

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
  double valsim;
  int nvar = getNvar();
  int nact = getSampleRankNumber();
  int neq  = nvar * nact;

  /* Print the title */

  if (OptDbg::query(EDbg::CONVERGE))
    mestitle(1,"Iterative Conditional Expectation (GS:%d - Simu:%d)",
             ipgs+1,isimu+1);

  /* Loop on the target */

  for (int ivar = 0, iecr = 0; ivar < nvar; ivar++)
  {
    int icase = getRank(ipgs,ivar);
    for (int iact = 0; iact < nact; iact++, iecr++)
    {
      if (!isConstraintTight(ipgs, ivar, iact, &valsim))
      {

        /* Loop on the Data */

        double yk = 0.;
        double vark = 1. / COVMAT(iecr, iecr);
        for (int jvar = 0, jecr = 0; jvar < nvar; jvar++)
        {
          int jcase = getRank(ipgs, jvar);
          for (int jact = 0; jact < nact; jact++, jecr++)
          {
            if (ivar != jvar || iact != jact)
              yk -= y[jcase][jact] * COVMAT(iecr, jecr);
          }
        }
        yk *= vark;
        valsim = getSimulate(y, yk, sqrt(vark), ipgs, ivar, iact, iter);
      }
      y[icase][iact] = valsim;
    }
  }

  // Update statistics (optional)

  updateStats(y, ipgs, iter);
}
