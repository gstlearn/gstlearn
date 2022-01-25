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

#include "Gibbs/GibbsMulti.hpp"
#include "Gibbs/AGibbs.hpp"
#include "Model/Model.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Law.hpp"
#include "Basic/OptDbg.hpp"
#include "Db/Db.hpp"
#include "geoslib_old_f.h"
#include "geoslib_define.h"

#include <math.h>

#define COVMAT(i,j)              (covmat[(i) * neq + (j)])

GibbsMulti::GibbsMulti()
    : AGibbs(),
      _model()
{
}

GibbsMulti::GibbsMulti(Db* db, Model * model)
    : AGibbs(db),
      _model(model)
{
}

GibbsMulti::GibbsMulti(const GibbsMulti &r)
  : AGibbs(r)
  , _model(r._model)
{
}

GibbsMulti& GibbsMulti::operator=(const GibbsMulti &r)
{
  if (this != &r)
  {
    AGibbs::operator=(r);
    _model = r._model;
  }
  return *this;
}

GibbsMulti::~GibbsMulti()
{
}

/****************************************************************************/
/*!
**  Initializes the Gibbs sampler for a set of inequalities
**
** \return  Error return code
**
** \param[in]  y             Gaussian vector
** \param[in]  isimu         Rank of the simulation
** \param[in]  ipgs          Rank of the GS
**
*****************************************************************************/
int GibbsMulti::calculInitialize(VectorVectorDouble& y,
                                 int isimu,
                                 int ipgs)
{
  const Model* model = getModel();
  int nact = getSampleRankNumber();
  int nvar = getNvar();

  /* Print the title */

  if (OptDbg::query(EDbg::CONVERGE))
    mestitle(1,"Initial Values for Gibbs Sampler (Simu:%d - GS:%d)",
             isimu+1,ipgs+1);

  /* Loop on the variables */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int icase = getRank(ipgs,ivar);

    /* Loop on the samples */

    double sk = sqrt(model->getTotalSill(ivar,ivar));
    for (int iact = 0; iact < nact; iact++)
    {
      double vmin, vmax;
      if (_boundsCheck(ipgs, ivar, iact, &vmin, &vmax)) return 1;

      /* Compute the median value of the interval */

      double pmin = (FFFF(vmin)) ? 0. : law_cdf_gaussian(vmin);
      double pmax = (FFFF(vmax)) ? 1. : law_cdf_gaussian(vmax);
      y[icase][iact] = sk * law_invcdf_gaussian((pmin + pmax) / 2.);
    }
  }
  return(0);
}

/**
 * Generate a simulated value
 * @param y     : Gaussian vector
 * @param yk    : Kriged value
 * @param sk    : Standard deviation
 * @param ipgs  : Rank of the current GS
 * @param ivar  : Rank of the current Variable
 * @param iact  : Rank of the target sample (relative)
 * @param iter  : Rank of the iteration
 * @return Simulated value
 */
double GibbsMulti::getSimulate(VectorVectorDouble& /*y*/,
                               double yk,
                               double sk,
                               int ipgs,
                               int ivar,
                               int iact,
                               int iter)
{
  // Define the environment

  int icase = getRank(ipgs, ivar);
  int iech  = getSampleRank(iact);

  // Read the Bounds

  const Db* db = getDb();
  double vmin = db->getLowerBound(iech, icase);
  double vmax = db->getUpperBound(iech, icase);

  // Apply optional decay

  getBoundsDecay(iter, &vmin, &vmax);

  /* Update the definition interval */

  if (!FFFF(vmin)) vmin = (vmin - yk) / sk;
  if (!FFFF(vmax)) vmax = (vmax - yk) / sk;

  /* Draw an authorized normal value */

  double value;
  if (FFFF(vmin) && FFFF(vmax))
    value = yk + sk * law_gaussian();
  else
    value = yk + sk * law_gaussian_between_bounds(vmin, vmax);
  return value;
}

/****************************************************************************/
/*!
**  Check/Show the facies against gaussian at wells
**
** \return Error return code
**
** \param[in]  y          Gaussian vector
** \param[in]  isimu      Rank of the simulation
** \param[in]  ipgs       Rank of the GS
**
*****************************************************************************/
int GibbsMulti::checkGibbs(const VectorVectorDouble& y, int isimu, int ipgs)
{
  Db* db = getDb();
  int nact = getSampleRankNumber();
  int nvar = getNvar();
  mestitle(1,"Checking gaussian values from Gibbs vs. bounds (PGS=%d Simu=%d)",
           ipgs+1,isimu+1);

  int nerror = 0;

  /* Loop on the variables */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int icase   = getRank(ipgs,ivar);

    /* Loop on the data */

    for (int iact=0; iact<nact; iact++)
    {
      int iech = getSampleRank(iact);
      double vmin = db->getLowerBound(iech,icase);
      double vmax = db->getUpperBound(iech,icase);
      if (FFFF(vmin)) vmin = -1.e30;
      if (FFFF(vmax)) vmax =  1.e30;

      /* Read the gaussian value */

      double gaus = y[icase][iact];

      /* Check inconsistency */

      if ((! FFFF(vmin) && gaus < vmin) ||
          (! FFFF(vmax) && gaus > vmax))
      {
        message("- Sample (#%d):",iech+1);
        message(" Simu#%d of Y%d=%lf",isimu+1,ivar+1,gaus);
        message(" does not lie within [");
        if (FFFF(vmin))
          message("NA,");
        else
          message("%lf",vmin);
        message(";");
        if (FFFF(vmax))
         message("NA");
        else
          message("%lf",vmax);
        message("]\n");
        nerror++;
      }
    }
  }
  if (nerror <= 0) message("No problem found\n");

  return nerror;
}

