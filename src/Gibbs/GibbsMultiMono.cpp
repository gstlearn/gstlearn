/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "geoslib_old_f.h"
#include "geoslib_define.h"
#include "geoslib_enum.h"

#include "Gibbs/GibbsMultiMono.hpp"
#include "Gibbs/AGibbs.hpp"
#include "Model/Model.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Law.hpp"
#include "Basic/OptDbg.hpp"
#include "Db/Db.hpp"

#include <math.h>

GibbsMultiMono::GibbsMultiMono()
    : AGibbs(),
      _models(),
      _rho(0.)
{
}

GibbsMultiMono::GibbsMultiMono(Db* db, std::vector<Model *> models, double rho)
    : AGibbs(db),
      _models(models),
      _rho(rho)
{
}

GibbsMultiMono::GibbsMultiMono(const GibbsMultiMono &r)
  : AGibbs(r)
  , _models(r._models)
  , _rho(r._rho)
{
}

GibbsMultiMono& GibbsMultiMono::operator=(const GibbsMultiMono &r)
{
  if (this != &r)
  {
    AGibbs::operator=(r);
    _models = r._models;
    _rho = r._rho;
  }
  return *this;
}

GibbsMultiMono::~GibbsMultiMono()
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
int GibbsMultiMono::calculInitialize(VectorVectorDouble& y,
                                     int isimu,
                                     int ipgs)
{
  int nact = getSampleRankNumber();
  int nvar = getNvar();

  /* Print the title */

  if (OptDbg::query(EDbg::CONVERGE))
    mestitle(1,"Initial Values for Gibbs Sampler (Simu:%d - GS:%d)",
             isimu+1,ipgs+1);

  /* Loop on the variables */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int icase   = getRank(ipgs,ivar);

    /* Loop on the samples */

    double sk = sqrt(getModels(ivar)->getTotalSill(0,0));
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
double GibbsMultiMono::getSimulate(VectorVectorDouble& y,
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
  double vmin = db->getLocVariable(ELoc::L,iech, icase);
  double vmax = db->getLocVariable(ELoc::U,iech, icase);

  // Apply optional decay

  getBoundsDecay(iter, &vmin, &vmax);

  // In multi-mono case, correct from the previously (linked) variable

  double yval;
  double sval;
  if (ivar > 0)
  {
    int icase0 = getRank(ipgs, 0);
    double rho = getRho();
    double sqr = sqrt(1. - rho * rho);
    yval = yk * sqr + rho * y[icase0][iact];
    sval = sk * sqr;
  }
  else
  {
    yval = yk;
    sval = sk;
  }

  /* Update the definition interval */

  if (!FFFF(vmin)) vmin = (vmin - yval) / sval;
  if (!FFFF(vmax)) vmax = (vmax - yval) / sval;

  /* Draw an authorized normal value */

  if (FFFF(vmin) && FFFF(vmax))
    return (yk + sk * law_gaussian());
  else
    return (yk + sk * law_gaussian_between_bounds(vmin, vmax));
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
int GibbsMultiMono::checkGibbs(const VectorVectorDouble& y, int isimu, int ipgs)
{
  Db* db   = getDb();
  int nact = getSampleRankNumber();
  int nvar = getNvar();
  mestitle(1,"Checking gaussian values from Gibbs vs. bounds (PGS=%d Simu=%d)",
           ipgs+1,isimu+1);

  int nerror = 0;
  double sqr = sqrt(1. - _rho * _rho);

  /* Loop on the variables */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int icase   = getRank(ipgs,ivar);
    int icase0  = getRank(ipgs,0);

    /* Loop on the data */

    for (int iact=0; iact<nact; iact++)
    {
      int iech = getSampleRank(iact);
      double vmin = db->getLocVariable(ELoc::L,iech,icase);
      double vmax = db->getLocVariable(ELoc::U,iech,icase);
      if (FFFF(vmin)) vmin = -1.e30;
      if (FFFF(vmax)) vmax =  1.e30;

      /* Read the gaussian value */

      double gaus = y[icase][iact];
      if (ivar > 0)
        gaus = sqr * gaus + _rho * y[icase0][iact];

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
         message(STRING_NA);
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

