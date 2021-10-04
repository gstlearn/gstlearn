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
#include "Gibbs/GibbsMulti.hpp"
#include "Gibbs/AGibbs.hpp"
#include "Model/Model.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Law.hpp"
#include "geoslib_old_f.h"
#include "geoslib_define.h"

#define MEAN(ivar,iech)          (mean[(ivar) * nech + (iech)])
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
** \param[in]  verbose       Verbose flag
**
*****************************************************************************/
int GibbsMulti::calculInitialize(VectorVectorDouble& y,
                             int isimu,
                             int ipgs,
                             bool verbose)
{
  Db* db = getDb();
  const Model* model = getModel();
  int nactive = db->getActiveSampleNumber();
  int nvar    = getNvar();

  /* Print the title */

  if (debug_query("converge"))
    mestitle(1,"Initial Values for Gibbs Sampler (Simu:%d - GS:%d)",
             isimu+1,ipgs+1);

  /* Loop on the variables */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int icase   = getRank(ipgs,ivar);

    /* Loop on the samples */

    double sk = sqrt(model->getTotalSill(ivar,ivar));
    for (int iact = 0; iact < nactive; iact++)
    {
      int iech = getSampleRank(iact);
      double vmin, vmax;
      if (_boundsCheck(iech, ipgs, ivar, &vmin, &vmax)) return 1;

      /* Compute the median value of the interval */

      double pmin = (FFFF(vmin)) ? 0. : law_cdf_gaussian(vmin);
      double pmax = (FFFF(vmax)) ? 1. : law_cdf_gaussian(vmax);
      y[icase][iact] = sk * law_invcdf_gaussian((pmin + pmax) / 2.);
    }
  }

  return(0);
}

