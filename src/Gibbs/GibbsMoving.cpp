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

#define MEAN(ivar,iech)          (mean[(ivar) * nech + (iech)])
#define COVMAT(i,j)              (covmat[(i) * neq + (j)])

GibbsMoving::GibbsMoving()
  : AGibbs()
  , _neigh(nullptr)
  , _wgt()
{
}

GibbsMoving::GibbsMoving(Db* db, Model* model, Neigh* neigh)
  : AGibbs(db, model)
  , _neigh(neigh)
  , _wgt()
{
}

GibbsMoving::GibbsMoving(const GibbsMoving &r)
  : AGibbs(r)
  , _neigh(r._neigh)
  , _wgt(r._wgt)
{
}

GibbsMoving& GibbsMoving::operator=(const GibbsMoving &r)
{
  if (this != &r)
  {
    AGibbs::operator=(r);
    _neigh = r._neigh;
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
** \param[in]  verbose     Verbose flag
**
*****************************************************************************/
int GibbsMoving::covmatAlloc(bool verbose)
{
  VectorInt ivars, iechs;
  double *covmat;

  // Initialization

  if (verbose) mestitle(1,"Gibbs using Moving Neighborhood");
  int error = 1;
  Db* db = getDb();
  Model* model = getModel();
  Neigh* neigh = getNeigh();
  int nvar = model->getVariableNumber();
  int nech = db->getSampleNumber();
  int nactive = db->getActiveSampleNumber();
  int nvardb = db->getVariableNumber();
  bool flag_var_defined = nvardb > 0;
  covmat = (double *) NULL;
  if (defineGeneralNeigh(1, db, model, neigh)) return 1;

  // Consistency check

  if (flag_var_defined && nvar != nvardb)
  {
    messerr("Inconsistency in Number of Variables between Model (%d) and Db (%d)",
            nvar,nvardb);
    return 1;
  }

  // Clear the set of weight vectors

  if (verbose) message("Establish Weights\n");
  _wgt.resize(nactive);

  // Loop on the active samples

  int iiech = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    if (! db->isActive(iech)) continue;
    GibbsWeights& ww = _wgt[iiech++];

    // Select the Neighborhood
    ww._ll.clear();
    ww._ranks = getGeneralNeigh(db, neigh, iech);
    int nsize = ww._ranks.size();
    int neq   = nvar * nsize;
    ww._neq   = neq;

    // Establishing the (moving) Covariance matrix
    covmat = model_covmat_by_varranks(model, db, ww._ranks, neq, 0, 1);
    if (covmat == (double *) NULL) goto label_end;

    // Inverting the (moving) Covariance matrix
    if (matrix_invert(covmat, neq, 0)) goto label_end;

    // Store the weights per variable for the target sample

    for (int ivar = 0; ivar < nvar; ivar++)
    {

      // Find the rank of the pivot

      int found = -1;
      for (int jech = 0; jech < nsize; jech++)
        if (iech == ww._ranks[jech]) found = jech;
      ww._pivot = found;

      // Dimension the vector of weights

      VectorDouble ll(neq);
      for (int i = 0; i < ww._neq; i++)
        ll[i] = COVMAT(i, found);
      ww._ll.push_back(ll);
    }
    covmat = (double *) mem_free((char *) covmat);
  }
  error = 0;

  label_end:
  (void) defineGeneralNeigh(-1, db, model, neigh);
  covmat = (double *) mem_free((char * ) covmat);
  return error;
}

/****************************************************************************/
/*!
**  Perform one update of the Gibbs sampler
**
** \param[in]  y           Gaussian Vector
** \param[in]  isimu       Rank of the simulation
** \param[in]  ipgs        Rank of the GS
** \param[in]  iter        Rank of the iteration
**
*****************************************************************************/
void GibbsMoving::update(VectorVectorDouble& y,
                         int isimu,
                         int ipgs,
                         int iter)
{
  Db* db = getDb();
  Model* model = getModel();
  int nactive = db->getActiveSampleNumber();
  int nvar    = model->getVariableNumber();

  /* Print the title */

  if (debug_query("converge"))
    mestitle(1, "Gibbs Sampler (Simu:%d - GS:%d)", isimu + 1, ipgs + 1);

  /* Loop on the variables */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int icase   = getRank(ipgs,ivar);

    /* Loop on the samples */

    for (int iact = 0; iact < nactive; iact++)
    {
      const GibbsWeights& ww = _wgt[iact];
      int nsize = ww._ranks.size();
      int neq = ww._neq;
      int pivot = ww._pivot;

      for (int ivar = 0; ivar < nvar; ivar++)
      {
        double sk = 1. / ww._ll[ivar][pivot];
        int rank = pivot + nsize * ivar;
        double yk = 0.;
        for (int i2 = 0; i2 < neq; i2++)
        {
          if (i2 != rank) yk -= y[icase][i2] * ww._ll[ivar][i2];
        }
        yk *= sk;

        /* Draw an authorized normal value */

        y[icase][iact] = getSimulate(y, yk, sk, iact, ipgs, ivar, iter);
      }
    }
  }

  // Update statistics (optional)

  updateStats(y, ipgs, iter);
}
