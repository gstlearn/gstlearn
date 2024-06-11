/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
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
  int nact = _getSampleRankNumber();
  int neq  = nvar * nact;

  // Establish Covariance Matrix

  if (verbose) message("Establish Covariance matrix\n");

  /* Establish the covariance matrix and invert it */

  _covmat = model->evalCovMatrixV(db, db, -1, -1);

  // Invert Covariance Matrix

  if (verbose) message("Invert Covariance matrix\n");
  if (matrix_invert(_covmat.data(),neq,-1))
  {
    messerr("Error during the covariance matrix inversion");
    return 1;
  }

  // Initialize the statistics (optional)

  _statsInit();

  return 0;
}

int GibbsUMulti::_getSize() const
{
  int nact = _getSampleRankNumber();
  int nvar = getNvar();
  return nact * nvar;
}

double GibbsUMulti::_getVariance(int iecr) const
{
  int neq  = _getSize();
  return (1. / COVMAT(iecr, iecr));
}

double GibbsUMulti::_getEstimate(int ipgs, int iecr, VectorVectorDouble& y)
{
  int nvar = getNvar();
  int nact = _getSampleRankNumber();
  int neq  = _getSize();

  double yk = 0.;
  for (int jvar = 0, jecr = 0; jvar < nvar; jvar++)
  {
    int jcase = getRank(ipgs, jvar);
    for (int jact = 0; jact < nact; jact++, jecr++)
    {
      yk -= y[jcase][jact] * COVMAT(iecr, jecr);
    }
  }
  return yk;
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
  double valsim, yk, vk;
  int nvar = getNvar();
  int nact = _getSampleRankNumber();

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
      if (!_isConstraintTight(icase, iact, &valsim))
      {
        // The term of y corresponding to the current (variable, sample)
        // is set to 0 in order to avoid testing it next.
        y[icase][iact] = 0.;

        // Calculate the estimate and the variance of estimation
        vk = _getVariance(iecr);
        yk = _getEstimate(ipgs, iecr, y) * vk;

        // Simulate the new value
        valsim = getSimulate(y, yk, sqrt(vk), icase, ipgs, ivar, iact, iter);
      }
      y[icase][iact] = valsim;
    }
  }

  // Update statistics (optional)

  _updateStats(y, ipgs, iter);
}
