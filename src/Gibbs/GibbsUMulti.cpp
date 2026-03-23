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
#include "Gibbs/GibbsUMulti.hpp"
#include "Basic/OptDbg.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "geoslib_old_f.h"

#include <cmath>

#define COVMAT(i, j) (_covmat[(i) * neq + (j)])

namespace gstlrn
{
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

  GibbsUMulti::GibbsUMulti(const GibbsUMulti& r)
    : GibbsMulti(r)
    , _covmat(r._covmat)
  {
  }

  GibbsUMulti& GibbsUMulti::operator=(const GibbsUMulti& r)
  {
    if (this != &r)
    {
      GibbsMulti::operator=(r);
      _covmat = r._covmat;
    }
    return *this;
  }

  GibbsUMulti::~GibbsUMulti() {}

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
  Id GibbsUMulti::covmatAlloc(bool verbose, bool /*verboseTimer*/)
  {
    // Initialization

    if (verbose) mestitle(1, "Gibbs using Unique Neighborhood");
    Db* db = getDb();
    Model* model = getModel();

    // Establish Covariance Matrix

    if (verbose) message("Establish Covariance matrix\n");

    /* Establish the covariance matrix and invert it */

    _covmat = model->evalCovMat(db, db, -1, -1);

    // Invert Covariance Matrix

    if (verbose) message("Invert Covariance matrix\n");
    _covmat.invert();

    // Initialize the statistics (optional)

    _statsInit();

    return 0;
  }

  Id GibbsUMulti::_getSize() const
  {
    auto nact = _getSampleRankNumber();
    auto nvar = getNvar();
    return nact * nvar;
  }

  double GibbsUMulti::_getVariance(Id iecr) const
  {
    return (1. / _covmat.getValue(iecr, iecr));
  }

  double GibbsUMulti::_getEstimate(Id ipgs, Id iecr, VectorVectorDouble& y)
  {
    auto nvar = getNvar();
    auto nact = _getSampleRankNumber();

    double yk = 0.;
    for (Id jvar = 0, jecr = 0; jvar < nvar; jvar++)
    {
      auto jcase = getRank(ipgs, jvar);
      for (Id jact = 0; jact < nact; jact++, jecr++)
      {
        yk -= y[jcase][jact] * _covmat.getValue(iecr, jecr);
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
  void GibbsUMulti::update(VectorVectorDouble& y, Id isimu, Id ipgs, Id iter)
  {
    double valsim, yk, vk;
    auto nvar = getNvar();
    auto nact = _getSampleRankNumber();

    /* Print the title */

    if (OptDbg::query(EDbg::CONVERGE))
      mestitle(1,
               "Iterative Conditional Expectation (GS:%d - Simu:%d)",
               ipgs + 1,
               isimu + 1);

    /* Loop on the target */

    for (Id ivar = 0, iecr = 0; ivar < nvar; ivar++)
    {
      auto icase = getRank(ipgs, ivar);
      for (Id iact = 0; iact < nact; iact++, iecr++)
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
} // namespace gstlrn
