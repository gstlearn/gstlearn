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
#include "Gibbs/GibbsUMultiMono.hpp"
#include "Basic/OptDbg.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "geoslib_old_f.h"

#include <cmath>

#define COVMAT(ivar, i, j) (_covmat[ivar][(i) * nact + (j)])

namespace gstlrn
{
  GibbsUMultiMono::GibbsUMultiMono()
    : GibbsMultiMono()
    , _covmat()
  {
  }

  GibbsUMultiMono::GibbsUMultiMono(
    Db* db,
    const std::vector<Model*>& models,
    double rho)
    : GibbsMultiMono(db, models, rho)
    , _covmat()
  {
  }

  GibbsUMultiMono::GibbsUMultiMono(const GibbsUMultiMono& r)
    : GibbsMultiMono(r)
    , _covmat(r._covmat)
  {
  }

  GibbsUMultiMono& GibbsUMultiMono::operator=(const GibbsUMultiMono& r)
  {
    if (this != &r)
    {
      GibbsMultiMono::operator=(r);
      _covmat = r._covmat;
    }
    return *this;
  }

  GibbsUMultiMono::~GibbsUMultiMono() {}

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
  Id GibbsUMultiMono::covmatAlloc(bool verbose, bool /*verboseTimer*/)
  {
    Db* db = getDb();

    // Initialization

    if (verbose)
      mestitle(1, "Gibbs using Unique Neighborhood in MultiMono case");
    auto nvar = getNVar();
    _covmat.resize(nvar);

    // Loop on the variables

    for (Id ivar = 0; ivar < nvar; ivar++)
    {
      Model* model = getModels(ivar);

      // Establish Covariance Matrix (always based on the first variable in MultiMono case)

      if (verbose) message("Establish Covariance matrix (Var=%d)\n", ivar + 1);
      _covmat[ivar] = model->evalCovMat(db, db, 0, 0);

      // Invert Covariance Matrix

      if (verbose) message("Invert Covariance matrix (Var=%d)\n", ivar + 1);
      _covmat[ivar].invert();
    }

    // Initialize the statistics (optional)

    _statsInit();

    return 0;
  }

  double GibbsUMultiMono::_getVariance(Id ivar, Id iact) const
  {
    return 1. / _covmat[ivar].getValue(iact, iact);
  }

  double GibbsUMultiMono::_getEstimate(
    Id icase,
    Id ivar,
    Id iact,
    VectorVectorDouble& y) const
  {
    auto nact = _getSampleRankNumber();

    double yk = 0.;
    for (Id jact = 0; jact < nact; jact++)
    {
      yk -= y[icase][jact] * _covmat[ivar].getValue(iact, jact);
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
  void
    GibbsUMultiMono::update(VectorVectorDouble& y, Id isimu, Id ipgs, Id iter)
  {
    double valsim;
    auto nact = _getSampleRankNumber();
    auto nvar = getNvar();

    /* Print the title */

    if (OptDbg::query(EDbg::CONVERGE))
      mestitle(
        1, "Iterative Conditional Expectation (PGS=%d - Simu:%d - Iter=%d)",
        ipgs + 1, isimu + 1, iter + 1);

    /* Loop on the target */

    for (Id ivar = 0; ivar < nvar; ivar++)
    {
      auto icase = getRank(ipgs, ivar);
      for (Id iact = 0; iact < nact; iact++)
      {
        if (!_isConstraintTight(icase, iact, &valsim))
        {
          // The term of y corresponding to the current (variable, sample)
          // is set to 0 in order to avoid testing it next.
          y[icase][iact] = 0.;

          // Calculate the estimate and the variance of estimation
          double vk = _getVariance(ivar, iact);
          double yk = _getEstimate(icase, ivar, iact, y) * vk;

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
