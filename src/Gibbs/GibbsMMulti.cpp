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

#include "Gibbs/GibbsMMulti.hpp"
#include "Gibbs/AGibbs.hpp"
#include "Model/Model.hpp"
#include "Basic/Law.hpp"
#include "Basic/Timer.hpp"
#include "Basic/HDF5format.hpp"
#include "Basic/OptDbg.hpp"
#include "Morpho/Morpho.hpp"
#include "Db/Db.hpp"
#include "Covariances/CovAniso.hpp"
#include "Matrix/MatrixSparse.hpp"

#include <math.h>

#define WEIGHTS(ivar, jvar, iact)  (weights[iact + nact * (jvar + ivar * nvar)])

GibbsMMulti::GibbsMMulti()
  : GibbsMulti()
  , _Cmat(nullptr)
  , _eps(EPSILON6)
  , _flagStoreInternal(true)
  , _areas()
  , _hdf5("Gibbs.hdf5","GibbsSet")
{
}

GibbsMMulti::GibbsMMulti(Db* db, Model* model)
  : GibbsMulti(db, model)
  , _Cmat(nullptr)
  , _eps(EPSILON6)
  , _flagStoreInternal(true)
  , _areas()
  , _hdf5("Gibbs.hdf5","GibbsSet")
{
}

GibbsMMulti::GibbsMMulti(const GibbsMMulti &r)
  : GibbsMulti(r)
  , _Cmat(r._Cmat)
  , _eps(r._eps)
  , _flagStoreInternal(r._flagStoreInternal)
  , _areas(r._areas)
  , _hdf5(r._hdf5)
{
}

GibbsMMulti& GibbsMMulti::operator=(const GibbsMMulti &r)
{
  if (this != &r)
  {
    GibbsMulti::operator=(r);
    _Cmat  = r._Cmat;
    _eps = r._eps;
    _flagStoreInternal = r._flagStoreInternal;
    _areas = r._areas;
    _hdf5 = r._hdf5;
  }
  return *this;
}

GibbsMMulti::~GibbsMMulti()
{
  if (_Cmat != nullptr)
    delete _Cmat;
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
int GibbsMMulti::covmatAlloc(bool verbose, bool verboseTimer)
{
  // Initialization

  if (verboseTimer) verbose = true;
  if (verbose) mestitle(1,"Gibbs using Moving Neighborhood");
  Db* db = getDb();
  Model* model = getModel();
  int nvar   = _getVariableNumber();
  int nact   = getSampleRankNumber();
  int nvardb = db->getLocNumber(ELoc::Z);
  bool flag_var_defined = nvardb > 0;

  // Consistency check

  if (flag_var_defined && nvar != nvardb)
  {
    messerr("Inconsistency in Number of Variables between Model (%d) and Db (%d)",
            nvar,nvardb);
    return 1;
  }

  // Establish the covariance matrix as sparse (hopefully)

  if (verbose)
    message("Building Covariance Sparse Matrix (Dimension = %d)\n",nact);
  Timer timer;
  _Cmat = model_covmat_by_ranks_Mat(model,db,nact,_getRanks(),db,nact,_getRanks(),-1,-1,
                                    nullptr, EPSILON3, verbose);
  if (_Cmat == nullptr) return 1;
  if (verboseTimer)
    timer.displayIntervalMilliseconds("Building Covariance");

  // Cholesky decomposition

  if (verbose)
    message("Cholesky Decomposition of Covariance Matrix\n");
  if (_Cmat->computeCholesky())
  {
    messerr("Fail to perform Cholesky decomposition");
    return 1;
  }
  if (verboseTimer)
    timer.displayIntervalMilliseconds("Cholesky Decomposition");

  // Evaluate storage capacity and store weights

  if (verbose)
    message("Calculating and storing the weights\n");
  if (_storeAllWeights(verbose)) return 1;
  if (verboseTimer)
    timer.displayIntervalMilliseconds("Calculating and storing weights");

  // Initialize the statistics (optional)

  _statsInit();
  return 0;
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
void GibbsMMulti::update(VectorVectorDouble &y, int isimu, int ipgs, int iter)
{
  double valsim, yk, vark;

  int nvar = _getVariableNumber();
  int nact = getSampleRankNumber();
  VectorDouble weights(nact * nvar * nvar);

  /* Print the title */

  if (OptDbg::query(EDbg::CONVERGE))
    mestitle(1, "Gibbs Sampler (Simu:%d - GS:%d)", isimu + 1, ipgs + 1);

  /* Loop on the target */

  for (int iact0 = 0; iact0 < nact; iact0++)
  {
    // Load the vector of weights
    _getWeights(iact0, weights);

    for (int ivar0 = 0; ivar0 < nvar; ivar0++)
    {
      int icase = getRank(ipgs,ivar0);
      if (! _isConstraintTight(ipgs, ivar0, iact0, &valsim))
      {
        // The term of y corresponding to the current (variable, sample)
        // is set to 0 in order to avoid testing it in the next loop
        y[icase][iact0] = 0.;

        // Calculate the estimate and the variance of estimation
        vark = 1. / WEIGHTS(ivar0, ivar0, iact0);
        yk = _getEstimate(ipgs, ivar0, iact0, y, weights) * vark;

        // Simulate the new value
        valsim = getSimulate(y, yk, sqrt(vark), ipgs, ivar0, iact0, iter);
      }

      y[icase][iact0] = valsim;
    }
  }

  // Update statistics (optional)

  _updateStats(y, ipgs, iter);
}

int GibbsMMulti::_getVariableNumber() const
{
  Model* model = getModel();
  return model->getVariableNumber();
}

/**
 * Calculate the set of (multivariate) weights for one given sample
 * @param iact0   Rank of the sample
 * @param b       Right-hand side vector
 * @param x       Vector of solution
 * @param weights Vector of weights
 * @param tol     Tolerance below which weights are set to 0
 * @return Number of zero elements in the current vector
 */
int GibbsMMulti::_calculateWeights(int iact0,
                                   VectorDouble &b,
                                   VectorDouble &x,
                                   VectorDouble &weights,
                                   double tol) const
{
  int nvar  = _getVariableNumber();
  int nact  = getSampleRankNumber();
  int n     = nact * nvar;
  int nzero = 0;

  // Loop on the variables

  for (int ivar0 = 0; ivar0 < nvar; ivar0++)
  {
    b.fill(0.);
    b[iact0 + ivar0 * nact] = 1.;

    // Solve the linear system and returns the result in 'x'
    _Cmat->solveCholesky(b, x);

    // Discarding the values leading to small vector of weights
    for (int i = 0; i < n; i++)
    {
      double wloc = ABS(x[i]);
      if (wloc < tol) x[i] = 0.;
      if (ABS(x[i]) < EPSILON10) nzero++;
    }

    // Storing the weights for the current sample and the current variable
    double *w_loc = &WEIGHTS(ivar0,0,0);
    double *x_loc = &x[0];
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int iact = 0; iact < nact; iact++)
      {
        *w_loc = *x_loc;
        w_loc++;
        x_loc++;
      }
  }
  return nzero;
}

int GibbsMMulti::_storeAllWeights(bool verbose)
{
  int nvar = _getVariableNumber();
  int nact = getSampleRankNumber();
  int n    = nact * nvar;
  VectorDouble b(n);
  VectorDouble x(n);
  VectorDouble weights(n * nvar);
  _areas.clear();

  // Decide if weights are stored internally or not

  long ntotal = n * n;
  if (verbose)
  {
    if (! _flagStoreInternal)
      message("Weights are stored externally (HDF5 format)\n");
    else
      message("Weights are stored internally\n");
    message("- Total core needs       = %ld\n",ntotal);
  }

  // Create the HDF5 file (optional)

  if (! _flagStoreInternal)
  {
#ifdef _USE_HDF5
    std::vector<hsize_t> dims(2);
    dims[0] = nact;
    dims[1] = nvar * n;
    _hdf5.openNewFile("h5data3.h5");
    _hdf5.openNewDataSetDouble("Set3", 2, dims.data());
#endif
  }

  // Loop on the samples

  int nzero = 0;
  for (int iact = 0; iact < nact; iact++)
  {
    nzero += _calculateWeights(iact, b, x, weights);
    _storeWeights(iact, weights);
  }

  // Optional printout of the number of zero weights

  if (verbose)
  {
    message("- Number of zero weights = %d\n",nzero);
    double reduc = 100. * (double) (ntotal - nzero) / (double) ntotal;
    message("- Percentage             = %6.2lf\n",reduc);
  }
  return 0;
}

/**
 * Storing the weights when processing the current sample
 * @param iact0 Rank of the current sample
 * @param weights Vector of weights to be stored
 */
void GibbsMMulti::_storeWeights(int iact0, const VectorDouble& weights)
{
  if (_flagStoreInternal)
  {
    // Store internally
    _areas.push_back(weights);
  }
  else
  {
    // Store in hdf5 file
    _hdf5.writeDataDoublePartial(iact0, weights);
  }
}

void GibbsMMulti::_getWeights(int iact0, VectorDouble& weights) const
{
  if (_flagStoreInternal)
  {
    // Load from the internal storage
    weights = _areas[iact0];
  }
  else
  {
    // Read from the external file
    weights = _hdf5.getDataDoublePartial(iact0);
  }
}

void GibbsMMulti::cleanup()
{
  _hdf5.closeDataSet();
  _hdf5.closeFile();
  _hdf5.deleteFile();
}

void GibbsMMulti::setFlagStoreInternal(bool flagStoreInternal)
{
#ifndef _USE_HDF5
  if (!flagStoreInternal)
    messerr("No HDF5 support: Cannot use External Storing of weights option!");
  _flagStoreInternal = true;
#else
  _flagStoreInternal = flagStoreInternal;
#endif
}

double GibbsMMulti::_getEstimate(int ipgs0,
                                 int ivar0,
                                 int iact0,
                                 const VectorVectorDouble &y,
                                 const VectorDouble& weights) const
{
  int nvar = _getVariableNumber();
  int nact = getSampleRankNumber();

  double yk   = 0.;
  const double *w_loc = &WEIGHTS(ivar0,0,0);
  for (int jvar = 0; jvar < nvar; jvar++)
  {
    int jcase = getRank(ipgs0, jvar);
    const VectorDouble& yloc = y[jcase];
    for (int jact = 0; jact < nact; jact++)
    {
      yk -= yloc[jact] * (*w_loc);
      w_loc++;
    }
  }
  return yk;
}
