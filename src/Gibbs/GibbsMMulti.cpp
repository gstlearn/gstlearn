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

#define WEIGHTS(ivar, jvar, iact)  (_weights[iact + nact * (jvar + ivar * nvar)])

GibbsMMulti::GibbsMMulti()
  : GibbsMulti()
  , _Cmat(nullptr)
  , _eps(EPSILON6)
  , _hdf5("Gibbs.hdf5","GibbsSet")
  , _flagStoreInternal(true)
  , _b()
  , _x()
  , _weights()
  , _areas()
{
}

GibbsMMulti::GibbsMMulti(Db* db, Model* model)
  : GibbsMulti(db, model)
  , _Cmat(nullptr)
  , _eps(EPSILON6)
  , _hdf5("Gibbs.hdf5","GibbsSet")
  , _flagStoreInternal(true)
  , _b()
  , _x()
  , _weights()
  , _areas()
{
}

GibbsMMulti::GibbsMMulti(const GibbsMMulti &r)
  : GibbsMulti(r)
  , _Cmat(r._Cmat)
  , _eps(r._eps)
  , _hdf5(r._hdf5)
  , _flagStoreInternal(r._flagStoreInternal)
  , _b(r._b)
  , _x(r._x)
  , _weights(r._weights)
  , _areas(r._areas)
{
}

GibbsMMulti& GibbsMMulti::operator=(const GibbsMMulti &r)
{
  if (this != &r)
  {
    GibbsMulti::operator=(r);
    _Cmat  = r._Cmat;
    _eps = r._eps;
    _hdf5 = r._hdf5;
    _flagStoreInternal = r._flagStoreInternal;
    _b = r._b;
    _x = r._x;
    _weights = r._weights;
    _areas = r._areas;
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
  MatrixSparse*  Cmat = nullptr;
  int n;
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

  // Core allocation

  n = nact * nvar;
  _b.resize(n);
  _x.resize(n);
  _weights.resize(nact * nvar * nvar);

  // Establish the covariance matrix as sparse (hopefully)

  if (verbose)
    message("Building Covariance Sparse Matrix (Dimension = %d)\n",nact);
  Timer timer;
  _Cmat = model_covmat_by_ranks_cs(model,db,nact,_getRanks(),db,nact,_getRanks(),-1,-1);
  if (_Cmat == nullptr) return 1;
  if (verboseTimer)
    timer.displayIntervalMilliseconds("Building Covariance");

  // Cholesky decomposition

  if (verbose)
    message("Cholesky Decomposition of Covariance Matrix\n");
  if (Cmat->computeCholesky())
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

  statsInit();
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
void GibbsMMulti::update(VectorVectorDouble& y,
                         int isimu,
                         int ipgs,
                         int iter)
{
  double valsim, yk, vark;

  int nvar = _getVariableNumber();
  int nact = getSampleRankNumber();

  /* Print the title */

  if (OptDbg::query(EDbg::CONVERGE))
    mestitle(1, "Gibbs Sampler (Simu:%d - GS:%d)", isimu + 1, ipgs + 1);

  /* Loop on the target */

  for (int iact = 0; iact < nact; iact++)
  {
    // Load the vector of weights
    _getWeights(iact);

    for (int ivar = 0; ivar < nvar; ivar++)
    {
      int icase = getRank(ipgs,ivar);
      if (! isConstraintTight(ipgs, ivar, iact, &valsim))
      {
        // Calculate the estimate and the variance of estimation
        // at point 'íact' for target variable 'ivar'
        _getEstimate(ipgs, ivar, iact, icase, y, &yk, &vark);

        // Simulate the new value
        valsim = getSimulate(y, yk, sqrt(vark), ipgs, ivar, iact, iter);
      }
      y[icase][iact] = valsim;
    }
  }

  // Update statistics (optional)

  updateStats(y, ipgs, iter);
}

int GibbsMMulti::_getVariableNumber() const
{
  Model* model = getModel();
  return model->getVariableNumber();
}

/**
 * Calculate the set of (multivariate) weights for one given sample
 * @param iact0   Rank of the sample
 * @param tol     Tolerance below which weights are set to 0
 * @return Number of zero elements in the current vector
 */
int GibbsMMulti::_calculateWeights(int iact0, double tol) const
{
  int nvar  = _getVariableNumber();
  int nact  = getSampleRankNumber();
  int n     = nact * nvar;
  int nzero = 0;

  // Loop on the variables

  for (int ivar0 = 0; ivar0 < nvar; ivar0++)
  {
    for (int i = 0; i < n; i++) _b[i] = 0.;
    _b[iact0 + ivar0 * nact] = 1.;

    // Solve the linear system and returns the result in 'x'
    _Cmat->solveCholesky(_b, _x);

    // Discarding the values leading to small vector of weights

    for (int i = 0; i < n; i++)
    {
      double wloc = ABS(_b[i]);
      if (wloc < tol) _b[i] = 0.;
      if (ABS(_b[i]) < EPSILON10) nzero++;
    }

    // Storing the weights for the current sample and the current variable

    double *w_loc = &WEIGHTS(ivar0,0,0);
    double *b_loc = &_b[0];
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int iact = 0; iact < nact; iact++)
      {
        *w_loc = *b_loc;
        w_loc++;
        b_loc++;
      }
  }

  return nzero;
}

int GibbsMMulti::_storeAllWeights(bool verbose)
{
  VectorDouble weights;
  int nvar   = _getVariableNumber();
  int nact   = getSampleRankNumber();
  _areas.clear();

  // Decide if weights are stored internally or not

  long ntotal = nact * nvar * nact * nvar;
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
    dims[1] = nvar * nact * nvar;
    _hdf5.openNewFile("h5data3.h5");
    _hdf5.openNewDataSetDouble("Set3", 2, dims.data());
#endif
  }

  // Loop on the samples

  int nremain = 0;
  for (int iact = 0; iact < nact; iact++)
  {
    nremain += _calculateWeights(iact);

    if (_flagStoreInternal)
    {
      // Store internally
      _areas.push_back(_weights);
    }
    else
    {
      // Store in hdf5 file
      _hdf5.writeDataDoublePartial(iact, _weights);
    }
  }

  // Optional printout of the number of zero weights

  if (verbose)
  {
    message("- Number of zero weights = %d\n",nremain);
    double reduc = 100. * (double) (ntotal - nremain) / (double) ntotal;
    message("- Percentage             = %6.2lf\n",reduc);
  }
  return 0;
}

void GibbsMMulti::_getWeights(int iact0) const
{
  if (! _areas.empty())
  {
    // Load from the internal storage
    _weights = _areas[iact0];
  }
  else
  {
    // Read from the external file
    _weights = _hdf5.getDataDoublePartial(iact0);
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

void GibbsMMulti::_getEstimate(int ipgs0,
                               int ivar0,
                               int iact0,
                               int icase,
                               VectorVectorDouble& y,
                               double *yk_arg,
                               double *vark_arg) const
{
  int nvar = _getVariableNumber();
  int nact = getSampleRankNumber();

  double yk   = 0.;
  double vark = 1. / WEIGHTS(ivar0, ivar0, iact0);

  // The term of y corresponding to the current (variable, sample)
  // is set to 0 in order to avoid testing it in the next loop
  y[icase][iact0] = 0.;

  VectorDouble::iterator wloc = _weights.begin() + nact * nvar * ivar0;
  for (int jvar = 0; jvar < nvar; jvar++)
  {
    int jcase = getRank(ipgs0, jvar);
    const VectorDouble& yloc = y[jcase];
    for (int jact = 0; jact < nact; jact++)
    {
      yk -= yloc[jact] * (*wloc);
      wloc++;
    }
  }

  // Setting returned arguments
  *yk_arg = yk * vark;
  *vark_arg = vark;
}
