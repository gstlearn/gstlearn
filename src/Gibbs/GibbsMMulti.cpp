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
#include "Gibbs/GibbsMMulti.hpp"
#include "Gibbs/AGibbs.hpp"
#include "Model/Model.hpp"
#include "Basic/Timer.hpp"
#include "Basic/HDF5format.hpp"
#include "Basic/OptDbg.hpp"
#include "Db/Db.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/NF_Triplet.hpp"

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <math.h>

static bool storeSparse = true;

GibbsMMulti::GibbsMMulti()
  : GibbsMulti()
  , _Cmat(nullptr)
  , _eps(EPSILON6)
  , _flagStoreInternal(true)
  , _areas()
  , _hdf5("Gibbs.hdf5","GibbsSet")
  , _matWgt()
  , _weights()
{
  _allocate();
}

GibbsMMulti::GibbsMMulti(Db* db, Model* model)
  : GibbsMulti(db, model)
  , _Cmat(nullptr)
  , _eps(EPSILON6)
  , _flagStoreInternal(true)
  , _areas()
  , _hdf5("Gibbs.hdf5","GibbsSet")
  , _matWgt()
  , _weights()
{
  _allocate();
}

GibbsMMulti::GibbsMMulti(const GibbsMMulti &r)
  : GibbsMulti(r)
  , _Cmat(r._Cmat)
  , _eps(r._eps)
  , _flagStoreInternal(r._flagStoreInternal)
  , _areas(r._areas)
  , _hdf5(r._hdf5)
  , _matWgt(r._matWgt)
  , _weights()
{
  _allocate();
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
    _matWgt = r._matWgt;
    _weights = r._weights;
  }
  return *this;
}

GibbsMMulti::~GibbsMMulti()
{
  delete _Cmat;
  delete _matWgt;
}

void GibbsMMulti::_allocate()
{
  _weights.resize(_getSize());
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
  int nact   = _getSampleRankNumber();
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
  _Cmat = model->evalCovMatrixSparse(db, db, -1, -1, _getRanks(), _getRanks());
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

double GibbsMMulti::_getVariance(int icol) const
{
  if (storeSparse) return (1. / _matWgt->getValue(icol, icol));
  return (1. / _weights[icol]);
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
  double valsim, yk, vk;

  int nvar = _getVariableNumber();
  int nact = _getSampleRankNumber();

  /* Print the title */

  if (OptDbg::query(EDbg::CONVERGE))
    mestitle(1, "Gibbs Sampler (Simu:%d - GS:%d)", isimu + 1, ipgs + 1);

  /* Loop on the target */

  for (int ivar0 = 0; ivar0 < nvar; ivar0++)
  {
    int icase = getRank(ipgs,ivar0);
    for (int iact0 = 0; iact0 < nact; iact0++)
    {
      // Load the vector of weights
      int icol = _getColumn(iact0, ivar0);
      _getWeights(icol);

      if (! _isConstraintTight(icase, iact0, &valsim))
      {
        // The term of y corresponding to the current (variable, sample)
        // is set to 0 in order to avoid testing it next.
        y[icase][iact0] = 0.;

        // Calculate the estimate and the variance of estimation
        vk = _getVariance(icol);
        yk = _getEstimate(ipgs, icol, y) * vk;

        // Simulate the new value
        valsim = getSimulate(y, yk, sqrt(vk), icase, ipgs, ivar0, iact0, iter);
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

int GibbsMMulti::_getSize() const
{
  int nact = _getSampleRankNumber();
  int nvar = _getVariableNumber();
  return nact * nvar;
}

int GibbsMMulti::_getColumn(int iact, int ivar) const
{
  int nact = _getSampleRankNumber();
  return (iact + ivar * nact);
}

void GibbsMMulti::_splitCol(int icol, int *iact, int *ivar) const
{
  int nact = _getSampleRankNumber();
  *ivar = icol / nact;
  *iact = icol - nact * (*ivar);
}

/**
 * Calculate the set of (multivariate) weights for a given sample / variable
 * @param icol    Rank of the column of interest
 * @param b       Right-hand side vector
 * @param tol     Tolerance below which weights are set to 0
 */
void GibbsMMulti::_calculateWeights(int icol,
                                    VectorDouble &b,
                                    double tol) const
{
  b.fill(0.);
  b[icol] = 1.;

  // Solve the linear system and returns the result in 'x'
  _Cmat->solveCholesky(b, _weights);

  if (tol <= 0.) return;

  // Discarding the values leading to small vector of weights
  for (int irow = 0, nrow = _getSize(); irow < nrow; irow++)
  {
    double xloc = ABS(_weights[irow]);
    if (xloc < tol) _weights[irow] = 0.;
  }
}

void GibbsMMulti::_updateStatWeights(int* nzero)
{
  for (int irow = 0, nrow = _getSize(); irow < nrow; irow++)
  {
    double wgt = _weights[irow];
    if (isZero(wgt)) (*nzero)++;
  }
}

int GibbsMMulti::_storeAllWeights(bool verbose)
{
  int nrow = _getSize();
  VectorDouble b(nrow);

  // Decide if weights are stored internally or not

  long ntotal = nrow * nrow;
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
    dims[0] = nrow;
    dims[1] = nrow;
    _hdf5.openNewFile("h5data3.h5");
    _hdf5.openNewDataSetDouble("Set3", 2, dims.data());
#endif
  }

  // Loop on the samples

  int nzero = 0;
  if (storeSparse)
  {
    _matWgt = new MatrixSparse();
    NF_Triplet NF_T;
    for (int icol = 0, ncol = _getSize(); icol < ncol; icol++)
    {
      _calculateWeights(icol, b);
      _updateStatWeights(&nzero);
      _storeWeightsMS(icol, NF_T);
    }
    _matWgt->resetFromTriplet(NF_T);
  }
  else
  {
    _areas.clear();
    for (int icol = 0, ncol = _getSize(); icol < ncol; icol++)
    {
      _calculateWeights(icol, b);
      _updateStatWeights(&nzero);
      _storeWeights(icol);
    }
  }

  // Optional printout of the number of zero weights

  if (verbose)
  {
    message("- Number of zero weights = %d\n",  nzero);
    double reduc = 100. * (double) (ntotal - nzero) / (double) ntotal;
    message("- Percentage             = %6.2lf\n",reduc);
  }
  return 0;
}

/**
 * Storing the weights when processing the current sample
 * @param icol  Rank of the column of interest
 */
void GibbsMMulti::_storeWeights(int icol)
{
  if (_flagStoreInternal)
  {
    // Store internally
    _areas.push_back(_weights);
  }
  else
  {
    // Store in hdf5 file
#ifdef _USE_HDF5
    _hdf5.writeDataDoublePartial(icol, _weights);
#else
    DECLARE_UNUSED(icol);
#endif
  }
}

void GibbsMMulti::_storeWeightsMS(int icol, NF_Triplet& NF_T)
{
  for (int irow = 0, nrow = _getSize(); irow < nrow; irow++)
    if (ABS(_weights[irow]) > EPSILON10) NF_T.add(irow, icol, _weights[irow]);
}

void GibbsMMulti::_getWeights(int icol) const
{
  if (_flagStoreInternal)
  {
    if (storeSparse) return;
    // Load from the internal storage
    _weights = _areas[icol];
  }
  else
  {
    // Read from the external file
    _weights = HDF5format::getDataDoublePartial(icol);
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

double GibbsMMulti::_getEstimate(int ipgs,
                                 int icol,
                                 const VectorVectorDouble &y) const
{
  int jact, jvar, jcase;
  double yk = 0.;

  if (storeSparse)
  {
    if (_matWgt->isFlagEigen())
    {
      for (Eigen::SparseMatrix<double>::InnerIterator it(_matWgt->getEigenMatrix(),icol); it; ++it)
      {
        _splitCol(it.row(), &jact, &jvar);
        jcase = getRank(ipgs, jvar);
        yk -= y[jcase][jact] * it.value();
      }
    }
  }
  else
  {
    int nvar = _getVariableNumber();
    int nact = _getSampleRankNumber();
    int irow = 0;
    for (jvar = 0; jvar < nvar; jvar++)
    {
      jcase = getRank(ipgs, jvar);
      for (jact = 0; jact < nact; jact++, irow++)
      {
        yk -= y[jcase][jact] * _weights[irow];
      }
    }
  }
  return yk;
}
