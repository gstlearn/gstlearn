/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
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
#include "Matrix/LinkMatrixSparse.hpp"

#include <math.h>

// External library /// TODO : Dependency to csparse to be removed
#include "csparse_d.h"
#include "csparse_f.h"

#define WEIGHTS(ivar, jvar, iact)  (_weights[iact + nact * (jvar + ivar * nvar)])

GibbsMMulti::GibbsMMulti()
  : GibbsMulti()
  , _Ln(nullptr)
  , _Pn()
  , _eps(EPSILON6)
  , _storeTables(false)
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
  , _Ln(nullptr)
  , _Pn()
  , _eps(EPSILON6)
  , _storeTables(false)
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
  , _Ln(r._Ln)
  , _Pn(r._Pn)
  , _eps(r._eps)
  , _storeTables(r._storeTables)
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
    _Ln  = r._Ln;
    _Pn  = r._Pn;
    _eps = r._eps;
    _storeTables = r._storeTables;
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
  if (_Ln != nullptr)
    _Ln = cs_spfree(_Ln);
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
  cs*  Cmat = nullptr;
  css* S = nullptr;
  csn* N = nullptr;
  int n;
  int error = 1;
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
  _Pn.resize(nact);
  _b.resize(n);
  _x.resize(n);
  _weights.resize(nact * nvar * nvar);

  // Establish the covariance matrix as sparse (hopefully)

  if (verbose)
    message("Building Covariance Sparse Matrix (Dimension = %d)\n",nact);
  Timer timer;
  Cmat = model_covmat_by_ranks_cs(model,db,nact,_getRanks(),db,nact,_getRanks(),-1,-1);
  if (Cmat == nullptr)
  {
    messerr("Impossible to create the Total Precision Matrix");
    goto label_end;
  }
  if (verboseTimer)
    timer.displayIntervalMilliseconds("Building Covariance");

  // Cholesky decomposition

  if (verbose)
    message("Cholesky Decomposition of Covariance Matrix\n");
  S = cs_schol(Cmat, 1);
  N = cs_chol(Cmat, S);
  if (S == nullptr || N == nullptr)
  {
    messerr("Fail to perform Cholesky decomposition");
    goto label_end;
  }

  for (int iact = 0; iact < nact; iact++) _Pn[iact] = S->Pinv[iact];
  if (verboseTimer)
    timer.displayIntervalMilliseconds("Cholesky Decomposition");

  // Store the Initial Covariance in Neutral File (optional)

  if (_storeTables)
  {
    if (verbose)
      message("Storing Initial Covariance\n");
    _tableStore(1, Cmat);
    if (verboseTimer)
      timer.displayIntervalMilliseconds("Storing Initial Covariance");
  }

  // Stripping the Cholesky decomposition matrix

  _Ln = cs_strip(N->L, _eps, 3, verbose);

  // Store the reconstructed Covariance in Neutral File (optional)

  if (_storeTables)
  {
    if (verbose)
      message("Calculating Reconstructed Covariance\n");
    cs* Lt = cs_transpose(_Ln, 1);
    cs* Cmat2 = cs_multiply (_Ln, Lt);
    _tableStore(2, Cmat2);
    if (verboseTimer)
      timer.displayIntervalMilliseconds("Calculating Reconstructed Covariance");
    Lt = cs_spfree(Lt);
    Cmat2 = cs_spfree(Cmat2);
  }

  // Evaluate storage capacity and store weights

  if (verbose)
    message("Calculating and storing the weights\n");
  if (_storeAllWeights(verbose)) goto label_end;
  if (verboseTimer)
    timer.displayIntervalMilliseconds("Calculating and storing weights");

  // Initialize the statistics (optional)

  statsInit();
  error = 0;

 label_end:
  Cmat = cs_spfree(Cmat);
  S = cs_sfree(S);
  N = cs_nfree(N);
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

void GibbsMMulti::_tableStore(int mode, const cs* A)
{
  if (!A) return;
  Table table;
  const Db* db = getDb();
  const Model* model = getModel();

  // Getting maximum distance for covariance calculation

  double distmax = 1.5 * model->getCova(0)->getRange();
  double dx = db->getUnit();
  int npas = (int) ceil(distmax / dx);
  Table tabmod(npas,3);

  // Retrieve the elements from the sparse matrix

  int   n = cs_getncol(A) ;
  int* Ap = A->p ;
  int* Ai = A->i ;
  double* Ax = A->x ;

  /* Loop on the elements */

  int ecr = 0;
  for (int j = 0 ; j < n; j++)
    for (int p = Ap [j] ; p < Ap [j+1]; p++)
    {
      int iech = j;
      int jech = Ai[p];
      double value = Ax[p];

      double dist = db->getDistance(iech, jech);
      if (dist > distmax) continue;

      // Store the covariance cloud
      table.addRow();
      table.setValue(ecr, 0, dist);
      table.setValue(ecr, 1, value);
      ecr++;

      // Store the contribution to the covariance model
      int ipas = (int) floor(dist / dx + 0.5);
      if (ipas >= 0 && ipas < npas)
      {
        tabmod.add(ipas, 0, 1.);
        tabmod.add(ipas, 1, dist);
        tabmod.add(ipas, 2, value);
      }
    }

  // Serialize the table (as a Covariance Cloud)
  if (mode == 1)
    (void) table.dumpToNF("CovInit");
  else
    (void) table.dumpToNF("CovBuilt");

  // Serialize the table (as a Experimental Covariance)
  for (int ipas = 0; ipas < npas; ipas++)
  {
    double value = tabmod.getValue(ipas, 0);
    if (value > 0.)
    {
      tabmod.setValue(ipas, 1, tabmod.getValue(ipas, 1) / value);
      tabmod.setValue(ipas, 2, tabmod.getValue(ipas, 2) / value);
    }
  }

  if (mode == 1)
    (void) tabmod.dumpToNF("CovModInit");
  else
    (void) tabmod.dumpToNF("CovModBuilt");
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
    cs_ipvec(nact, _Pn.data(), _b.data(), _x.data()); /* x = P*b */
    cs_lsolve(_Ln,  _x.data()); /* x = L\x */
    cs_ltsolve(_Ln, _x.data()); /* x = L'\x */
    cs_pvec(nact, _Pn.data(), _x.data(), _b.data()); /* b = P'*x */

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
