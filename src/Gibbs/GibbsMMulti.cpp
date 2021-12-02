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
#include "Gibbs/GibbsMMulti.hpp"
#include "Gibbs/AGibbs.hpp"
#include "Model/Model.hpp"
#include "Basic/Law.hpp"
#include "Basic/Timer.hpp"
#include "Basic/HDF5format.hpp"
#include "Morpho/Morpho.hpp"
#include "Db/Db.hpp"
#include "Covariances/CovAniso.hpp"
#include "csparse_f.h"
#include "geoslib_f.h"
#include "geoslib_old_f.h"

#include <math.h>

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
** \param[in]  verbose     Verbose flag
**
*****************************************************************************/
int GibbsMMulti::covmatAlloc(bool verbose)
{
  // Initialization

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
  int nvardb = db->getVariableNumber();
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
  Cmat = model_covmat_by_ranks_cs(model,db,nact,_getRanks(),db,nact,_getRanks(),-1,-1,false,true);
  if (Cmat == nullptr)
  {
    messerr("Impossible to create the Total Precision Matrix");
    goto label_end;
  }
  timer.Interval("Building Covariance");

  // Cholesky decomposition

  if (verbose) message("Cholesky Decomposition of Covariance Matrix\n");
  n = Cmat->n;
  S = cs_schol(Cmat, 0);
  N = cs_chol(Cmat, S);
  if (S == nullptr || N == nullptr)
  {
    messerr("Fail to perform Cholesky decomposition");
    goto label_end;
  }
  for (int iact = 0; iact < nact; iact++) _Pn[iact] = S->Pinv[iact];
  timer.Interval("Cholesky Decomposition");

  // Store the Initial Covariance in Neutral File (optional)

  if (_storeTables)
  {
    _tableStore(1, Cmat);
    timer.Interval("Storing Initial Covariance");
  }

  // Stripping the Cholesky decomposition matrix

  _Ln = cs_strip(N->L, _eps, verbose);
  timer.Interval("Stripping Cholesky Matrix");

  // Evaluate storage capacity and store weights internally (or not)

  if (_storeAllWeights(verbose)) goto label_end;

  // Store the reconstructed Covariance in Neutral File (optional)

  if (_storeTables)
  {
    cs* Lt = cs_transpose(_Ln, 1);
    cs* Cmat2 = cs_multiply (_Ln, Lt);
    Lt = cs_spfree(Lt);
    _tableStore(2, Cmat2);
    Cmat2 = cs_spfree(Cmat2);
    timer.Interval("Storing Reconstructed Covariance");
  }

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
  double valsim;

  int nvar = _getVariableNumber();
  int nact = getSampleRankNumber();

  /* Print the title */

  if (debug_query("converge"))
    mestitle(1, "Gibbs Sampler (Simu:%d - GS:%d)", isimu + 1, ipgs + 1);

  /* Loop on the target */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int icase = getRank(ipgs,ivar);
    for (int iact = 0; iact < nact; iact++)
    {
      if (! isConstraintTight(ipgs, ivar, iact, &valsim))
      {
        _getWeights(iact);
        double yk   = 0.;
        double vark = 1. / WEIGHTS(ivar, ivar, iact);
        for (int jvar = 0; jvar < nvar; jvar++)
        {
          int jcase = getRank(ipgs, jvar);
          for (int jact = 0; jact < nact; jact++)
          {
            if (ivar != jvar || iact != jact)
              yk -= y[jcase][jact] * WEIGHTS(ivar, jvar, jact);
          }
        }
        yk *= vark;
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
  double dx = db->getDX(0);
  int npas = static_cast<int>(ceil(distmax / dx));
  Table tabmod(npas,2);

  // Retrieve the elements from the sparse matrix

  int   n = A->n ;
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
      table.resize(ecr+1, 2);
      table.update(ecr, 0, dist);
      table.update(ecr, 1, value);
      ecr++;

      // Store the contribution to the covariance model
      int ipas = static_cast<int>(floor(dist + dx/2.) / dx);
      if (ipas >= 0 && ipas < npas)
      {
        tabmod.increment(ipas, 0, 1.);
        tabmod.increment(ipas, 1, value);
      }
    }

  // Serialize the table (as a Covariance Cloud)
  if (mode == 1)
    table.serialize("CovInit", false);
  else
    table.serialize("CovBuilt", false);

  // Serialize the table (as a Experimental Covariance)
  for (int ipas = 0; ipas < npas; ipas++)
  {
    double value = tabmod.getValue(ipas, 0);
    if (value > 0.) value = tabmod.getValue(ipas, 1) / value;
    tabmod.update(ipas, 1, value);
  }

  if (mode == 1)
    tabmod.serialize("CovModInit", false);
  else
    tabmod.serialize("CovModBuilt", false);
}

/**
 * Calculate the set of (multivariate) weights for one given sample
 * @param iact0   Rank of the sample
 * @param tol     Tolerance below which weights are set to 0
 * @return
 */
int GibbsMMulti::_calculateWeights(int iact0, double tol) const
{
  int nvar = _getVariableNumber();
  int nact = getSampleRankNumber();
  int n    = nact * nvar;

  // Loop on the variables

  for (int ivar0 = 0; ivar0 < nvar; ivar0++)
  {
    for (int i = 0; i < n; i++) _b[i] = 0.;
    _b[iact0 + ivar0 * nact] = 1.;

    // Solve the linear system and returns the result in 'x'
    cs_ipvec(nact, _Pn.data(), _b.data(), _x.data()); /* x = P*b */
    cs_lsolve(_Ln, _x.data()); /* x = L\x */
    cs_ltsolve(_Ln, _x.data()); /* x = L'\x */
    cs_pvec(nact, _Pn.data(), _x.data(), _b.data()); /* b = P'*x */

    // Discarding the values leading to small vector of weights

    for (int i = 0; i < n; i++)
    {
      double wloc = ABS(_b[i]);
      if (wloc < tol) _b[i] = 0.;
    }

    // Storing the weights for the current sample and the current variable

    for (int ivar = 0; ivar < nvar; ivar++)
      for (int iact = 0; iact < nact; iact++)
        WEIGHTS(ivar0, ivar, iact) = _b[ivar * nact + iact];
  }

  return 0;
}

int GibbsMMulti::_storeAllWeights(bool verbose)
{
  VectorDouble weights;
  int nvar   = _getVariableNumber();
  int nact   = getSampleRankNumber();
  _areas.clear();

  // Check if the weight must be stored internally or externally

  long total = nact * nvar * nact * nvar;
  if (verbose)
  {
    message("Total core needs                = %ld (bytes)\n",total);
  }

  // Decide if weights are stored internally or not

  if (verbose)
  {
    if (! _flagStoreInternal)
      message("Weights are stored externally (HDF5 format)\n");
    else
      message("Weights are stored internally\n");
  }

  // Create the HDF5 file (optional)

  if (! _flagStoreInternal)
  {
    std::vector<hsize_t> dims(2);
    dims[0] = nact;
    dims[1] = nvar * nact * nvar;
    _hdf5.openNewFile("h5data3.h5");
    _hdf5.openNewDataSet("Set3", 2, dims.data(), H5::PredType::NATIVE_DOUBLE);
  }

  // Loop on the samples

  for (int iact = 0; iact < nact; iact++)
  {
    if (_calculateWeights(iact)) return 1;

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
