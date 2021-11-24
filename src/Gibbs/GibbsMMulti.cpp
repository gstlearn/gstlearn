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
#include "Model/Model.hpp"
#include "Basic/Law.hpp"
#include "Basic/Timer.hpp"
#include "Morpho/Morpho.hpp"
#include "csparse_f.h"
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "csparse_f.h"

#define COVMAT(i,j)              (covmat[(i) * neq + (j)])
#define QFLAG(iech,jech)         (QFlag[(iech) * nech + jech])

#define DEBUG 0
#define TOTAL_MAX 50000000

GibbsMMulti::GibbsMMulti()
  : GibbsMulti()
  , _Ln(nullptr)
  , _Pn()
  , _eps(EPSILON6)
  , _storeTables(false)
  , _b()
  , _x()
  , _areas()
{
}

GibbsMMulti::GibbsMMulti(Db* db, Model* model)
  : GibbsMulti(db, model)
  , _Ln(nullptr)
  , _Pn()
  , _eps(EPSILON6)
  , _storeTables(false)
  , _b()
  , _x()
  , _areas()
{
}

GibbsMMulti::GibbsMMulti(const GibbsMMulti &r)
  : GibbsMulti(r)
  , _Ln(r._Ln)
  , _Pn(r._Pn)
  , _eps(r._eps)
  , _storeTables(r._storeTables)
  , _b(r._b)
  , _x(r._x)
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
    _b = r._b;
    _x = r._x;
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
  int nech   = getSampleNumber();
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

  n = nech * nvar;
  _Pn.resize(nech);
  _b.resize(n);
  _x.resize(n);

  // Establish the covariance matrix as sparse (hopefully)

  if (verbose)
    message("Building Covariance Sparse Matrix (Dimension = %d)\n",nech);
  Timer timer;
  Cmat = model_covmat_by_ranks_cs(model,db,nech,nullptr,db,nech,nullptr,-1,-1,false,true);
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
  for (int i = 0; i < nech; i++) _Pn[i] = S->Pinv[i];
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

  _storeAllWeights(verbose);

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

  // Display the weights (conditional to DEBUG variable)

  if (DEBUG) _display();

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
  WgtVect area;
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
      int iech = getSampleRank(iact);
      if (! isConstraintTight(ipgs, ivar, iact, &valsim))
      {
        _getWeights(iech, area);
        int nbgh    = area._nbgh;
        int pivot   = area._pivot;
        double yk   = 0.;
        double vark = 1. / area._weights[ivar][pivot + nbgh * ivar];
        for (int jvar = 0; jvar < nvar; jvar++)
        {
          int jcase = getRank(ipgs, jvar);
          for (int jbgh = 0; jbgh < nbgh; jbgh++)
          {
            int jact = area._mvRanks[jbgh];
            if (ivar != jvar || iact != jact)
              yk -= y[jcase][jact] * area._weights[jvar][jbgh + nbgh * jvar];
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

void GibbsMMulti::_display(int iech) const
{
  WgtVect area;
  int nvar = _getVariableNumber();

  _getWeights(iech, area);

  message("Active sample #%d\n",iech);
  message("Position within the neighboring sample list: %d\n",area._pivot);

  print_ivector("Ranks",0,area._nbgh,area._mvRanks);
  for (int ivar=0; ivar < nvar; ivar++)
    print_vector("Weights", 0, area._nbgh, area._weights[ivar]);
}

void GibbsMMulti::_display() const
{
  int nact = getSampleRankNumber();
  mestitle(1,"Weight in Moving Gibbs");
  for (int iact = 0; iact < nact; iact++)
    _display(iact);
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
  int npas = ceil(distmax / dx);
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
      int ipas = floor(dist + dx/2.) / dx;
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

int GibbsMMulti::_calculateWeights(int iech, WgtVect& area, double tol) const
{
  int pivot = -1;
  int nvar  = _getVariableNumber();
  int nech  = getSampleNumber();
  int n = nech * nvar;
  for (int i = 0; i < n; i++) _b[i] = 0.;
  _b[iech] = 1.;

  // Solve the linear system and returns the result in 'x'
  cs_ipvec(nech, _Pn.data(), _b.data(), _x.data()); /* x = P*b */
  cs_lsolve(_Ln, _x.data()); /* x = L\x */
  cs_ltsolve(_Ln, _x.data()); /* x = L'\x */
  cs_pvec(nech, _Pn.data(), _x.data(), _b.data()); /* b = P'*x */

  // Discarding the values leading to small vector of weights

  int nbgh = 0;
  area._mvRanks.reserve(nech);
  for (int j = 0; j < nech; j++)
  {
    double wloc = ABS(_b[j]);
    if (wloc > tol)
    {
      area._mvRanks[nbgh] = j;
      if (iech == j) pivot = nbgh;
      nbgh++;
    }
  }
  area._mvRanks.resize(nbgh);
  if (area._mvRanks.empty()) return 1;

  // Storing the weights per variable for the target sample

  area._weights.clear();
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    VectorDouble ll(nbgh);
    for (int i = 0; i < nbgh; i++)
      ll[i] = _b[ivar * nbgh + area._mvRanks[i]];
    area._weights.push_back(ll);
    if (area._weights.empty()) return 1;
  }
  area._pivot = pivot;
  area._nbgh  = nbgh;
  area._size  = _getSizeOfArea(area);

  return 0;
}

bool GibbsMMulti::_checkForInternalStorage(bool verbose)
{
  WgtVect area;
  int nvar   = _getVariableNumber();
  int nact   = getSampleRankNumber();

  // Evaluate the sum of the number of weights per sample

  long sumWgt   = nact * nvar * nact * nvar;
  long sumRanks = nact * nact;
  long overHead = 1;
  for (int iact = 0; iact < nact ; iact++)
  {
    overHead += 2 + sizeof(VectorInt) + sizeof(VectorDouble);
    overHead += nvar * sizeof(VectorDouble);
  }

  long total = sumWgt * sizeof(double) + sumRanks * sizeof(int) + overHead;

  // Printout of the different core storage quanta

  if (verbose)
  {
    message("Total core for Weights          = %ld (bytes)\n",sumWgt * sizeof(double));
    message("Total core for Weights          = %ld (bytes)\n",sumRanks * sizeof(int));
    message("Core for Overhead               = %ld (bytes)\n",overHead);
    message("Total core needs                = %ld (bytes)\n",total);
    message("(To speed up, calculations are based on an upper bound of dimension of weights)\n");
  }

  // Decide if weights are stored internally or not

  bool flagStoreInternal;
  if (verbose)
    message("Decision Memory Threshold       = %ld (bytes)\n",TOTAL_MAX);
  if (total > TOTAL_MAX)
  {
    flagStoreInternal = false;
    if (verbose)
      message("Weights are not stored internally: they will be calculated on the fly\n");
  }
  else
  {
    flagStoreInternal = true;
    if (verbose)
      message("Weights are calculated only once and stored internally\n");
  }

  return flagStoreInternal;
}

void GibbsMMulti::_storeAllWeights(bool verbose)
{
  WgtVect area;
  int nact = getSampleRankNumber();
  _areas.clear();

  // Check if the weight must be stored internally or calculated on the fly
  if (! _checkForInternalStorage(verbose)) return;

  // Loop on the samples

  for (int iact = 0; iact < nact; iact++)
  {
    int iech = getSampleRank(iact);

    if (_calculateWeights(iech, area))
    {
      _areas.clear();
      return;
    }
    _areas.push_back(area);
  }
}

void GibbsMMulti::_getWeights(int iech, WgtVect& area) const
{
  if (! _areas.empty())
    area = _areas[iech];
  else
    _calculateWeights(iech, area);
}

/**
 * Calculate the size of the current WgtVect structure
 * @param area Current WgtVect structure
 * @return
 */
int GibbsMMulti::_getSizeOfArea(const WgtVect& area) const
{
  int s1 = sizeof(area._size) + sizeof(area._nbgh) + sizeof(area._pivot);
  int s2 = ut_vector_size(area._mvRanks);
  int s3 = ut_vector_size(area._weights);
  int total = s1 + s2 + s3;
  return total;
}
