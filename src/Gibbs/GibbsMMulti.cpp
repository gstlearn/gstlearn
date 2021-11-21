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
#include "Db/Db.hpp"
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

GibbsMMulti::GibbsMMulti()
  : GibbsMulti()
  , _Ln(nullptr)
  , _Pn()
  , _epsilon1(EPSILON6)
  , _epsilon2(EPSILON6)
  , _storeTables(false)
  , _ranks()
  , _b()
  , _x()
{
}

GibbsMMulti::GibbsMMulti(Db* db, Model* model)
  : GibbsMulti(db, model)
  , _Ln(nullptr)
  , _Pn()
  , _epsilon1(EPSILON6)
  , _epsilon2(EPSILON6)
  , _storeTables(false)
  , _ranks()
  , _b()
  , _x()
{
}

GibbsMMulti::GibbsMMulti(const GibbsMMulti &r)
  : GibbsMulti(r)
  , _Ln(r._Ln)
  , _Pn(r._Pn)
  , _epsilon1(r._epsilon1)
  , _epsilon2(r._epsilon2)
  , _storeTables(r._storeTables)
  , _ranks(r._ranks)
  , _b(r._b)
  , _x(r._x)
{
}

GibbsMMulti& GibbsMMulti::operator=(const GibbsMMulti &r)
{
  if (this != &r)
  {
    GibbsMulti::operator=(r);
    _Ln  = r._Ln;
    _Pn  = r._Pn;
    _epsilon1 = r._epsilon1;
    _epsilon2 = r._epsilon2;
    _storeTables = r._storeTables;
    _ranks = r._ranks;
    _b = r._b;
    _x = r._x;
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
  _ranks.resize(nech);
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

  _Ln = cs_strip(N->L, _epsilon2, verbose);
  timer.Interval("Stripping Cholesky Matrix");

  // Evaluate storage capacity

  _storageEvaluate();

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
  double valsim;

  int pivot, nbgh;

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
        VectorVectorDouble wgt = _getWeights(iech,&nbgh,&pivot);
        double yk   = 0.;
        double vark = 1. / wgt[ivar][pivot + nbgh * ivar];
        for (int jvar = 0; jvar < nvar; jvar++)
        {
          int jcase = getRank(ipgs, jvar);
          for (int jbgh = 0; jbgh < nbgh; jbgh++)
          {
            int jact = _ranks[jbgh];
            if (ivar != jvar || iact != jact)
              yk -= y[jcase][jact] * wgt[jvar][jbgh + nbgh * jvar];
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
  int nbgh, pivot;
  int nvar = _getVariableNumber();

  VectorVectorDouble wgt = _getWeights(iech,&nbgh,&pivot);

  message("Active sample #%d\n",iech);
  message("Position within the neighboring sample list: %d\n",pivot);

  print_ivector("Ranks",0,nbgh,_ranks);
  for (int ivar=0; ivar < nvar; ivar++)
    print_vector("Weights",0,nbgh,wgt[ivar]);
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

VectorVectorDouble GibbsMMulti::_getWeights(int  iech,
                                            int *nbgh_arg,
                                            int *pivot_arg) const
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

  // Discarding the values leading to small weights

  int nbgh = 0;
  for (int j = 0; j < nech; j++)
  {
    double wloc = ABS(_b[j]);
    if (wloc > _epsilon1)
    {
      _ranks[nbgh] = j;
      if (iech == j) pivot = nbgh;
      nbgh++;
    }
  }

  // Storing the weights per variable for the target sample

  VectorVectorDouble weights;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    VectorDouble ll(nbgh);
    for (int i = 0; i < nbgh; i++)
      ll[i] = _b[ivar * nbgh + _ranks[i]];
    weights.push_back(ll);
  }

  *pivot_arg = pivot;
  *nbgh_arg = nbgh;
  return weights;
}

void GibbsMMulti::_storageEvaluate()
{
  int nbgh, pivot;
  int nvar   = _getVariableNumber();
  int nech   = getSampleNumber();
  int nact   = getSampleRankNumber();

  // Evaluate the sum of the number of weights per sample

  int total = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    for (int iact = 0; iact < nact; iact++)
    {
      int iech = getSampleRank(iact);
      VectorVectorDouble wgt = _getWeights(iech,&nbgh,&pivot);
      total += nbgh;
    }
  }

  int total_size =
      nech * sizeof(std::vector<std::vector<VectorDouble> >) +
      nvar * nech * sizeof(std::vector<VectorDouble>) *
      nech * nvar * nech * sizeof(VectorDouble) +
      total * sizeof(double);

  message("Total size for storage = %d / %d\n",total_size,total);
}
