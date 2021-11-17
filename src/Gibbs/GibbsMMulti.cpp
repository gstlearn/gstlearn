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
  , _wgt()
  , _Q(nullptr)
  , _epsilon1(EPSILON6)
  , _epsilon2(EPSILON6)
{
}

GibbsMMulti::GibbsMMulti(Db* db, Model* model)
  : GibbsMulti(db, model)
  , _wgt()
  , _Q(nullptr)
  , _epsilon1(EPSILON6)
  , _epsilon2(EPSILON6)
{
}

GibbsMMulti::GibbsMMulti(const GibbsMMulti &r)
  : GibbsMulti(r)
  , _wgt(r._wgt)
  , _Q(r._Q)
  , _epsilon1(r._epsilon1)
  , _epsilon2(r._epsilon2)
{
}

GibbsMMulti& GibbsMMulti::operator=(const GibbsMMulti &r)
{
  if (this != &r)
  {
    GibbsMulti::operator=(r);
    _wgt = r._wgt;
    _Q = r._Q;
    _epsilon1 = r._epsilon1;
    _epsilon2 = r._epsilon2;
  }
  return *this;
}

GibbsMMulti::~GibbsMMulti()
{
  if (_Q != nullptr)
    _Q = cs_spfree(_Q);
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
  double* x = nullptr;
  double* b = nullptr;
  int n;
  int nbgh_max = ITEST;
  int nbgh_min = ITEST;
  int error = 1;
  int iact_max = -1;
  int nstrip = 0;
  double delta_max = 0.;
  double q_min = TEST;
  double q_max = 0.;
  Db* db = getDb();
  Model* model = getModel();
  int nvar = _getVariableNumber();
  int nact = getSampleRankNumber();
  int nech = db->getSampleNumber();
  int nvardb = db->getVariableNumber();
  bool flag_var_defined = nvardb > 0;

  // Consistency check

  if (flag_var_defined && nvar != nvardb)
  {
    messerr("Inconsistency in Number of Variables between Model (%d) and Db (%d)",
            nvar,nvardb);
    return 1;
  }

  // Establish the covariance matrix as sparse (hopefully)

  if (verbose) message("Building Complete Covariance Sparse Matrix\n");
  Cmat = model_covmat_by_ranks_cs(model,db,nech,nullptr,db,nech,nullptr,-1,-1,false,true);
  if (Cmat == nullptr)
  {
    messerr("Impossible to create the Total Precision Matrix");
    goto label_end;
  }


  if (verbose) message("Cholesky Decomposition of Covariance Matrix\n");
  n = Cmat->n;
  S = cs_schol(Cmat, 0);
  N = cs_chol(Cmat, S);
  if (S == nullptr || N == nullptr)
  {
    messerr("Fail to perform Cholesky decomposition");
    goto label_end;
  }

  // Core allocation

  x = (double *) mem_alloc(sizeof(double) * n,0);
  if (x == nullptr) goto label_end;
  b = (double *) mem_alloc(sizeof(double) * n,0);
  if (b == nullptr) goto label_end;

  // Kriging weights

  if (verbose) message("Establishing Sets of Kriging Weights\n");
  _wgt.resize(nact);

  // Optional printout

  if (verbose)
  {
    mestitle(2,"Gibbs Convergence characteristics");
    message("Epsilon (Cholesky) = %lf\n",_epsilon2);
    message("Epsilon (weight) = %lf\n",_epsilon1);
  }

  // Stripping the Cholesky decomposition matrix

  cs_strip(N->L, _epsilon2, verbose);

  // Store experimental and approximate covariances in the table

  _tableStore(db, Cmat, N, verbose);

  for (int iact = 0; iact < nact; iact++)
  {
    int iech = getSampleRank(iact);
    GibbsWeights& ww = _wgt[iact];

    for (int j = 0; j < n; j++) b[j] = 0.;
    b[iech] = 1.;

    // Solve the linear system and returns the result in 'x'
    cs_ipvec(nech, S->Pinv, b, x); /* x = P*b */
    cs_lsolve(N->L, x); /* x = L\x */
    cs_ltsolve(N->L, x); /* x = L'\x */
    cs_pvec(nech, S->Pinv, x, b); /* b = P'*x */

    // Discarding the values leading to small weights

    ww._ranks = VectorInt(nech);
    int nbgh = 0;
    for (int j = 0; j < nech; j++)
    {
      double wloc = ABS(b[j]);
      if (FFFF(q_min) || wloc < q_min) q_min = wloc;
      if (FFFF(q_max) || wloc > q_max) q_max = wloc;
      if (wloc > _epsilon1)
      {
        ww._ranks[nbgh] = j;
        if (iact == j) ww._pivot = nbgh;
        nbgh++;
      }
      else
      {
        b[j] = 0.;
        nstrip++;
      }
    }
    ww._ranks.resize(nbgh);

    // Check against the identify

    cs_mulvec(Cmat, nech, b, x);
    for (int j = 0; j < nech; j++)
    {
      double ref = (j == iact) ? 1. : 0.;
      double delta = ABS(x[j] - ref);
      if (delta > delta_max)
      {
        delta_max = delta;
        iact_max = iact;
      }
    }

    // Storing the weights per variable for the target sample

    if (IFFFF(nbgh_min) || nbgh > nbgh_max) nbgh_max = nbgh;
    if (IFFFF(nbgh_min) || nbgh < nbgh_min) nbgh_min = nbgh;
    int neq  = nvar * nbgh;
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      VectorDouble ll(neq);
      for (int i = 0; i < nbgh; i++)
        ll[i] = b[ivar * nbgh + ww._ranks[i]];
      ww._ll.push_back(ll);
    }
  }

  if (verbose)
  {
    message("Range of Q values = [%lf ; %lf]\n", q_min,q_max);
    message("Number of weights stripped off = %d\n",nstrip);
    message("Number of weights per neighborhood is in [%d ; %d]\n",nbgh_min,nbgh_max);
    message("Maximum departure from identity (at sample #%d) = %lg\n",iact_max,delta_max);
  }

  // Improve the conditioning of the Precision resulting matrix
  // Note: the return code is not tested on purpose, to let the rest
  // of the test to be performed.

  if (_testConditioning(verbose)) return 1;
  if (DEBUG) _display();

  // Initialize the statistics (optional)

  statsInit();
  error = 0;

 label_end:
  x = (double *) mem_free((char *) x);
  b = (double *) mem_free((char *) b);
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

  int nact = getSampleRankNumber();
  int nvar = _getVariableNumber();

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
        const GibbsWeights& ww = _wgt[iact];
        int nbgh  = static_cast<int>(ww._ranks.size());
        int pivot = ww._pivot;

        /* Loop on the Data */

        double yk = 0.;
        double vark = 1. / ww._ll[ivar][pivot + nbgh * ivar];
        for (int jvar = 0; jvar < nvar; jvar++)
        {
          int jcase = getRank(ipgs, jvar);
          for (int jbgh = 0; jbgh < nbgh; jbgh++)
          {
            int jact = ww._ranks[jbgh];
            if (ivar != jvar || iact != jact)
              yk -= y[jcase][jact] * ww._ll[jvar][jbgh + nbgh * jvar];
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

int GibbsMMulti::_testConditioning(bool verbose)
{

  // Construct the Precision matrix Q (symmetric by construction)

  if (_buildQ()) return 1;

  // Check that the matrix is symmetric

  if (! cs_isSymmetric(_Q, verbose)) return 1;

  // Check that Q is definite-positive

  if (! cs_isDefinitePositive(_Q,  verbose)) return 1;

  // Update the newly modified weights extracted from Q

  _extractWeightFromQ();

  return 0;
}

void GibbsMMulti::_display(int iact) const
{
  int nvar = _getVariableNumber();

  message("Active sample #%d\n",iact);
  const GibbsWeights& ww = _wgt[iact];
  message("Position within the neighboring sample list: %d\n",ww._pivot);

  int nbgh = static_cast<int>(ww._ranks.size());
  print_ivector("Ranks",0,nbgh,ww._ranks);
  for (int ivar=0; ivar < nvar; ivar++)
    print_vector("Weights",0,nbgh,ww._ll[ivar]);
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

void GibbsMMulti::_makeQSymmetric(cs* Q) const
{
  int n = Q->n;

  // Set the lower part of the matrix to symmetrical values
  for (int irow = 0; irow < n; irow++)
    for (int icol = 0; icol < irow; icol++)
    {
      double val1 = cs_get_value(Q, irow, icol);
      double val2 = cs_get_value(Q, icol, irow);
      cs_set_value(Q, irow, icol, (val1 + val2) / 2.);
    }

  // Copy the upper part by symmetry

  for (int irow = 0; irow < n; irow++)
    for (int icol = irow; icol < n; icol++)
    {
      double val = cs_get_value(Q, icol, irow);
      cs_set_value(Q, irow, icol, val);
    }
}

int GibbsMMulti::_buildQ()
{
  double* mat = nullptr;
  double *eigval = nullptr;
  double *eigvec = nullptr;
  int n;
  cs* T = nullptr;
  int nvar = _getVariableNumber();
  int nact = getSampleRankNumber();

  // Constitute the triplet

  T = cs_spalloc(0, 0, 1, 1, 1);
  if (T == nullptr) return 1;

  // Create partial precision matrix Q from the weights

  for (int ivar = 0; ivar < nvar; ivar++)
  for (int iact = 0; iact < nact; iact++)
  {
    const GibbsWeights& ww = _wgt[iact];
    int nbgh  = static_cast<int>(ww._ranks.size());
    int icol = iact + nbgh * ivar;

    for (int jbgh = 0; jbgh < nbgh; jbgh++)
    {
      int jact = ww._ranks[jbgh];
      for (int jvar = 0; jvar < nvar; jvar++)
      {
        int jcol = jact + nbgh * jvar;
        double wgt = ww._ll[jvar][jbgh + nbgh * jvar];
        if (! cs_entry(T, icol, jcol, wgt)) goto label_end;
      }
    }
  }

  // Convert from triplet to sparse matrix

  _Q = cs_triplet(T);
  n = _Q->n;

  // Get the Eigen values

  mat = cs_toArray(_Q);
  eigval = (double *) mem_alloc(sizeof(double) * n,1);
  eigvec = (double *) mem_alloc(sizeof(double) * n * n,1);
  matrix_eigen(mat, n, eigval, eigvec);
  message("Range of Eigen values = [%lf ; %lf]\n",eigval[n-1],eigval[0]);
  mat = (double *) mem_free((char * ) mat);
  eigval = (double *) mem_free((char *) eigval);
  eigvec = (double *) mem_free((char *) eigvec);

  label_end:
  T = cs_spfree(T);
  return 0;
}

void GibbsMMulti::_extractWeightFromQ()
{
  int nvar = _getVariableNumber();
  int nact = getSampleRankNumber();

  for (int ivar = 0; ivar < nvar; ivar++)
    for (int iact = 0; iact < nact; iact++)
    {
      GibbsWeights& ww = _wgt[iact];
      int nbgh = static_cast<int>(ww._ranks.size());
      int icol = iact + nbgh * ivar;

      for (int jbgh = 0; jbgh < nbgh; jbgh++)
      {
        int jact = ww._ranks[jbgh];
        for (int jvar = 0; jvar < nvar; jvar++)
        {
          int jcol = jact + nbgh * jvar;
          double value = cs_get_value(_Q, icol, jcol);
          ww._ll[jvar][jbgh + nbgh * jvar] = value;
        }
      }
    }
}

void GibbsMMulti::_tableStore(const Db* db, const cs* Cmat, const csn* N, bool verbose)
{
  int nech = db->getSampleNumber();
  Table table(nech * nech, 3);

  // Constitute the approximate covariance matrix (from modified Cholesky)

  cs* L = N->L;
  cs* Lt = cs_transpose(L, 1);
  cs* Cmat2 = cs_multiply (L, Lt);
  Lt = cs_spfree(Lt);

  int ecr = 0;
  double ecart = 0.;
  double ecart_max = 0.;
  int iech_max = -1;
  int jech_max = -1;
  for (int iech = 0; iech < nech; iech++)
    for (int jech = 0; jech < nech; jech++, ecr++)
    {
      table.update(ecr, 0, db->getDistance(iech, jech));
      table.update(ecr, 1, cs_get_value(Cmat, iech, jech));
      table.update(ecr, 2, cs_get_value(Cmat2, iech, jech));
      double delta = ABS(cs_get_value(Cmat,iech,jech) - cs_get_value(Cmat2, iech, jech));
      if (delta > ecart_max)
      {
        iech_max = iech;
        jech_max = jech;
        ecart_max = delta;
      }
      ecart += delta * delta;
    }

  ecart /= (double) (nech * nech);
  if (verbose)
  {
    message("Departure for Covariance Matrix between initial and reconstituted\n");
    message("- Maximum = %lf (Initial=%lf - Reconstituted=%lf\n",
            ecart_max,cs_get_value(Cmat,iech_max,jech_max),cs_get_value(Cmat2,iech_max,jech_max));
    message("- Average = %lf\n",ecart);
  }

  // Serialize the table
  table.serialize("Covariance", false);
}
