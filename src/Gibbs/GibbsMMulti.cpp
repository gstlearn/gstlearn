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

GibbsMMulti::GibbsMMulti()
  : GibbsMulti()
  , _neigh(nullptr)
  , _wgt()
  , _flagSymNeigh(false)
  , _flagSymQ(false)
  , _flagPrintQ(false)
  , _Q(nullptr)
  , _epsilon(EPSILON6)
{
}

GibbsMMulti::GibbsMMulti(Db* db, Model* model, Neigh* neigh)
  : GibbsMulti(db, model)
  , _neigh(neigh)
  , _wgt()
  , _flagSymNeigh(false)
  , _flagSymQ(false)
  , _flagPrintQ(false)
  , _Q(nullptr)
  , _epsilon(EPSILON6)
{
}

GibbsMMulti::GibbsMMulti(const GibbsMMulti &r)
  : GibbsMulti(r)
  , _neigh(r._neigh)
  , _wgt(r._wgt)
  , _flagSymNeigh(r._flagSymNeigh)
  , _flagSymQ(r._flagSymQ)
  , _flagPrintQ(r._flagPrintQ)
  , _Q(r._Q)
  , _epsilon(r._epsilon)
{
}

GibbsMMulti& GibbsMMulti::operator=(const GibbsMMulti &r)
{
  if (this != &r)
  {
    GibbsMulti::operator=(r);
    _neigh = r._neigh;
    _wgt = r._wgt;
    _flagSymNeigh = r._flagSymNeigh;
    _flagSymQ = r._flagSymQ;
    _flagPrintQ = r._flagPrintQ;
    _Q = r._Q;
    _epsilon = r._epsilon;
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
  cs*  Q = nullptr;
  css* S = nullptr;
  csn* N = nullptr;
  double* x = nullptr;
  double* b = nullptr;
  int n;
  int nbgh_max = 0;
  int nbgh_min = ITEST;
  int error = 1;
  int iact_max = -1;
  double  delta_max = 0.;

  Db* db = getDb();
  Model* model = getModel();
  Neigh* neigh = getNeigh();
  int nvar = _getVariableNumber();
  int nact = getSampleRankNumber();
  int nech = db->getSampleNumber();
  int nvardb = db->getVariableNumber();
  bool flag_var_defined = nvardb > 0;
  if (defineGeneralNeigh(1, db, model, neigh)) return 1;

  // Consistency check

  if (flag_var_defined && nvar != nvardb)
  {
    messerr("Inconsistency in Number of Variables between Model (%d) and Db (%d)",
            nvar,nvardb);
    return 1;
  }

  // Establish the covariance matrix as sparse (hopefully)

  if (verbose) message("Building Complete Covariance Sparse Matrix\n");
  Q = model_covmat_by_ranks_cs(model,db,nech,nullptr,db,nech,nullptr,-1,-1,false,true);
  if (Q == nullptr)
  {
    messerr("Impossible to create the Total Precision Matrix");
    goto label_end;
  }

  if (verbose) message("Cholesky Decomposition of Covariance Matrix\n");
  n = Q->n;
  S = cs_schol(Q, 0);
  N = cs_chol(Q, S);
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

    // Cleaning the vector of the precision matrix Q

    ww._pivot = iact;
    for (int j = 0; j < nech; j++)
    {
      double wloc = ABS(b[j]);
      if (wloc > _epsilon)
        ww._ranks.push_back(j);
      else
        b[j] = 0.;
    }

    // Check against the identify

    cs_mulvec(Q, nech, b, x);
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

    int nbgh = static_cast<int>(ww._ranks.size());
    if (nbgh > nbgh_max) nbgh_max = nbgh;
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
    mestitle(2,"Gibbs Convergence characteristics");
    message("Epsilon = %lf\n",_epsilon);
    message("Number of weights per neighborhood is in [%d ; %d]\n",nbgh_min,nbgh_max);
    message("Maximum departure from identity (at sample #%d) = %lg\n",iact_max,delta_max);
  }

  // Improve the conditioning of the Precision resulting matrix
  // Note: the return code is not tested on purpose, to let the rest
  // of the test to be performed.

  if (verbose) message("Beautifying Kriging Weights (Symmetry, Checks, ...)\n");
  if (_improveConditioning(verbose)) return 1;

  // Initialize the statistics (optional)

  statsInit();
  error = 0;

  label_end:
  x = (double *) mem_free((char *) x);
  b = (double *) mem_free((char *) b);
  Q = cs_spfree(Q);
  S = cs_sfree(S);
  N = cs_nfree(N);
  return error;
}

int GibbsMMulti::_covmatAllocMemo(bool verbose)
{
  VectorInt ivars, iechs;
  VectorBool QFlag;
  double *covmat;

  // Initialization

  if (verbose) mestitle(1,"Gibbs using Moving Neighborhood");
  int error = 1;
  Db* db = getDb();
  Model* model = getModel();
  Neigh* neigh = getNeigh();
  int nvar = _getVariableNumber();
  int nact = getSampleRankNumber();
  int nech = db->getSampleNumber();
  int nvardb = db->getVariableNumber();
  bool flag_var_defined = nvardb > 0;
  covmat = nullptr;
  if (defineGeneralNeigh(1, db, model, neigh)) return 1;

  // Consistency check

  if (flag_var_defined && nvar != nvardb)
  {
    messerr("Inconsistency in Number of Variables between Model (%d) and Db (%d)",
            nvar,nvardb);
    return 1;
  }

  // Clear the set of weight vectors

  if (verbose) message("Establishing Neighborhoods\n");
  _wgt.resize(nact);

  QFlag.resize(nech * nech, false);
  for (int iact = 0; iact < nact; iact++)
  {
    int iech = getSampleRank(iact);

    // Select the Neighborhood
    VectorInt ranks = getGeneralNeigh(db, neigh, iech);

    // Update the neighborhood contingency matrix
    _setQFlag(QFlag, nech, iech, ranks);
  }

  // Kriging weights

  if (verbose) message("Establishing Sets of Kriging Weights\n");
  for (int iact = 0; iact < nact; iact++)
  {
    int iech = getSampleRank(iact);
    GibbsWeights& ww = _wgt[iact];

    // Read the neighborhood from the contingency table

    ww._ranks = _getQFlag(QFlag, nech, iech);
    int nbgh  = static_cast<int>(ww._ranks.size());
    int neq   = nvar * nbgh;

    // Establishing the (moving) Covariance matrix
    covmat = model_covmat_by_ranks(model,
                                   db, neq, ww._ranks.data(),
                                   db, neq, ww._ranks.data(),
                                   -1, -1, 0, 1);
    if (covmat == nullptr) goto label_end;

    // Inverting the (moving) Covariance matrix
    if (matrix_invert(covmat, neq, 0)) goto label_end;

    // Find the rank of the pivot (in absolute sample number)

    int found = -1;
    for (int jech = 0; jech < nbgh && found < 0; jech++)
      if (iech == ww._ranks[jech]) found = jech;
    ww._pivot = found;

    // Modify ranks into internal numbering

    for (int jech = 0; jech < nbgh; jech++)
    {
      ww._ranks[jech] = getRelativeRank(ww._ranks[jech]);
      if (ww._ranks[jech] < 0) goto label_end;
    }

    // Store the weights per variable for the target sample

    for (int ivar = 0; ivar < nvar; ivar++)
    {
      VectorDouble ll(neq);
      for (int i = 0; i < neq; i++)
      {
        ll[i] = COVMAT(i, ww._pivot + ivar * nbgh);
      }
      ww._ll.push_back(ll);
    }
    covmat = (double *) mem_free((char *) covmat);
  }

  // Improve the conditioning of the Precision resulting matrix
  // Note: the return code is not tested on purpose, to let the rest
  // of the test to be performed.

  if (verbose) message("Beautifying Kriging Weights (Symmetry, Checks, ...)\n");
  if (_improveConditioning(verbose)) return 1;

  // Initialize the statistics (optional)

  statsInit();
  error = 0;

  label_end:
  (void) defineGeneralNeigh(-1, db, model, neigh);
  covmat = (double *) mem_free((char * ) covmat);
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

int GibbsMMulti::_improveConditioning(bool verbose)
{
  bool err_def, err_sym;

  // Construct the Precision matrix Q

  if (_buildQ()) return 1;

  // Make the Matrix symmetric

  if (_flagSymQ)
    _makeQSymmetric(_Q);

  // Optional printout

  if (_flagPrintQ)
    cs_print_nice("Reconstructed Precision Matrix", _Q, -1, -1);

  // Check that the matrix is symmetric

  err_sym = ! cs_isSymmetric(_Q, verbose);

  // Check that Q is definite-positive

  err_def = ! cs_isDefinitePositive(_Q,  verbose);

  // Summarize errors

  if (err_sym || err_def) return 1;

  // Update the newly modified weights extracted from Q

  _extractWeightFromQ();

  return 0;
}

void GibbsMMulti::_print(int iact) const
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

int GibbsMMulti::_getVariableNumber() const
{
  Model* model = getModel();
  return model->getVariableNumber();
}

void GibbsMMulti::_setQFlag(VectorBool& QFlag,
                            int nech,
                            int iech,
                            const VectorInt& ranks) const
{
  int nbgh = static_cast<int>(ranks.size());
  for (int ibgh = 0; ibgh < nbgh; ibgh++)
  {
    QFLAG(iech, ranks[ibgh]) = 1;
    if (_flagSymNeigh) QFLAG(ranks[ibgh], iech) = 1;
  }
}

VectorInt GibbsMMulti::_getQFlag(VectorBool& QFlag, int nech, int iech)
{
  VectorInt ranks;

  for (int jech = 0; jech < nech; jech++)
  {
    if (QFLAG(iech,jech)) ranks.push_back(jech);
  }
  return ranks;
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
