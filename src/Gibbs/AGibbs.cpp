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
#include "Gibbs/AGibbs.hpp"
#include "Model/Model.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Law.hpp"
#include "geoslib_old_f.h"
#include "geoslib_define.h"
#include "geoslib_enum.h"

AGibbs::AGibbs()
    : _npgs(1),
      _nvar(1),
      _nburn(1),
      _niter(1),
      _flagOrder(true),
      _flagMultiMono(true),
      _flagDecay(true),
      _flagStats(false),
      _rho(1.),
      _sqr(0.),
      _eps(EPSILON3),
      _ranks(),
      _db(nullptr),
      _model(nullptr),
      _stats()
{
}

AGibbs::AGibbs(Db* db, Model* model)
    : _npgs(1),
      _nvar(1),
      _nburn(1),
      _niter(1),
      _flagOrder(true),
      _flagMultiMono(true),
      _flagDecay(true),
      _flagStats(false),
      _rho(1.),
      _sqr(0.),
      _eps(EPSILON3),
      _ranks(),
      _db(db),
      _model(model),
      _stats()
{
}

AGibbs::AGibbs(Db* db, Model* model,
               int npgs, int nvar, int nburn, int niter,
               int flag_order, bool flag_multi_mono, bool flag_decay,
               double rho, double eps)
    : _npgs(1),
      _nvar(1),
      _nburn(1),
      _niter(1),
      _flagOrder(false),
      _flagMultiMono(true),
      _flagDecay(true),
      _flagStats(false),
      _rho(1.),
      _sqr(0.),
      _eps(eps),
      _ranks(),
      _db(db),
      _model(model),
      _stats()
{
  init(npgs, nvar, nburn, niter,
       flag_order, flag_multi_mono, flag_decay, rho, eps);
}

AGibbs::AGibbs(const AGibbs &r)
    : _npgs(r._npgs),
      _nvar(r._nvar),
      _nburn(r._nburn),
      _niter(r._niter),
      _flagOrder(r._flagOrder),
      _flagMultiMono(r._flagMultiMono),
      _flagDecay(r._flagDecay),
      _flagStats(r._flagStats),
      _rho(r._rho),
      _sqr(r._sqr),
      _eps(r._eps),
      _ranks(r._ranks),
      _db(r._db),
      _model(r._model),
      _stats(r._stats)
{
}

AGibbs& AGibbs::operator=(const AGibbs &r)
{
  if (this != &r)
  {
    _npgs = r._npgs;
    _nvar = r._nvar;
    _nburn = r._nburn;
    _niter = r._niter;
    _flagOrder = r._flagOrder;
    _flagMultiMono = r._flagMultiMono;
    _flagDecay = r._flagDecay;
    _flagStats = r._flagStats;
    _rho = r._rho;
    _sqr = r._sqr;
    _eps = r._eps;
    _ranks = r._ranks;
    _db = r._db;
    _model = r._model;
    _stats = r._stats;
  }
  return *this;
}

AGibbs::~AGibbs()
{
}

void AGibbs::init(int npgs,
                  int nvar,
                  int nburn,
                  int niter,
                  int flag_order,
                  bool flag_multi_mono,
                  bool flag_decay,
                  double rho,
                  double eps)
{
  _npgs = npgs;
  _nvar = nvar;
  _nburn = nburn;
  _niter = niter;
  _flagOrder = flag_order;
  _flagMultiMono = flag_multi_mono;
  _flagDecay = flag_decay;
  _rho = rho;
  _eps = eps;

  // In the real multivariate scheme, the value of '_rho' is not significant
  if (_flagMultiMono) _rho = 1.;
  _sqr = sqrt(1. - _rho * rho);

  // Evaluate the array of active sample ranks
  _ranks = calculateSampleRanks();
}

/****************************************************************************/
/*!
**  Check/Show the facies against gaussian at wells
**
** \return Error return code
**
** \param[in]  y          Gaussian vector
** \param[in]  isimu      Rank of the simulation
** \param[in]  ipgs       Rank of the GS
**
*****************************************************************************/
int AGibbs::checkGibbs(const VectorVectorDouble& y,
                       int isimu,
                       int ipgs)
{
  int nactive = _db->getActiveSampleNumber();
  int nvar    = getNvar();
  mestitle(1,"Checking gaussian values from Gibbs vs. bounds (PGS=%d Simu=%d)",
           ipgs+1,isimu+1);

  int nerror = 0;

  /* Loop on the variables */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int icase   = getRank(ipgs,ivar);
    int icase0  = getRank(ipgs,0);

    /* Loop on the data */

    for (int iact=0; iact<nactive; iact++)
    {
      int iech = getSampleRank(iact);
      double vmin = _db->getLowerBound(iech,icase);
      double vmax = _db->getUpperBound(iech,icase);
      if (FFFF(vmin)) vmin = -1.e30;
      if (FFFF(vmax)) vmax =  1.e30;

      /* Read the gaussian value */

      double gaus = y[icase][iact];
      if (ivar > 0)
        gaus = _sqr * gaus + _rho * y[icase0][iact];

      /* Check inconsistency */

      if ((! FFFF(vmin) && gaus < vmin) ||
          (! FFFF(vmax) && gaus > vmax))
      {
        message("- Sample (#%d):",iech+1);
        message(" Simu#%d of Y%d=%lf",isimu+1,ivar+1,gaus);
        message(" does not lie within [");
        if (FFFF(vmin))
          message("NA,");
        else
          message("%lf",vmin);
        message(";");
        if (FFFF(vmax))
         message("NA");
        else
          message("%lf",vmax);
        message("]\n");
        nerror++;
      }
    }
  }

  if (nerror <= 0) message("No problem found\n");

  return nerror;
}

/****************************************************************************/
/*!
**  Correct the bounds according to the order relationship
**
** \return  Error code: 1 if there is no solution; 0 otherwise
**
** \param[in]  iech          Absolute rank of the sample
** \param[in]  ipgs          Rank of the GS
** \param[in]  ivar          Rank of the variable
**
** \param[out]  vmin_arg   Output minimum bound
** \param[out]  vmax_arg   Output maximum bound
**
** \remark Attributes LOC_GAUSFAC are mandatory
**
*****************************************************************************/
int AGibbs::_boundsCheck(int iech,
                         int ipgs,
                         int ivar,
                         double *vmin_arg,
                         double *vmax_arg)
{
  const Db* db = getDb();
  int icase = getRank(ipgs, ivar);
  double vmin = db->getLowerBound(iech,icase);
  double vmax = db->getUpperBound(iech,icase);

  if (! FFFF(vmin) && ! FFFF(vmax) && (vmin) > (vmax))
  {
    messerr("Sample %d: Bounds are wrongly ordered: Vmin(%lf) > Vmax(%lf)",
            iech+1,vmin,vmax);
    return 1;
  }

  *vmin_arg = vmin;
  *vmax_arg = vmax;
  return 0;
}

/****************************************************************************/
/*!
**  Print the inequality
**
** \param[in]  iact      Relative rank of the sample
** \param[in]  ivar      Rank of the variable
** \param[in]  nfois     Rank of the iteration (<=0 for bootstrap)
** \param[in]  flag_cv   1 to print the convergence criterion
** \param[in]  simval    Simulated value
** \param[in]  vmin      Lower threshold
** \param[in]  vmax      Upper threshold
**
*****************************************************************************/
void AGibbs::_printInequalities(int iact,
                                int ivar,
                                int nfois,
                                int flag_cv,
                                double simval,
                                double vmin,
                                double vmax) const
{
  int flag_min,flag_max,idim;

  /* Initializations */

  int nech = _db->getSampleNumber();
  int iech = getSampleRank(iact);
  int nvar = getNvar();
  flag_min = flag_max = 1;
  if (FFFF(vmin)) flag_min = 0;
  if (FFFF(vmax)) flag_max = 0;

  /* Print the simulated value */

  message("Sample (%3d/%3d) - Variable (%3d/%3d) = %8.4lf in ",
          iech+1,nech,ivar+1,nvar,simval);

  /* Print the bounds */

  if (! flag_min)
    message("[      NA,");
  else
    message("[%8.4lf,",vmin);
  if (! flag_max)
    message("      NA]");
  else
    message("%8.4lf]",vmax);

  /* Print the coordinates */

  message(" at point (");
  for (idim=0; idim<_db->getNDim(); idim++)
  {
    if (idim != 0) message(",");
    message("%8.4lf",_db->getCoordinate(iech,idim));
  }
  message(")");

  message("\n");
}

/****************************************************************************/
/*!
 **  Print the initial status for Gibbs iterations
 **
 ** \param[in]  flag_init   TRUE for the Initial printout
 ** \param[in]  y           Gaussian vector
 ** \param[in]  isimu       Rank of the simulation
 ** \param[in]  ipgs        Rank of the GS
 **
 *****************************************************************************/
void AGibbs::print(bool flag_init,
                   const VectorVectorDouble& y,
                   int isimu,
                   int ipgs) const
{
  int nactive = _db->getActiveSampleNumber();
  int nvar    = getNvar();


  if (flag_init)
  {
    mestitle(1, "Gibbs Initial Status (Simu:%d - GS=%d)", isimu + 1, ipgs + 1);
  }
  else
  {
    mestitle(1, "Gibbs Results (Simu:%d - GS=%d)", isimu + 1, ipgs + 1);
    message("Number of bootstrap iterations = %d\n", _nburn);
    message("Total number of iterations     = %d\n", _niter);
  }

  /* Loop on the variables */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    mestitle(2,"Variable %d",ivar+1);
    int icase = getRank(ipgs, ivar);

    /* Loop on the samples */

    for (int iact = 0; iact < nactive; iact++)
    {
      int iech = getSampleRank(iact);
      double vmin = _db->getLowerBound(iech, icase);
      double vmax = _db->getUpperBound(iech, icase);
      _printInequalities(iact,ivar,-1,0,y[icase][iact], vmin, vmax);
    }
  }
}

int AGibbs::getDimension() const
{
  int nsize = _npgs * _nvar;
  return nsize;
}

int AGibbs::getRank(int ipgs, int ivar) const
{
  int rank = ivar + _nvar * (ipgs);
  return rank;
}

VectorVectorDouble AGibbs::allocY() const
{
  VectorVectorDouble y;

  int nsize = getDimension();
  int nactive = _db->getActiveSampleNumber();
  y.resize(nsize);
  for (int i = 0; i < nsize; i++)
    y[i].resize(nactive);
  return y;
}

/**
 * Store the Gaussian array in LOC_GAUS variable.
 * This should be performed once for all GS and all variables
 *
 * @param y The Gaussian vector to be stored
 * @param isimu Rank of the simulation
 * @param ipgs  Rank of the GS
 */
void AGibbs::storeResult(const VectorVectorDouble& y,
                         int isimu,
                         int ipgs)
{
  int nsize = getDimension();
  int nactive = _db->getActiveSampleNumber();
  int nvar    = getNvar();

  /* Loop on the variables */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int icase = getRank(ipgs, ivar);
    int rank  = icase + nsize * isimu;

      /* Loop on the samples */

    for (int iact = 0; iact < nactive; iact++)
    {
      int iech = getSampleRank(iact);
      _db->setFromLocator(LOC_GAUSFAC, iech,  rank,  y[icase][iact]);
    }
  }

  // In case of Statistics, process this information
  if (_flagStats) _stats.plot(isimu);
}

/****************************************************************************/
/*!
**  Initializes the Gibbs sampler for a set of inequalities
**
** \return  Error return code
**
** \param[in]  y             Gaussian vector
** \param[in]  isimu         Rank of the simulation
** \param[in]  ipgs          Rank of the GS
** \param[in]  verbose       Verbose flag
**
*****************************************************************************/
int AGibbs::calculInitialize(VectorVectorDouble& y,
                             int isimu,
                             int ipgs,
                             bool verbose)
{
  Db* db = getDb();
  Model* model = getModel();
  int nactive = db->getActiveSampleNumber();
  int nvar    = getNvar();

  /* Print the title */

  if (debug_query("converge"))
    mestitle(1,"Initial Values for Gibbs Sampler (Simu:%d - GS:%d)",
             isimu+1,ipgs+1);

  /* Loop on the variables */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int icase   = getRank(ipgs,ivar);

    /* Loop on the samples */

    double sk = sqrt(model->getTotalSill(ivar,ivar));
    for (int iact = 0; iact < nactive; iact++)
    {
      int iech = getSampleRank(iact);
      double vmin, vmax;
      if (_boundsCheck(iech, ipgs, ivar, &vmin, &vmax)) return 1;

      /* Compute the median value of the interval */

      double pmin = (FFFF(vmin)) ? 0. : law_cdf_gaussian(vmin);
      double pmax = (FFFF(vmax)) ? 1. : law_cdf_gaussian(vmax);
      y[icase][iact] = sk * law_invcdf_gaussian((pmin + pmax) / 2.);
    }
  }

  // Re-initialize the statistics (optional)

  if (_flagStats) _stats.clear();

  return(0);
}

/**
 * Generate a simulated value
 * @param y     : Gaussian vector
 * @param yk    : Kriged value
 * @param sk    : Standard deviation
 * @param iact  : Rank of the target sample (relative)
 * @param ipgs  : Rank of the current GS
 * @param ivar  : Rank of the current Variable
 * @param iter  : Rank of the iteration
 * @return Simulated value
 */
double AGibbs::getSimulate(VectorVectorDouble& y,
                           double yk,
                           double sk,
                           int iact,
                           int ipgs,
                           int ivar,
                           int iter)
{
  // Define the environment

  int icase = getRank(ipgs, ivar);
  int iech = getSampleRank(iact);

  // Read the Bounds

  double vmin = _db->getLowerBound(iech, icase);
  double vmax = _db->getUpperBound(iech, icase);

  // Apply decay

  if (_nburn > 0 && _flagDecay && iter <= _nburn)
  {
    double ratio = (double) iter / (double) _nburn;
    if (!FFFF(vmin))
      vmin = THRESH_INF + (vmin - THRESH_INF) * ratio;
    if (!FFFF(vmax))
      vmax = THRESH_SUP + (vmax - THRESH_SUP) * ratio;
  }

  // In multi-mono case, correct from the previously (linked) variable

  double yval;
  double sval;
  if (_flagMultiMono && ivar > 0)
  {
    int icase0 = getRank(ipgs, 0);
    double sqr = getSqr();
    yval = yk * sqr + getRho() * y[icase0][iact];
    sval = sk * sqr;
  }
  else
  {
    yval = yk;
    sval = sk;
  }

  /* Update the definition interval */

  if (!FFFF(vmin)) vmin = (vmin - yval) / sval;
  if (!FFFF(vmax)) vmax = (vmax - yval) / sval;

  /* Draw an authorized normal value */

  return (yk + sk * law_gaussian_between_bounds(vmin, vmax));
}

VectorInt AGibbs::calculateSampleRanks() const
{
  VectorInt ranks;
  if (! _db->hasSelection()) return ranks;
  for (int iech = 0; iech < _db->getSampleNumber(); iech++)
  {
    if (_db->isActive(iech)) ranks.push_back(iech);
  }
  return ranks;
}

int AGibbs::getSampleRankNumber() const
{
  if (_ranks.empty())
    return _db->getSampleNumber();
  else
    return _ranks.size();
}

int AGibbs::getSampleRank(int i) const
{
  if (_ranks.empty())
    return i;
  else
    return _ranks[i];
}

void AGibbs::updateStats(const VectorVectorDouble& y, int ipgs, int niter)
{
  if (! _flagStats) return;

  // Calculate the number of columns
  int ncols = getDimension();

  // Resize the statistics array
  _stats.resize(niter+1, ncols);

  // Loop on the columns

  for (int ivar = 0; ivar < getNvar(); ivar++)
  {
    int icol = getRank(ipgs, ivar);

    // Update statistics
    double value = _stats.getValue(niter-1, icol);
    value = (value * niter + ut_vector_mean(y[icol])) / (niter+1.);

    // Update the statistics array
    _stats.update(niter, icol, value);
  }
}

