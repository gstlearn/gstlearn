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

#include "Basic/Timer.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Law.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "geoslib_old_f.h"
#include "geoslib_define.h"
#include "geoslib_enum.h"

#include <math.h>

AGibbs::AGibbs()
    : AStringable(),
      _npgs(1),
      _nvar(1),
      _nburn(1),
      _niter(1),
      _flagOrder(true),
      _flagDecay(true),
      _optionStats(0),
      _ranks(),
      _db(nullptr),
      _stats()
{
}

AGibbs::AGibbs(Db* db)
    : AStringable(),
      _npgs(1),
      _nvar(1),
      _nburn(1),
      _niter(1),
      _flagOrder(true),
      _flagDecay(true),
      _optionStats(0),
      _ranks(),
      _db(db),
      _stats()
{
}

AGibbs::AGibbs(Db* db,
               int npgs, int nvar, int nburn, int niter, int seed,
               int flag_order, bool flag_decay)
    : AStringable(),
      _npgs(1),
      _nvar(1),
      _nburn(1),
      _niter(1),
      _flagOrder(false),
      _flagDecay(true),
      _optionStats(0),
      _ranks(),
      _db(db),
      _stats()
{
  init(npgs, nvar, nburn, niter, seed, flag_order, flag_decay);
}

AGibbs::AGibbs(const AGibbs &r)
    : AStringable(r),
      _npgs(r._npgs),
      _nvar(r._nvar),
      _nburn(r._nburn),
      _niter(r._niter),
      _flagOrder(r._flagOrder),
      _flagDecay(r._flagDecay),
      _optionStats(r._optionStats),
      _ranks(r._ranks),
      _db(r._db),
      _stats(r._stats)
{
}

AGibbs& AGibbs::operator=(const AGibbs &r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _npgs = r._npgs;
    _nvar = r._nvar;
    _nburn = r._nburn;
    _niter = r._niter;
    _flagOrder = r._flagOrder;
    _flagDecay = r._flagDecay;
    _optionStats = r._optionStats;
    _ranks = r._ranks;
    _db = r._db;
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
                  int seed,
                  int flag_order,
                  bool flag_decay)
{
  _npgs = npgs;
  _nvar = nvar;
  _nburn = nburn;
  _niter = niter;
  _flagOrder = flag_order;
  _flagDecay = flag_decay;

  // Evaluate the array of active sample ranks
  _ranks = calculateSampleRanks();

  // Initialize the series of Random Numbers
  law_set_random_seed(seed);
}

/****************************************************************************/
/*!
**  Correct the bounds according to the order relationship
**
** \return  Error code: 1 if there is no solution; 0 otherwise
**
** \param[in]  ipgs        Rank of the GS
** \param[in]  ivar        Rank of the variable
** \param[in]  iact        Internal Rank of the sample
**
** \param[out]  vmin_arg   Output minimum bound
** \param[out]  vmax_arg   Output maximum bound
**
** \remark Attributes ELoc::GAUSFAC are mandatory
**
*****************************************************************************/
int AGibbs::_boundsCheck(int ipgs,
                         int ivar,
                         int iact,
                         double *vmin_arg,
                         double *vmax_arg)
{
  const Db* db = getDb();
  int icase = getRank(ipgs, ivar);
  int iech  = getSampleRank(iact);
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
** \param[in]  simval    Simulated value
** \param[in]  vmin      Lower threshold
** \param[in]  vmax      Upper threshold
**
*****************************************************************************/
void AGibbs::_printInequalities(int iact,
                                int ivar,
                                double simval,
                                double vmin,
                                double vmax) const
{
  int flag_min,flag_max,idim;

  /* Initializations */

  int iech = getSampleRank(iact);
  int nvar = getNvar();
  flag_min = flag_max = 1;
  if (FFFF(vmin)) flag_min = 0;
  if (FFFF(vmax)) flag_max = 0;

  /* Print the simulated value */

  message("Sample (%3d/%3d) - Variable (%3d/%3d) = %8.4lf in ",
          iech+1,_db->getSampleNumber(),ivar+1,nvar,simval);

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
void AGibbs::_displayCurrentVector(bool flag_init,
                                   const VectorVectorDouble& y,
                                   int isimu,
                                   int ipgs) const
{
  int nact = getSampleRankNumber();
  int nvar = getNvar();

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

    for (int iact = 0; iact < nact; iact++)
    {
      int iech = getSampleRank(iact);
      double vmin = _db->getLowerBound(iech, icase);
      double vmax = _db->getUpperBound(iech, icase);
      _printInequalities(iact,ivar,y[icase][iact], vmin, vmax);
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
  int nact  = getSampleRankNumber();
  y.resize(nsize);
  for (int i = 0; i < nsize; i++)
    y[i].resize(nact);
  return y;
}

/**
 * Store the Gaussian array in ELoc::GAUS variable.
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
  int nact  = getSampleRankNumber();
  int nvar  = getNvar();

  /* Loop on the variables */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int icase = getRank(ipgs, ivar);
    int rank  = icase + nsize * isimu;

      /* Loop on the samples */

    for (int iact = 0; iact < nact; iact++)
    {
      int iech = getSampleRank(iact);
      _db->setFromLocator(ELoc::GAUSFAC, iech,  rank,  y[icase][iact]);
    }
  }

  // In case of Statistics, process this information
  if (_optionStats == 0)
    return;
  else if (_optionStats == 1)
    _stats.display();
  else if (_optionStats == 2)
    _stats.plot(isimu);
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
    return static_cast<int>(_ranks.size());
}

int AGibbs::getSampleRank(int i) const
{
  if (_ranks.empty())
    return i;
  else
    return _ranks[i];
}

int AGibbs::getSampleNumber() const
{
  if (_db == nullptr)
    return 0;
  else
    return _db->getSampleNumber();
}

void AGibbs::updateStats(const VectorVectorDouble& y,
                         int ipgs,
                         int jter,
                         double amort)
{
  if (_optionStats == 0) return;
  if (jter < _nburn) return;
  int iter = jter - _nburn;

  // Loop on the columns

  for (int ivar = 0; ivar < getNvar(); ivar++)
  {

    // Update statistics

    int jcol;
    double result;
    int icol = getRank(ipgs, ivar);
    double residu = 1. - amort;
    double oldw   = (1. - pow(amort, (double) iter))   / residu;
    double neww   = (1. - pow(amort, (double) iter+1)) / residu;

    // The mean
    jcol = _getColRankStats(ipgs, ivar, 0);
    double oldmean = _stats.getValue(iter-1, jcol);
    double newmean = VH::mean(y[icol]);
    result = (oldmean * oldw * amort + newmean) / neww;
    _stats.update(iter, jcol, result);

    // The standard deviation
    jcol = _getColRankStats(ipgs, ivar, 1);
    double oldvar = _stats.getValue(iter-1,jcol);
    double newvar = VH::variance(y[icol]);
    result = (oldvar * oldw * amort + newvar) / neww;
    _stats.update(iter, jcol, result);
  }
}

/**
 * Returns the number of Rows for storing the statistics
 * This number is based on the number of iterations, exclusing the burnout
 * @return
 */
int AGibbs::_getRowNumberStats() const
{
  int nrows = _niter - _nburn;
  return nrows;
}

/**
 * Returns the number of Columns for storing the statistics
 * This number is based on:
 * - the number of GS
 * - the number of variables (or GRF)
 * - the storage of mean and standard deviation
 * @return
 */
int AGibbs::_getColNumberStats() const
{
  int ncols = 2 * getDimension();
  return ncols;
}

/**
 * Return the column for the storage of a value in Statistics
 * @param ipgs  Rank of the GS
 * @param ivar  Rank of the Variable or GRF
 * @param mode  0 for Mean and 1 for Standard Deviation
 * @return
 */
int AGibbs::_getColRankStats(int ipgs, int ivar, int mode) const
{
  int rank = getRank(ipgs, ivar);
  if (mode == 0)
    return (2 * rank);
  else
    return (2 * rank + 1);
}

/**
 * Test wheter a constraint is tight at a sample (data is a hard data)
 * @param ipgs  Rank of the GS
 * @param ivar  Rank of the variable
 * @param iact  Rank of the sample (within internal _ranks)
 * @param value Constraining value (if sample is an active constraint)
 * @return
 */
bool AGibbs::isConstraintTight(int ipgs,
                               int ivar,
                               int iact,
                               double* value) const
{
  int icase = getRank(ipgs, ivar);
  int iech = getSampleRank(iact);

  double vmin = _db->getLowerBound(iech, icase);
  double vmax = _db->getUpperBound(iech, icase);

  bool isActive = !FFFF(vmin) && !FFFF(vmax) && vmin == vmax;
  if (isActive)
    *value = vmin;
  else
    *value = TEST;
  return isActive;
}

void AGibbs::statsInit()
{
  if (_optionStats == 0) return;
  _stats.init(_getRowNumberStats(), _getColNumberStats());
}

/**
 * Given the initial bounds, calculate the modified bounds
 * during the Burning stage
 * @param iter: Rank of the iteration
 * @param vmin: Lower bound (in/out)
 * @param vmax: Upper bound (in/out)
 */
void AGibbs::getBoundsDecay(int iter, double *vmin, double *vmax) const
{
  // Do not modify the bounds if no Decay is defined
  if (! _flagDecay) return;
  // Do not modify the bounds after the burning stage
  if (iter > _nburn) return;

  double ratio = (double) iter / (double) _nburn;
  if (!FFFF(*vmin))
    *vmin = THRESH_INF + ((*vmin) - THRESH_INF) * ratio;
  if (!FFFF(*vmax))
    *vmax = THRESH_SUP + ((*vmax) - THRESH_SUP) * ratio;
}

/**
 * Find the relative rank of an absolute sample rank
 * @param iech Absolute sample rank
 * @return The rank within the vector if relative ranks (or -1)
 */
int AGibbs::getRelativeRank(int iech)
{
  int nact = getSampleRankNumber();
  for (int iact = 0; iact < nact; iact++)
  {
    if (getSampleRank(iact) == iech) return iact;
  }
  return -1;
}

int AGibbs::run(VectorVectorDouble& y, int ipgs, int isimu, bool verbose, bool flagCheck)
{
  if (calculInitialize(y, isimu, ipgs)) return 1;
  if (flagCheck)
    _displayCurrentVector(true, y, isimu, ipgs);

  /* Iterations of the Gibbs sampler */

  Timer timer;
  for (int iter = 0; iter < getNiter(); iter++)
    update(y, isimu, ipgs, iter);
  if (verbose) timer.displayIntervalMilliseconds("Gibbs iterations");

  /* Check the validity of the Gibbs results (optional) */

  if (flagCheck)
  {
    checkGibbs(y, isimu, ipgs);
    _displayCurrentVector(false, y, isimu, ipgs);
  }

  // Store the results

  storeResult(y, isimu, ipgs);

  // Clean the spurious elements

  cleanup();

  return 0;
}

String AGibbs::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << toTitle(0, "Gibbs Characteristics");

  sstr << "Number of Gaussian Systems" << _npgs;
  sstr << "Number of Variables" << _nvar;
  sstr << "Number of Gibbs Iterations" << _niter;
  sstr << "Number of Burning Iterations" << _nburn;
  if (_flagDecay)
    sstr << "Decay option is switched ON" << std::endl;
  if (_flagOrder == 1)
    sstr << "Variables are ordered sequentially upwards" << std::endl;
  if (_flagOrder == -1)
    sstr << "Variables are ordered sequentially downwards" << std::endl;
  if (_optionStats == 1)
    sstr << "Statistics on Trajectories are stored for print out" << std::endl;
  if (_optionStats == 2)
    sstr << "Statistics on Trajectories are stored in Neutral File" << std::endl;

  return sstr.str();

}
