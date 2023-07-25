/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_define.h"

#include "Enum/EPostUpscale.hpp"

#include "Db/Db.hpp"
#include "Basic/Grid.hpp"
#include "Calculators/CalcSimuPost.hpp"

CalcSimuPost::CalcSimuPost()
    : ACalcDbToDb(),
      _verbose(false),
      _flagMatch(false),
      _rankCheck(0),
      _upscale(EPostUpscale::UNKNOWN),
      _stats(),
      _names(),
      _iechout(0),
      _iattOut(0),
      _niter(0),
      _nvarOut(0),
      _nfact(),
      _iuids()
{
}

CalcSimuPost::~CalcSimuPost()
{
}

bool CalcSimuPost::_check()
{
  /************************************************************/
  /* Both Files are compulsory: the output one must be a Grid */
  /************************************************************/

  if (! hasDbin())   return false;
  if (! hasDbout())  return false;
  if (! isGridOut()) return false;

  /**************************************************/
  /* Cross-checking the Space Dimension consistency */
  /**************************************************/

  if (getDbin()->getNDim() > getDbout()->getNDim())
  {
    messerr("The Space Dimension of Dbin(%d) must not be larger than the one of Dbout(%d)",
            getDbin()->getNDim(), getDbout()->getNDim());
    return false;
  }

  // Identify the variables from the input file
  if (_defineNames()) return false;

  // Define the number of iterations
  _defineIterations();

  // Define the final number of output variables
  if (_defineVaroutNumber()) return false;

  // Final printout (optional)
  _environPrint();

  return true;
}

void CalcSimuPost::_environPrint() const
{
  if (!_verbose) return;

  mestitle(1, "Simulation Post-Processing");

  // Count of iterations
  message("Multiplicity order for all variables\n");
  for (int ivar = 0, nvar = _getNVar(); ivar < nvar; ivar++)
    message("- Variable %d (%s) = %d\n", ivar + 1, _names[ivar].c_str(),
            _nfact[ivar]);
  message("Number of Iterations: %d", _getNiter());
  if (_flagMatch)
    message("(using the 'matching' criterion)\n");
  else
    message("(using the 'product' criterion)\n");

  message("Number of Statistics: %d\n", _getNStats());
  if (getTransfoNvar() > 0)
    message("Number of Transform Variables: %d\n", getTransfoNvar());
  message("Number of Output Variables: %d\n", _getNVarout());
}

bool CalcSimuPost::_mustBeChecked() const
{
  if (_rankCheck == 0) return false;
  if (_rankCheck < 0)  return true;
  return (_rankCheck == (_iechout+1));
}

bool CalcSimuPost::_preprocess()
{
  _iattOut = _addVariableDb(2, 1, ELoc::UNKNOWN, 0, _getNVarout(), 0.);
  if (_iattOut < 0) return false;

  return true;
}

bool CalcSimuPost::_postprocess()
{
  int ecr = 0;
  if (getTransfoNvar() > 0)
  {
    for (int ivar = 0, nvar = getTransfoNvar(); ivar < nvar; ivar++)
      for (int istat = 0, nstat = _getNStats(); istat < nstat; istat++)
      {
        std::ostringstream name;
        name << "Trans" << ivar+1 << "." << _stats[istat].getDescr();
        _renameVariable(2, 0, _iattOut + ecr, name.str(), 1);
        ecr++;
      }
  }
  else
  {
    for (int ivar = 0, nvar = _getNVar(); ivar < nvar; ivar++)
      for (int istat = 0, nstat = _getNStats(); istat < nstat; istat++)
      {
        std::ostringstream name;
        name << _names[ivar] << "." << _stats[istat].getDescr();
        _renameVariable(2, 0, _iattOut + ecr, name.str(), 1);
        ecr++;
      }
  }
  return true;
}

void CalcSimuPost::_rollback()
{
  _cleanVariableDb(1);
}

bool CalcSimuPost::_run()
{
  return _process() == 0;
}

int CalcSimuPost::_defineVaroutNumber()
{
  if (getTransfoNvar() <= 0)
  {
    messerr("The count of final variable is only possible after the number of Output variables is given");
    return 1;
  }
  if (_getNStats() <= 0)
  {
    messerr("The argument 'stats' should not be left empty");
    return 1;
  }

  // Loop on the statistic options

  int nvarin = getTransfoNvar();
  _nvarOut = 0;
  for (int ioption = 0, noption = _getNStats(); ioption < noption; ioption++)
  {
    if (_stats[ioption] != EPostStat::UNKNOWN)
      _nvarOut += nvarin;
  }

  if (_nvarOut <= 0)
  {
    messerr("The number of output variables cannot be zero");
    return 1;
  }
  return 0;
}

/**
 * Read a multivariate vector for a set of variable indices and a target sample
 * @param iech  Rank of the target sample
 * @param indices Vector of variable indices
 * @return
 */
VectorDouble CalcSimuPost::_readIn(int iech, const VectorInt& indices) const
{
  int nvar = _getNVar();
  VectorDouble tab(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
    tab[ivar] = getDbin()->getArray(iech, _iuids[ivar][indices[ivar]]);

  if (_mustBeChecked())
  {
    std::ostringstream phrase;
    phrase << "    Sample Rank #" << iech;
    VH::display(phrase.str(), tab, false);
  }
  return tab;
}

void CalcSimuPost::_writeOut(int iech, const VectorDouble& tabout) const
{
  for (int ivar = 0; ivar < _getNVarout(); ivar++)
    getDbout()->setArray(iech, _iattOut + ivar, tabout[ivar]);
}

/**
 * Upscale the multivariate array defined for several samples into a multivariate array (in place)
 * @param Y_p_k_s VectorVectorDouble containing information for several samples
 */
VectorDouble CalcSimuPost::_upscaleFunction(const VectorVectorDouble& Y_p_k_s) const
{
  int nsample = (int) Y_p_k_s.size();
  int nvar = (int) Y_p_k_s[0].size();
  VectorDouble Y_p_k(nvar);

  // Initialization values
  double valinit;
  if (_upscale == EPostUpscale::MINI)
    valinit = 1.e30;
  else if (_upscale == EPostUpscale::MAXI)
    valinit = -1.e30;
  else
    valinit = 0.;

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int ndef = 0;
    double result = valinit;
    for (int ip = 0; ip < nsample; ip++)
    {
      double value = Y_p_k_s[ip][ivar];
      if (FFFF(value)) continue;
      ndef++;

      // Perform the upscaling rule
      if (_upscale == EPostUpscale::NUM)
        result += 1;
      else if (_upscale == EPostUpscale::MEAN)
        result += value;
      else if (_upscale == EPostUpscale::MINI)
      {
        if (value < result) result = value;
      }
      else if (_upscale == EPostUpscale::MAXI)
      {
        if (value > result) result = value;
      }
      else
        messageAbort("Unknown Upscale mode");
    }

    if (ndef <= 0)
      Y_p_k[ivar] = TEST;
    else
    {
      if (_upscale == EPostUpscale::MEAN)
        Y_p_k[ivar] = result / ndef;
      else
        Y_p_k[ivar] = result;
    }
  }

  if (_mustBeChecked())
  {
    std::ostringstream phrase;
    phrase << "    Upscale (" << _upscale.getDescr() << ")";
    VH::display(phrase.str(), Y_p_k, false);
  }

  return Y_p_k;
}

/**
 * Calculate statistics on the various iterations contained in 'tabPileOut' and return the results
 * as a multivariate array 'tabout'
 * @param Y_p VectorVectorDouble containing information for several iterations
 * @param tabout Vector for multivariate statistics
 */
void CalcSimuPost::_statisticsFunction(const VectorVectorDouble &Y_p,
                                       VectorDouble &tabout) const
{
  int niter   = _getNiter();
  int nvarout = getTransfoNvar();
  int nstat   = _getNStats();

  int ecr = 0;
  for (int ivar = 0; ivar < nvarout; ivar++)
  {
    // Loop on the iterations
    int ndef = 0;
    double mm = 0.;
    double vv = 0.;
    double mini =  1.e30;
    double maxi = -1.e30;
    if (!Y_p.empty())
    {
      for (int iter = 0; iter < niter; iter++)
      {
        double value = Y_p[iter][ivar];
        if (FFFF(value)) continue;

        ndef++;
        mm += value;
        vv += value * value;
        if (value < mini) mini = value;
        if (value > maxi) maxi = value;
      }
    }

    if (ndef > 0)
    {
      mm = mm / ndef;
      vv = vv / ndef - mm * mm;
    }
    else
    {
      mm = vv = mini = maxi = TEST;
    }

    for (int istat = 0; istat < nstat; istat++)
    {
      if (_stats[istat] == EPostStat::MEAN)
        tabout[ecr++] = mm;
      else if (_stats[istat] == EPostStat::VAR)
        tabout[ecr++] = vv;
      else if (_stats[istat] == EPostStat::MINI)
        tabout[ecr++] = mini;
      else if (_stats[istat] == EPostStat::MAXI)
        tabout[ecr++] = maxi;
    }
  }

  if (_mustBeChecked() && ! Y_p.empty())
  {
    std::ostringstream phrase;
    phrase << "  Statistics (";
    for (int istat = 0, nstat = _getNStats(); istat < nstat; istat++)
    {
      if (istat > 0) phrase << ", ";
      phrase << _stats[istat].getDescr();
    }
    phrase << ")";
    VH::display(phrase.str(), tabout, false);
  }
}

void CalcSimuPost::_printIndices(int rank, const VectorInt &indices) const
{
  int nvar = _getNVar();
  message("  Iteration %3d/%3d -> Indices (1-based):", rank, _niter);
  for (int ivar = 0; ivar < nvar; ivar++)
    message(" %d/%d", indices[ivar] + 1, _nfact[ivar]);
  message("\n");
}

VectorInt CalcSimuPost::_getIndices(int rank) const
{
  if (rank < 0 || rank >= _niter)
  {
    messerr("Argument 'rank' cannot be negative or exceed _niter(%d)", _niter);
    return VectorInt();
  }

  int nvar = _getNVar();
  VectorInt indices(nvar);

  // Dispatch according to the multiplicity mode
  if (_flagMatch)
  {
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      indices[ivar] = rank;
    }
  }
  else
  {
    int local = rank;
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      int jvar = nvar - ivar - 1;
      int divid = local / _nfact[jvar];
      indices[jvar] = local - divid * _nfact[jvar];
      local = divid;
    }
  }
  return indices;
}

void CalcSimuPost::_defineIterations()
{
  // Switch according to the "match" flag

  if (_flagMatch)
  {
    _niter = VH::minimum(_nfact);
  }
  else
  {
    _niter = VH::product(_nfact);
  }
}

int CalcSimuPost::_defineNames()
{
  int nvar = _getNVar();
  if (getDbin() == nullptr)
  {
    messerr("The input Db must be defined beforehand");
    return 1;
  }
  if (nvar <= 0)
  {
    messerr("Some variables must be defined in the input Db");
    return 1;
  }

  // For each name, find the multiplicity nvar for each variable in the input Db

  _nfact.clear();
  _nfact.resize(nvar,0);
  _iuids.clear();
  _iuids.resize(nvar,0);
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    // Expand each filename
    VectorString subnames = getDbin()->expandNameList(_names[ivar]);

    // Get the multiplicity factor
    int nfois = (int) subnames.size();
    if (nfois <= 0)
    {
      messerr("The variable (%s) does not seem to exist in the input Db", _names[ivar].c_str());
      return 1;
    }
    _nfact[ivar] = nfois;

    // Identify the UID for each one of the expanded variable
    _iuids[ivar].resize(nfois);

    for (int ifois = 0; ifois < nfois; ifois++)
    {
      _iuids[ivar][ifois] = getDbin()->getUID(subnames[ifois]);
      if (_iuids[ivar][ifois] < 0)
      {
        messerr("The variable (%s) does not have a propoer UID", subnames[ifois].c_str());
        return 1;
      }
    }
  }
  return 0;
}

int CalcSimuPost::_process()
{
  int nechin = getDbin()->getSampleNumber();
  int nechout = getDbout()->getSampleNumber();
  int niter = _getNiter();
  VectorDouble statres(_getNVarout());

  // Get the indices of the samples within the Grid
  // There is no need to check that 'dbout' is a grid (see _check)
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
  VectorInt indblock = dbgrid->locateDataInGrid(getDbin(), VectorInt(), true);

  // Loop on the samples of the Output File
  for (_iechout = 0; _iechout < nechout; _iechout++)
  {
    if (! getDbout()->isActive(_iechout)) continue;
    VectorVectorDouble Y_p;

    // Loop on the iterations
    for (int iter = 0; iter < niter; iter++)
    {
      // Getting the indices of the multivariate simulation
      VectorInt indices = _getIndices(iter);

      // Loop on the samples contained in the target cell
      int ninside = 0;
      VectorVectorDouble Z_n_k_s;
      for (int iechin = 0; iechin < nechin; iechin++)
      {
        // We skip the masked samples and those which do not correspond to the target block
        if (!getDbin()->isActive(iechin)) continue;
        if (indblock[iechin] != _iechout) continue;

        if (_mustBeChecked())
        {
          if (ninside == 0)
          {
            if (iter == 0) message("\n== Cell #%d/%d\n", _iechout+1, nechout);
            _printIndices(iter, indices);
          }
        }

        // Reading the variables for the current input sample rank
        Z_n_k_s.push_back(_readIn(iechin, indices));
        ninside++;
      }
      if (ninside <= 0) continue;

      // Perform the transformation

      VectorVectorDouble Y_p_k_s;
      if (getTransfoNvar() <= 0)
      {
        Y_p_k_s = Z_n_k_s;
      }
      else
      {
        Y_p_k_s = transformFunction(Z_n_k_s);
        if (_mustBeChecked())
          VH::display("    After Transform  ",VH::flatten(Y_p_k_s), false);
      }

      // Upscale among all the samples of the target cell

      VectorDouble Y_p_k = _upscaleFunction(Y_p_k_s);

      // Stack the results for the current iteration

      Y_p.push_back(Y_p_k);
    }

    // Calculate the statistics on the intermediate vector

    _statisticsFunction(Y_p, statres);

    // Save the results

    _writeOut(_iechout, statres);
  }
  return 0;
}

/**
 * Returns the number of variables after transformation
 * In the current version, this number is set to 2
 * @return
 */
int CalcSimuPost::getTransfoNvar() const
{
  return 2;
}

/**
 * Perform the Transformation to convert the multivariate input vector 'tabin'
 * into a multivariate output vector.
 * @param Z_n_k_s  Input information (Dimension: 's' x 'n')
 *
 * @return Output information (Dimension: 's' x 'p')
 *
 * @remark In the current version, the Transformation consists in calculating:
 * - the sum of all values
 * - the standard deviation of all values
 */
VectorVectorDouble CalcSimuPost::transformFunction(const VectorVectorDouble& Z_n_k_s) const
{
  int nvarout = getTransfoNvar();
  int nsample = (int) Z_n_k_s.size();
  VectorVectorDouble Y_p_k_s(nsample);

  for (int is = 0; is < nsample; is++)
    Y_p_k_s[is].resize(nvarout);

  for (int is = 0; is < nsample; is++)
    Y_p_k_s[is][0] = VH::cumul(Z_n_k_s[is]);

  for (int is = 0; is < nsample; is++)
    Y_p_k_s[is][1] = VH::stdv(Z_n_k_s[is]);

  return Y_p_k_s;
}

/**
 *
 * @param dbin Input data base
 * @param dbout Output data base (must be a Grid)
 * @param names Vector of simulation names
 * @param flag_match True if the ranks of simulations must match; otherwise: product
 * @param upscale Option within EPostUpscale
 * @param stats Vector of options within EPostStat
 * @param verbose Verbose flag
 * @param rank_check Rank (1-based) of the target element to be checked (0: None; -1: All)
 * @param namconv Naming convention
 * @return Error code
 *
 * @remark
 * - N : number of variables in 'dbin' defined by 'names' (index 'n')
 * - k_1, ..., k_N : number of outcomes for each variable
 * - K : number of multivariate simulations according to 'flag_match' (index 'k')
 * - Transformation function 'F' from Z in R_N to Y in R_P (index 'p')
 *
 * Description of the Flow Chart:
 *
 * -# Loop over cells of 'dbout' (index 'C')
 * -# Loop over the simulations (index 'k')
 * -# Find the active samples of 'dbin' in 'C' (index 's') and build table of Z_n^{k}(s)
 * -# Apply the Transform function to table: Y_p^{k}(s)
 * -# Upscale to the target cell according to upscaling rule '_upscale': up_Y_p^{k}(C)
 * -# Compute statistics according to stat rule '_stats'
 */
int simuPost(Db *dbin,
             DbGrid *dbout,
             const VectorString &names,
             bool flag_match,
             const EPostUpscale &upscale,
             const std::vector<EPostStat>& stats,
             bool verbose,
             int rank_check,
             const NamingConvention& namconv)
{
  CalcSimuPost calcul;
  calcul.setDbin(dbin);
  calcul.setDbout(dbout);
  calcul.setNames(names);
  calcul.setUpscale(upscale);
  calcul.setStats(stats);
  calcul.setFlagMatch(flag_match);
  calcul.setVerbose(verbose);
  calcul.setRankCheck(rank_check);
  calcul.setNamingConvention(namconv);

  int error = (calcul.run()) ? 0 : 1;
  return error;
}
