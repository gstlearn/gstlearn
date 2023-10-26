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
#include "geoslib_define.h"

#include "Enum/EPostUpscale.hpp"

#include "Db/Db.hpp"
#include "Basic/Grid.hpp"
#include "Matrix/Table.hpp"
#include "Calculators/CalcSimuPost.hpp"

CalcSimuPost::CalcSimuPost()
    : ACalcDbToDb(),
      _verbose(false),
      _flagMatch(false),
      _flagUpscale(false),
      _checkLevel(0),
      _checkTargets(),
      _upscale(EPostUpscale::UNKNOWN),
      _stats(),
      _names(),
      _iechout(0),
      _iter(-1),
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

  if (_flagUpscale)
  {
    if (! hasDbout())  return false;
    if (! isGridOut()) return false;
  }
  else
  {
    // Set the Dbout equal to Dbin
    setDbout(getDbin());
  }

  /**************************************************/
  /* Cross-checking the Space Dimension consistency */
  /**************************************************/

  if (getDbin()->getNDim() > getDbout()->getNDim())
  {
    messerr("The Space Dimension of Dbin(%d) must not be larger than the one of Dbout(%d)",
            getDbin()->getNDim(), getDbout()->getNDim());
    return false;
  }

  // Cross-checking options
  if (_flagUpscale)
  {
    if (_upscale == EPostUpscale::UNKNOWN)
    {
      messerr("When 'dbout' is specified, some Upscaling is required");
      messerr("Therefor the 'upscale' option must be defined");
      return false;
    }
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
    message(" (using the 'matching' criterion)\n");
  else
    message(" (using the 'product' criterion)\n");

  message("Number of Statistics: %d\n", _getNStats());
  if (_getTransfoNvar() > 0)
    message("Number of Transform Variables: %d\n", _getTransfoNvar());
  message("Number of Output Variables: %d\n", _getNVarout());
}

bool CalcSimuPost::_mustBeChecked(int level) const
{
  if (_checkTargets.empty()) return false;
  if (level > _checkLevel) return false;
  return (VH::isInList(_checkTargets, _iechout+1));
}

bool CalcSimuPost::_preprocess()
{
  if (_flagUpscale)
    _iattOut = _addVariableDb(2, 1, ELoc::UNKNOWN, 0, _getNVarout(), 0.);
  else
    _iattOut = _addVariableDb(1, 1, ELoc::UNKNOWN, 0, _getNVarout(), 0.);
  if (_iattOut < 0) return false;

  return true;
}

bool CalcSimuPost::_postprocess()
{
  int ecr = 0;
  for (int ivar = 0, nvar = _getNEff(); ivar < nvar; ivar++)
    for (int istat = 0, nstat = _getNStats(); istat < nstat; istat++)
    {
      std::ostringstream oper;
      oper << "Var" << ivar + 1 << "." << _stats[istat].getDescr();
      if (_flagUpscale)
        _renameVariable(2, VectorString(), ELoc::UNKNOWN, 0, _iattOut + ecr, oper.str(), 1);
      else
        _renameVariable(1, VectorString(), ELoc::UNKNOWN, 0, _iattOut + ecr, oper.str(), 1);
      ecr++;
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
  if (_getNStats() <= 0)
  {
    messerr("The argument 'stats' should not be left empty");
    return 1;
  }

  // Loop on the statistic options

  int nvarin = _getNEff();
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
 * @param indices Vector of variable indices (Dimension: 'n')
 * @param tabin Output Vector of variables (Dimension: 'n')
 */
void CalcSimuPost::_readIn(int iech, const VectorInt& indices, VectorDouble& tabin) const
{
  int nvar = _getNVar();
  for (int ivar = 0; ivar < nvar; ivar++)
    tabin[ivar] = getDbin()->getArray(iech, _iuids[ivar][indices[ivar]]);

  if (_mustBeChecked(1))
  {
    message("    Sample Rank #%d - Coordinates:",iech);
    for (int idim = 0, ndim = getDbin()->getNDim(); idim < ndim; idim++)
      message(" %lf", getDbin()->getCoordinate(iech,  idim));
    message("\n");
  }
  if (_mustBeChecked(2))
    VH::display("    Initial    ", tabin, false);
}

void CalcSimuPost::_writeOut(int iech, const VectorDouble& tabout) const
{
  for (int ivar = 0; ivar < _getNVarout(); ivar++)
    getDbout()->setArray(iech, _iattOut + ivar, tabout[ivar]);
}

/**
 * Upscale the multivariate array defined for several samples into a multivariate array (in place)
 * @param Y_p_k_s Information for several samples (Dimension: 's' x 'd') where 'd' is 'n' or 'p'
 * @param tabout Vector of upscaled values (Dimension: 'd')
 */
void CalcSimuPost::_upscaleFunction(const VectorVectorDouble& Y_p_k_s, VectorDouble& tabout) const
{
  int nsample = (int) Y_p_k_s.size();
  int nvar = (int) Y_p_k_s[0].size();

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
      tabout[ivar] = TEST;
    else
    {
      if (_upscale == EPostUpscale::MEAN)
        tabout[ivar] = result / ndef;
      else
        tabout[ivar] = result;
    }
  }

  if (_mustBeChecked(2))
  {
    std::ostringstream string;
    string << "    Upscaled (" << nsample << ")";
    VH::display(string.str(), tabout, false);
  }
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
  int nvarout = _getNEff();
  int nstat   = _getNStats();

  int ecr = 0;
  for (int ivar = 0; ivar < nvarout; ivar++)
  {
    VectorDouble local(niter, TEST);
    if (!Y_p.empty())
    {
      int ndef = 0;
      for (int jter = 0; jter < niter; jter++)
      {
        double value = Y_p[jter][ivar];
        if (FFFF(value)) continue;
        local[ndef++] = value;
      }
    }

    for (int istat = 0; istat < nstat; istat++)
    {
      if (_stats[istat] == EPostStat::MEAN)
        tabout[ecr++] = VH::mean(local);
      else if (_stats[istat] == EPostStat::VAR)
        tabout[ecr++] = VH::variance(local);
      else if (_stats[istat] == EPostStat::VARP)
        tabout[ecr++] = VH::variance(local, true);
      else if (_stats[istat] == EPostStat::STD)
        tabout[ecr++] = VH::stdv(local);
      else if (_stats[istat] == EPostStat::STDP)
        tabout[ecr++] = VH::stdv(local, true);
      else if (_stats[istat] == EPostStat::MINI)
        tabout[ecr++] = VH::minimum(local);
      else if (_stats[istat] == EPostStat::MAXI)
        tabout[ecr++] = VH::maximum(local);
      else if (_stats[istat] == EPostStat::MED)
        tabout[ecr++] = VH::median(local);
    }
  }

  if (_mustBeChecked(0) && ! Y_p.empty())
  {
    Table stat_table = Table(_getNEff(),_getNStats());
    stat_table.setTitle("Statistics");
    stat_table.setSkipTitle(true);
    stat_table.setSkipDescription(true);

    for (int istat = 0; istat < nstat; istat++)
      stat_table.setColumnName(istat, _stats[istat].getDescr());

    for (int ivar = 0; ivar < _getNEff(); ivar++)
    {
      std::ostringstream name;
      name << "Var " << ivar+1;
      stat_table.setRowName(ivar, name.str());
    }

    int lec = 0;
    for (int ivar = 0; ivar < _getNEff(); ivar++)
      for (int istat = 0; istat < nstat; istat++)
        stat_table.setValue(ivar, istat, tabout[lec++]);
    stat_table.display();
  }
}

void CalcSimuPost::_printIndices(const VectorVectorInt &indices) const
{
  int nvar = _getNVar();
  message("  Iteration (1-based) %3d/%3d -> Indices:", _iter+1, _niter);
  for (int ivar = 0; ivar < nvar; ivar++)
    message(" %d/%d", indices[_iter][ivar] + 1, _nfact[ivar]);
  message("\n");
}

VectorVectorInt CalcSimuPost::_getIndices() const
{
  int nvar = _getNVar();
  int niter = _getNiter();
  VectorVectorInt indices(niter);

  // Loop on the iterations

  for (int jter = 0; jter < niter; jter++)
  {
    indices[jter].resize(nvar,0);

    // Dispatch according to the multiplicity mode
    if (_flagMatch)
    {
      for (int ivar = 0; ivar < nvar; ivar++)
      {
        indices[jter][ivar] = jter;
      }
    }
    else
    {
      int local = jter;
      for (int ivar = 0; ivar < nvar; ivar++)
      {
        int jvar = nvar - ivar - 1;
        int divid = local / _nfact[jvar];
        indices[jter][jvar] = local - divid * _nfact[jvar];
        local = divid;
      }
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

/**
 * Returns the number of output variables, which is equal to
 * - the number of variables after the transformation step (if defined)
 * - otherwise the number of input variables
 * @return
 */
int CalcSimuPost::_getNEff() const
{
  if (_getTransfoNvar() > 0)
    return _getTransfoNvar();
  else
    return _getNVar();
}

VectorInt CalcSimuPost::_samplesInCellIdenticalSpaceDimension(const VectorInt& indblock) const
{
  VectorInt local;
  for (int iechin = 0, nechin = getDbin()->getSampleNumber(); iechin < nechin; iechin++)
  {
    if (!getDbin()->isActive(iechin)) continue;
    if (indblock[iechin] != _iechout) continue;
    local.push_back(iechin);
  }
  return local;
}

VectorInt CalcSimuPost::_samplesInCellDifferentSpaceDimension() const
{
  VectorInt local;
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
  for (int iechin = 0, nechin = getDbin()->getSampleNumber(); iechin < nechin; iechin++)
  {
    if (!getDbin()->isActive(iechin)) continue;
    VectorDouble coor = getDbin()->getSampleCoordinates(iechin);
    if (dbgrid->sampleBelongsToCell(coor, _iechout))
      local.push_back(iechin);
  }
  return local;
}

/**
 * Indicates which type of pre-sorting of the sample information must be chosen
 * It can be:
 * - 0: if one Target item corresponds to One Data item (no Upscaling phase)
 * - 1: when each Datum item belongs to one Target item at most
 *      (case when 'dbin' and 'dbout' have same space dimension)
 * - 2: when a Datum item may belong to several Target items
 *      (case when 'dbin' space dimension is smaller than 'dbout' space dimension)
 * @return
 */
int CalcSimuPost::_getSortingCase() const
{
  if (! _flagUpscale)
    return 0;
  else if (getDbin()->getNDim() == getDbout()->getNDim())
    return 1;
  else
    return 2;
}

int CalcSimuPost::_process()
{
  int nechout = getDbout()->getSampleNumber();
  int niter = _getNiter();
  VectorDouble sampleIn(_getNVar());
  VectorDouble sampleOut;
  if (_getTransfoNvar() > 0)
    sampleOut.resize(_getTransfoNvar());
  VectorDouble tabUpscaled;
  if (_flagUpscale)
    tabUpscaled.resize(_getNEff());
  VectorDouble statres(_getNVarout());

  // Get the indices of the samples within the Grid
  // There is no need to check that 'dbout' is a grid (see _check)
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
  VectorInt indblock;
  int icase = _getSortingCase();
  if (icase == 1)
    indblock = dbgrid->locateDataInGrid(getDbin(), VectorInt(), true);

  // Getting the indices of the multivariate simulation
  VectorVectorInt indices = _getIndices();

  // Loop on the samples of the Output File
  for (_iechout = 0; _iechout < nechout; _iechout++)
  {
    if (! getDbout()->isActive(_iechout)) continue;
    VectorVectorDouble Y_p;

    // Get the vector of samples contained in the target cell
    VectorInt local;
    if (icase == 0)
      local.push_back(_iechout);
    else if (icase == 1)
      local = _samplesInCellIdenticalSpaceDimension(indblock);
    else
      local = _samplesInCellDifferentSpaceDimension();
    if (local.empty()) continue;

    // Modif DR

    if (_mustBeChecked(0))
    {
      if (_flagUpscale)
        message("\n== Cell #%d/%d (regrouping %d samples)\n", _iechout+1, nechout, (int) local.size());
      else
        message("\n== Cell #%d/%d\n", _iechout+1, nechout);
    }

    // Loop on the iterations
    for (_iter = 0; _iter < niter; _iter++)
    {
      if (_mustBeChecked(2))
         _printIndices(indices);

      // Loop on the samples contained in the target cell
      VectorVectorDouble Z_n_k_s;
      for (int is = 0, nlocal = (int) local.size(); is < nlocal; is++)
      {
        int iechin = local[is];

        // Reading the variables for the current input sample rank and current iteration
        _readIn(iechin, indices[_iter], sampleIn);

        // Transformation (optional)
        if (_getTransfoNvar() <= 0)
          Z_n_k_s.push_back(sampleIn);
        else
        {
          _transformFunction(sampleIn, sampleOut);
          if (_mustBeChecked(2))
            VH::display("    Transformed",sampleOut, false);
          Z_n_k_s.push_back(sampleOut);
        }
      }

      if (_flagUpscale)
      {
        // Upscale
        _upscaleFunction(Z_n_k_s, tabUpscaled);
        Y_p.push_back(tabUpscaled);
      }
      else
      {
        // Simply add this iteratio to the stack (flatten VVD to VD)
        Y_p.push_back(VH::flatten(Z_n_k_s));
      }
    }

    // Calculate the statistics on the resulting multivariate vector

    _statisticsFunction(Y_p, statres);

    // Save the results

    _writeOut(_iechout, statres);
  }
  return 0;
}

/**
 * @param dbin Input data base
 * @param dbout Output data base (must be a Grid)
 * @param names Vector of simulation names
 * @param flag_match True if the ranks of simulations must match; otherwise: product
 * @param upscale Option within EPostUpscale
 * @param stats Vector of options within EPostStat
 * @param verbose Verbose flag
 * @param check_targets Rank (1-based) of the target element to be checked (0: None; -1: All)
 * @param check_level 0: Statistics; 1: Sample Selection; 2: Iteration definition
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
             const VectorInt& check_targets,
             int check_level,
             const NamingConvention& namconv)
{
  CalcSimuPost calcul;
  calcul.setDbin(dbin);
  if (dbout != nullptr)
  {
    calcul.setFlagUpscale(true);
    calcul.setDbout(dbout);
  }
  calcul.setNames(names);
  calcul.setUpscale(upscale);
  calcul.setStats(stats);
  calcul.setFlagMatch(flag_match);
  calcul.setVerbose(verbose);
  calcul.setCheckTargets(check_targets);
  calcul.setCheckLevel(check_level);
  calcul.setNamingConvention(namconv);

  int error = (calcul.run()) ? 0 : 1;
  return error;
}
