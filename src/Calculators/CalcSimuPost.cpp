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

#include "Db/Db.hpp"
#include "Basic/Grid.hpp"
#include "Calculators/CalcSimuPost.hpp"

CalcSimuPost::CalcSimuPost()
    : ACalcDbToDb(),
      _verbose(false),
      _flagMatch(false),
      _niter(0),
      _names(),
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

  return true;
}

bool CalcSimuPost::_preprocess()
{
  return true;
}

bool CalcSimuPost::_postprocess()
{
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

/**
 * Operate the statistical operation to convert the multivariate input vector 'tabin'
 * into a multivariate output vector
 * @param tabin  Input vector
 *
 * @remark In the current version, the statistical operation consists in calculating:
 * - the sum of all values
 * - the product of all values
 * The resulting vector is dimensioned to 2.
 */
VectorDouble CalcSimuPost::_statisticsOperate(const VectorDouble &tabin) const
{
  int nvar = _getNVar();
  VectorDouble tabout(2);

  double sum = 0.;
  for (int ivar = 0; ivar < nvar; ivar++)
    sum += tabin[ivar];
  tabout[0] = sum;

  double prod = 1.;
  for (int ivar = 0; ivar < nvar; ivar++)
    prod *= tabin[ivar];
  tabout[1] = prod;

  if (_verbose) VH::display("", tabout);

  return tabout;
}

VectorDouble CalcSimuPost::_loadIn(int iech, const VectorInt& indices) const
{
  int nvar = _getNVar();
  VectorDouble tab(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
    tab[ivar] = getDbin()->getArray(iech, _iuids[ivar][indices[ivar]]);

  message("Sample #%d\n", iech);
  if (_verbose) VH::display("", tab);
  return tab;
}

/**
 * Merging the array defined for several samples into a multivariate array (in place)
 * @param tabPileIn VectorVectorDouble containing information for several samples
 * @param tabin Multivariate sample vector
 *
 * @remark: in the current version, we simply average the (defined) values per variable
 */
void CalcSimuPost::_mergePileIn(const VectorVectorDouble& tabPileIn, VectorDouble& tabin) const
{
  int nsample = (int) tabPileIn.size();
  int nvar = _getNVar();

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int ndef = 0;
    double result = 0.;
    for (int ip = 0; ip < nsample; ip++)
    {
      double value = tabPileIn[ip][ivar];
      if (FFFF(value)) continue;
      result += value;
      ndef++;
    }
    tabin[ivar] = (ndef > 0) ? result / ndef : TEST;
  }
}

/**
 * Merge the output arrays defined for several iterations into a multivariate array (in place)
 * @param tabPileOut VectorVectorDouble containing information for several iterations
 * @param tabout Multivariate output vector
 *
 * @remark: in the current version, we simply average the (defined) values per variable
 */
void CalcSimuPost::_mergePileOut(const VectorVectorDouble& tabPileOut, VectorDouble& tabout) const
{
  if (tabPileOut.empty()) return;
  int niter = (int) tabPileOut.size();
  int nvarout = (int) tabPileOut[0].size();

  for (int ivar = 0; ivar < nvarout; ivar++)
  {
    int ndef = 0;
    double result = 0;
    for (int iter = 0; iter < niter; iter++)
    {
      double value = tabPileOut[iter][ivar];
      if (FFFF(value)) continue;
      result += value;
      ndef++;
    }

    tabout[ivar] = (ndef > 0) ? result / ndef : TEST;
  }

  if (_verbose) VH::display("Average over the cell", tabout);
}

int CalcSimuPost::_process()
{
  int nvar = _getNVar();
  int nechin = getDbin()->getSampleNumber();
  int nechout = getDbout()->getSampleNumber();
  VectorDouble tabin(nvar);
  VectorDouble tabout(2); // This will be modified when the statistics Operator is defined

  // Get the indices of the samples within the Grid
  // There is no need to check that 'dbout' is a grid (see _check)
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
  VectorInt indblock = dbgrid->locateDataInGrid(getDbin());

  // Loop on the samples of the Output File
  for (int iechout = 0; iechout < nechout; iechout++)
  {

    // Loop on the combinations
    VectorVectorDouble tabPileOut;
    for (int iter = 0; iter < _niter; iter++)
    {
      // Getting the indices of the simulation outcome for each variable
      VectorInt indices = _getIndices(iter);

      // Loop on the samples contained in the target cell
      int ninside = 0;
      VectorVectorDouble tabPileIn;
      for (int iechin = 0; iechin < nechin; iechin++)
      {
        // We skip the masked samples and those which do not correspond to the target block
        if (!getDbin()->isActive(iechin)) continue;
        if (indblock[iechin] != iechout) continue;

        if (_verbose)
        {
          if (ninside == 0)
          {
            if (iter == 0) message("== Cell #%d/%d\n", iechout, nechout);
            _printIndices(iter, indices);
          }
        }

        // Loading the variables for the current input sample rank
        tabPileIn.push_back(_loadIn(iechin, indices));
        ninside++;
      }
      if (ninside <= 0) continue;

      // Merging the various samples contained in the same cell

      _mergePileIn(tabPileIn, tabin);

      // Operate the multivariate processing on the input vector

      tabPileOut.push_back(_statisticsOperate(tabin));
    }

    // Merging the various output vectors per cell

    _mergePileOut(tabPileOut, tabout);
  }
  return 0;
}

void CalcSimuPost::_printIndices(int rank, const VectorInt &indices) const
{
  if (!_verbose) return;
  int nvar = _getNVar();
  message("Iteration %3d/%3d -> Indices (1-based):", rank, _niter);
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

  // Optional verbose control
  if (_verbose)
  {
    int nvar = _getNVar();
    mestitle(1,"Simulation Post-Processing");
    message("Multiplicity order for all variables\n");
    for (int ivar = 0; ivar < nvar; ivar++)
      message("- Variable %d (%s) = %d\n", ivar+1, _names[ivar].c_str(), _nfact[ivar]);
    if (_flagMatch)
      message("Using the 'matching' criterion, ");
    else
      message("Using the 'product' criterion, ");
    message("the final number of outcomes is %d\n", _niter);
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

int simuPost(Db *dbin,
             DbGrid *dbout,
             const VectorString &names,
             bool flag_match,
             bool verbose,
             const NamingConvention &namconv)
{
  CalcSimuPost calcul;
  calcul.setDbin(dbin);
  calcul.setDbout(dbout);
  calcul.setNames(names);
  calcul.setFlagMatch(flag_match);
  calcul.setVerbose(verbose);
  calcul.setNamingConvention(namconv);

  int error = (calcul.run()) ? 0 : 1;
  return error;
}
