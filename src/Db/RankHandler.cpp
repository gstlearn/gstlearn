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
#include "Db/RankHandler.hpp"

#include "Db/Db.hpp"
#include "Basic/VectorHelper.hpp"

RankHandler::RankHandler(const Db* db,
                         bool useSel,
                         bool useZ,
                         bool useVerr,
                         bool useExtD)
  : _useSel(useSel)
  , _useZ(useZ)
  , _useVerr(useVerr)
  , _useExtD(useExtD)
  , _nvar(0)
  , _nExtD(0)
  , _iptrSel(-1)
  , _iptrZ()
  , _iptrVerr()
  , _iptrExtD()
  , _nbgh()
  , _index()
  , _Zflatten()
  , _db(db)
  , _workNbgh()
{
  _nvar = _db->getNLoc(ELoc::Z);
  if (_nvar <= 0) _nvar = 1;

  // Prepare the vector of sample ranks per variable
  _index.resize(_nvar);

  // Column index for selection (if 'useSel' and if present)
  _iptrSel = (_useSel) ? _db->getColIdxByLocator(ELoc::SEL, 0) : -1;

  // Column indices for variables (if 'useZ' and if present)
  _iptrZ = VectorInt();
  if (useZ && _db->hasLocator(ELoc::Z))
  {
    _iptrZ.resize(_nvar);
    for (int ivar = 0; ivar < _nvar; ivar++)
      _iptrZ[ivar] = _db->getColIdxByLocator(ELoc::Z, ivar);
  }

  // Column indices for variance of measurement error (if 'useVerr' and if present)
  _iptrVerr = VectorInt();
  if (_useVerr && _db->hasLocator(ELoc::V))
  {
    _iptrVerr.resize(_nvar);
    for (int ivar = 0; ivar < _nvar; ivar++)
      _iptrVerr[ivar] = _db->getColIdxByLocator(ELoc::V, ivar);
  }

  // Column indices for external drifts (if 'useExtD' and if present)
  _nExtD = 0;
  _iptrExtD = VectorInt();
  if (_useExtD && _db->hasLocator(ELoc::F))
  {
    _nExtD = _db->getNLoc(ELoc::F);
    _iptrExtD.resize(_nExtD);
    for (int iext = 0; iext < _nExtD; iext++)
      _iptrExtD[iext] = _db->getColIdxByLocator(ELoc::F, iext);
  }
}

RankHandler::RankHandler(const RankHandler& r)
  : _useSel(r._useSel)
  , _useZ(r._useZ)
  , _useVerr(r._useVerr)
  , _useExtD(r._useExtD)
  , _nvar(r._nvar)
  , _nExtD(r._nExtD)
  , _iptrSel(r._iptrSel)
  , _iptrZ(r._iptrZ)
  , _iptrVerr(r._iptrVerr)
  , _iptrExtD(r._iptrExtD)
  , _nbgh(r._nbgh)
  , _index(r._index)
  , _Zflatten(r._Zflatten)
  , _db(r._db)
  , _workNbgh(r._workNbgh)
{
}

RankHandler& RankHandler::operator=(const RankHandler& r)
{
  if (this != &r)
  {
    _useSel   = r._useSel;
    _useZ     = r._useZ;
    _useVerr  = r._useVerr;
    _useExtD  = r._useExtD;
    _nvar     = r._nvar;
    _nExtD    = r._nExtD;
    _iptrSel  = r._iptrSel;
    _iptrZ    = r._iptrZ;
    _iptrVerr = r._iptrVerr;
    _iptrExtD = r._iptrExtD;
    _nbgh     = r._nbgh;
    _index    = r._index;
    _Zflatten = r._Zflatten;
    _db       = r._db;
    _workNbgh = r._workNbgh;
  }
  return *this;
}

RankHandler::~RankHandler()
{
}

/**
 * Defines the list of indices 'index' for valid samples for the set of
 * variables (Z locator)
 *
 * @param nbgh Vector giving the ranks of the elligible samples (optional)
 */
void RankHandler::defineSampleRanks(const VectorInt& nbgh)
{
  // Define the ranks of the samples
  if (nbgh.empty())
  {
    int nech_init = _db->getNSample();
    VH::sequenceInPlace(nech_init, _workNbgh);
    _nbgh = _workNbgh;
  }
  else
  {
    _nbgh = nbgh;
  }
  int nech = (int)_nbgh.size();

  double value;
  _Zflatten.clear();

  // Loop on the variables
  for (int ivar = 0; ivar < _nvar; ivar++)
  {
    VectorInt ranks;

    // Loop on the elligible sample ranks
    for (int irel = 0; irel < nech; irel++)
    {
      int iabs = _nbgh[irel];

      // Check against a possible selection
      if (_iptrSel >= 0)
      {
        value = _db->getValueByColIdx(iabs, _iptrSel);
        if (value <= 0) continue;
      }
 
      // Check against validity of the Variance of Measurement Error variable
      if (!_iptrVerr.empty())
      {
        value = _db->getValueByColIdx(iabs, _iptrVerr[ivar]);
        if (FFFF(value) || value < 0) continue;
      }

      // Check against the validity of ALL external drifts
      if (!_iptrExtD.empty())
      {
        bool valid = true;
        for (int iext = 0; iext < _nExtD && valid; iext++)
        {
          value = _db->getValueByColIdx(iabs, _iptrExtD[iext]);
          if (FFFF(value)) valid = false;
        }
        if (!valid) continue;
      }

      // Check against the existence of a target variable
      if (!_iptrZ.empty())
      {
        value = _db->getValueByColIdx(iabs, _iptrZ[ivar]);
        if (FFFF(value)) continue;

        _Zflatten.push_back(value);
      }

      // The sample is finally accepted: its ABSOLUTE index is stored
      ranks.push_back(iabs);
    }

    // The vector of ranks is stored for the current variable
    _index[ivar] = ranks;
  }
}

/**
 * @brief Get the number of active samples for a given variable
 * 
 * @param ivar Rank of the target variable
 * @return int 
 */
int RankHandler::getCount(int ivar) const
{
  if (ivar < 0 || ivar >= _nvar)
  {
    messerr("RankHandler::getCount: invalid variable index %d", ivar);
    return -1;
  }
  return (int)_index[ivar].size();
}

/**
 * @brief Get the total number of active samples for all variables
 * 
 * @return int 
 */
int RankHandler::getTotalCount() const
{
  int count = 0;
  for (int ivar = 0; ivar < _nvar; ivar++)
    count += getCount(ivar);
  return count;
}

/**
 * @brief Return the total number of samples 
 * 
 * @return int 
 */
int RankHandler::getNumber() const
{
  return (int) _nbgh.size();
}

int RankHandler::identifyVariableRank(int ipos) const
{
  int ntotal = 0;
  for (int ivar = 0; ivar < _nvar; ivar++)
  {
    int size = getCount(ivar);
    if (ipos < ntotal + size)
      return ivar; // Found the variable
    ntotal += size;
  }
  return -1; // Not found
}

int RankHandler::identifySampleRank(int ipos) const
{
  int ntotal = 0;
  for (int ivar = 0; ivar < _nvar; ivar++)
  {
    int size = getCount(ivar);
    if (ipos < ntotal + size)
    {
      int jpos = ipos - ntotal;
      return _index[ivar][jpos]; // Found the sample
    }
    ntotal += size;
  }
  return -1; // Not found
}

void RankHandler::dump(bool flagFull) const
{
  mestitle(0, "Rank Handler");
  message("Use Selection: %d\n", (int)_useSel);
  message("Use Z-variable: %d\n", (int)_useZ);
  message("Use Variance of Measurement Error: %d\n", (int)_useVerr);
  message("Use External Drift: %d\n", (int)_useExtD);
  message("\n");
  message("Number of Variables: %d\n", _nvar);
  message("Number of External Drifts: %d\n", _nExtD);

  if (! flagFull) return;

  mestitle(1,"Variable and Sample ranks - Variable values");
  for (int ivar = 0, lec = 0; ivar < _nvar; ivar++)
  {
    message("Variable= %d: \n", ivar);
    int size = getCount(ivar);
    for (int i = 0; i < size; i++, lec++)
      message("- Sample= %2d : Variable value= %lf\n", _index[ivar][i], _Zflatten[lec]);
    message("\n");
  }
}