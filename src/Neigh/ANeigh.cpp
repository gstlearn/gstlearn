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
#include "geoslib_old_f.h"

#include "Neigh/ANeigh.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighBench.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Faults/Faults.hpp"
#include "Basic/AException.hpp"
#include "Basic/OptDbg.hpp"

#include <math.h>
#include <algorithm>
#include <set>

ANeigh::ANeigh(const Db *dbin, const ANeighParam *neighparam, const Db* dbout)
    : _dbin(dbin),
      _dbout(dbout),
      _dbgrid(nullptr),
      _neighParam(neighparam),
      _rankColCok(),
      _iechMemo(-1),
      _flagSimu(false),
      _flagIsUnchanged(false),
      _ptOut(),
      _nbghMemo()
{
  initialize(dbin, neighparam, dbout);
}

ANeigh::ANeigh(const ANeigh &r)
    : _dbin(r._dbin),
      _dbout(r._dbout),
      _dbgrid(r._dbgrid),
      _neighParam(r._neighParam),
      _rankColCok(r._rankColCok),
      _iechMemo(r._iechMemo),
      _flagSimu(r._flagSimu),
      _flagIsUnchanged(r._flagIsUnchanged),
      _ptOut(r._ptOut),
      _nbghMemo(r._nbghMemo)
{
}

ANeigh& ANeigh::operator=(const ANeigh &r)
{
  if (this != &r)
  {
    _dbin = r._dbin;
    _dbout = r._dbout;
    _dbgrid = r._dbgrid;
    _neighParam = r._neighParam;
    _rankColCok = r._rankColCok;
    _iechMemo = r._iechMemo;
    _flagSimu = r._flagSimu;
    _flagIsUnchanged = r._flagIsUnchanged;
    _ptOut = r._ptOut;
    _nbghMemo = r._nbghMemo;
  }
  return *this;
}

ANeigh::~ANeigh()
{
}

int ANeigh::initialize(const Db *dbin,
                       const ANeighParam *neighparam,
                       const Db *dbout)
{
  if (neighparam == nullptr || dbin == nullptr) return 1;
  _neighParam = neighparam;
  _dbin  = dbin;
  _dbout = dbout;

  // Check if the output Db is defined and is a grid
  if (_dbout != nullptr)
    _dbgrid = dynamic_cast<const DbGrid*>(_dbout);
  return 0;
}

void ANeigh::setIsChanged()
{
  _flagIsUnchanged = false;
  _nbghMemo.clear();
};

VectorInt ANeigh::select(int iech_out)
{
  if (_dbout == nullptr) return VectorInt();
  if (! _dbout->isSampleIndexValid(iech_out)) return VectorInt();
  VectorInt ranks;

  // Check if the current target coincides with the previous one
  // Then do not do anything (even in presence of colocation)
  if (_isSameTarget(iech_out)) return _nbghMemo;

  // Select the neighboring samples

  if (hasChanged(iech_out))
    ranks = getNeigh(iech_out);
  else
    ranks = _nbghMemo;

  // Stop the neighborhood search if not enough point is available
  if ((int) ranks.size() <= 0) return ranks;

  // Set the flag telling if neighborhood has changed or not
  // and memorize the new set of ranks
  _checkUnchanged(iech_out, ranks);

  // Update in case of Colocated option

  _updateColCok(ranks, iech_out);
  return ranks;
}

/**
 * Checks if the current target matches the target previously treated
 * in the same procedure. If match is reached, then there is no need
 * to compute a new neighborhood: use the previous Vector of sample ranks.
 * Store the references of the new 'dbout' and 'iech_out' for next optimizations
 * @param iech_out Rank of the current target sample
 * @return
 */
bool ANeigh::_isSameTarget(int iech_out)
{
  // Check if the target remained unchanged
  bool flagSame = true;
  if (_iechMemo < 0 || iech_out != _iechMemo) flagSame = false;
  if (_nbghMemo.empty()) flagSame = false;

  _flagIsUnchanged = flagSame;
  return flagSame;
}

void ANeigh::_checkUnchanged(int iech_out, const VectorInt &ranks)
{
  VectorInt rsorted = ranks;
  std::sort(rsorted.begin(), rsorted.end());

  // Check if Neighborhood has changed

  if (_nbghMemo.size() != ranks.size())
    // Two series do not share the same dimension

    _flagIsUnchanged = false;
  else
  {
    // Two (sorted) series have the same size: check if they are equal

    _flagIsUnchanged = (rsorted == _nbghMemo);
  }

  // Store the vector of sample ranks for the current neighborhood search

  _iechMemo = iech_out;
  _nbghMemo = rsorted;
}

/**
 * Update the set of selected samples in case of colocated option
 * This is done only if:
 * - the colocation option is ON (vector of colocated variable is defined)
 * - at least one of the colocated variables at the target is valid
 * - the target does not coincide with a sample already selected
 * If the colocation option is validated, an additional member is added to 'ranks':
 * its value is conventionally set to -1.
 * @param ranks      Vector of samples already selected
 * @param iech_out   Rank of the target site (in dbout)
 */
void ANeigh::_updateColCok(VectorInt &ranks, int iech_out)
{
  if (_rankColCok.empty()) return;
  int nvarin = (int) _rankColCok.size();

  /* Do not add the target if no variable is defined */
  bool found = false;
  for (int ivar = 0; ivar < nvarin && !found; ivar++)
  {
    int jvar = _rankColCok[ivar];
    if (jvar < 0) continue;
    if (!FFFF(_dbout->getArray(iech_out, jvar))) found = true;
  }
  if (! found) return;

  /* Do not add the target if it coincides with an already selected sample */
  int nsel = (int) ranks.size();
  for (int iech = 0; iech < nsel; iech++)
  {
    if (distance_inter(_dbin, _dbout, ranks[iech], iech_out, NULL) <= 0.)
      return;
  }

  /* Add the target */

  ranks.push_back(-1);
  _flagIsUnchanged = false;
  return;
}
