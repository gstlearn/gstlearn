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
#include "geoslib_old_f.h"

#include "Neigh/ANeigh.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighBench.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"

#include <math.h>
#include <algorithm>

ANeigh::ANeigh(const ASpace* space)
  : ASpaceObject(space)
  , ASerializable()
  , _dbin(nullptr)
  , _dbout(nullptr)
  , _dbgrid(nullptr)
  , _rankColCok()
  , _iechMemo(-1)
  , _flagSimu(false)
  , _flagXvalid(false)
  , _flagKFold(false)
  , _useBallSearch(false)
  , _ballLeafSize(10)
  , _flagIsUnchanged(false)
  , _nbghMemo()
  , _ball()
{
}

ANeigh::ANeigh(const ANeigh& r)
  : ASpaceObject(r)
  , ASerializable(r)
  , _dbin(r._dbin)
  , _dbout(r._dbout)
  , _dbgrid(r._dbgrid)
  , _rankColCok(r._rankColCok)
  , _iechMemo(r._iechMemo)
  , _flagSimu(r._flagSimu)
  , _flagXvalid(r._flagXvalid)
  , _flagKFold(r._flagKFold)
  , _useBallSearch(r._useBallSearch)
  , _ballLeafSize(r._ballLeafSize)
  , _flagIsUnchanged(r._flagIsUnchanged)
  , _nbghMemo(r._nbghMemo)
{
}

ANeigh& ANeigh::operator=(const ANeigh &r)
{
  if (this != &r)
  {
    ASpaceObject::operator=(r);
    ASerializable::operator=(r);
    _dbin = r._dbin;
    _dbout = r._dbout;
    _dbgrid = r._dbgrid;
    _rankColCok = r._rankColCok;
    _iechMemo = r._iechMemo;
    _flagSimu = r._flagSimu;
    _flagXvalid = r._flagXvalid;
    _flagKFold  = r._flagKFold;
    _useBallSearch = r._useBallSearch;
    _ballLeafSize = r._ballLeafSize;
    _flagIsUnchanged = r._flagIsUnchanged;
    _nbghMemo = r._nbghMemo;
  }
  return *this;
}

ANeigh::~ANeigh()
{
}

int ANeigh::attach(const Db *dbin, const Db *dbout)
{
  if (dbin == nullptr || dbout == nullptr) return 1;
  _dbin  = dbin;
  _dbout = dbout;

  // Check if the output Db is defined and is a grid
  if (_dbout != nullptr) _dbgrid = dynamic_cast<const DbGrid*>(_dbout);

  // Attach the Ball Tree search (if relevant)

  attachBall();

  setIsChanged();
  return 0;
}

void ANeigh::attachBall(double (*dist_function)(const double* x1,
                                                const double* x2,
                                                int size))
{
  // Attach the Ball only if the option is switched ON
  if (!_useBallSearch) return;

  // Nothing can be done unless the Input Db is specifiied
  if (_dbin == nullptr) return;

  _ball.init(_dbin, dist_function, _ballLeafSize);
}

void ANeigh::setIsChanged(bool status)
{
  _flagIsUnchanged = status;
  _nbghMemo.clear();
};

void ANeigh::reset()
{
  _flagIsUnchanged = false;
  _nbghMemo.clear();
  _rankColCok.clear();
  _iechMemo   = -1;
  _flagSimu   = false;
  _flagXvalid = false;
  _flagKFold  = false;
}

/**
 * Generic function for performing neighborhood selection.
 * This function ALWAYS modifies (and resizes) the returned array 'ranks'
 * @param iech_out Rank of the target point (in 'dbout')
 * @param ranks Input / Output vector of neighboring sample ranks
 */
void ANeigh::select(int iech_out, VectorInt& ranks)
{
  // Validating the call with respect to the input and output Dbs
  if (_dbin == nullptr || _dbout == nullptr)
  {
    messageAbort("'dbin' and 'dbout' must have been attached beforehand");
    return;
  }
  if (! _dbout->isSampleIndexValid(iech_out))
  {
    ranks.clear();
    return;
  }

  // Check if the current target coincides with the previous one
  // Then do not do anything (even in presence of colocation)
  if (_isSameTarget(iech_out))
  {
    ranks = _nbghMemo;
    _flagIsUnchanged = true;
    return;
  }

  // Performing the selection of the neighboring samples

  if (hasChanged(iech_out))
  {
    getNeigh(iech_out, ranks);

    // Set the flag telling if neighborhood has changed or not
    // and memorize the new set of ranks
    _checkUnchanged(iech_out, ranks);
  }
  else
  {
    ranks = _nbghMemo;
    _flagIsUnchanged = true;
  }

  // Stop the neighborhood search if not enough point is available
  if ((int) ranks.size() <= 0) return;

  // Update in case of Colocated option

  _updateColCok(ranks, iech_out);
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
}

void ANeigh::_neighCompress(VectorInt& ranks)
{
  int necr = 0;
  int number = (int) ranks.size();
  for (int i = 0; i < number; i++)
    if (ranks[i] >= 0) ranks[necr++] = i;
  ranks.resize(necr);
}

/****************************************************************************/
/*!
 **  Print the information selected in the neighborhood
 **
 ** \param[in]  ranks     Array of the data ranks
 ** \li                   -1 if not selected
 ** \li                   >=0 gives the angular sector in ENeigh::MOVING
 **
 *****************************************************************************/
void ANeigh::_display(const VectorInt& ranks)
{
  String string;
  int ndim = _dbin->getNDim();
  int nech = _dbin->getSampleNumber();
  bool flag_ext = _dbin->getLocNumber(ELoc::BLEX) > 0;

  /* Neighborhood data */

  mestitle(1, "Data selected in neighborhood");
  tab_prints(NULL, "Rank");
  tab_prints(NULL, "Sample");
  if (_dbin->hasLocVariable(ELoc::C)) tab_prints(NULL, "Code");
  for (int idim = 0; idim < ndim; idim++)
  {
    string = getLocatorName(ELoc::X, idim);
    tab_prints(NULL, string.c_str());
  }
  if (flag_ext)
  {
    for (int idim = 0; idim < ndim; idim++)
    {
      string = getLocatorName(ELoc::BLEX, idim);
      tab_prints(NULL, string.c_str());
    }
  }
  if (getType() == ENeigh::MOVING) tab_prints(NULL, "Sector");
  message("\n");

  /* Loop on the sample points */

  int nsel = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    if (ranks[iech] < 0) continue;

    tab_printi(NULL, nsel + 1);
    tab_printi(NULL, iech + 1);
    if (_dbin->hasLocVariable(ELoc::C))
      tab_printi(NULL, static_cast<int>(_dbin->getLocVariable(ELoc::C,iech,0)));
    for (int idim = 0; idim < ndim; idim++)
      tab_printg(NULL, _dbin->getCoordinate(iech, idim));
    if (flag_ext)
    {
      for (int idim = 0; idim < ndim; idim++)
        tab_printg(NULL, _dbin->getLocVariable(ELoc::BLEX,iech, idim));
    }
    if (getType() == ENeigh::MOVING) tab_printi(NULL, ranks[iech] + 1);
    message("\n");
    nsel++;
  }
}

/****************************************************************************/
/*!
 **  Discard a sample for which all variables are undefined
 **
 **  Returns 1 if all variables are undefined; 0 otherwise
 **
 ** \param[in]  iech      Rank of the sample
 **
 ** \remarks When the Neighborhood is performed in the case of Simulations
 ** \remarks checking for all variables being undefined is performed
 ** \remarks on ELoc::SIMU rather than on ELoc::Z
 **
 *****************************************************************************/
bool ANeigh::_discardUndefined(int iech)
{
  if (_dbin->getLocNumber(ELoc::Z) <= 0) return 0;

  if (! _flagSimu)
  {
    if (_dbin->isAllUndefined(iech)) return 0;
  }
  else
  {
    // Here the check is performed on all variables belonging to ELoc::SIMU
    // The number of variables is equal to 'nbsimu' * 'nvar'
    if (_dbin->isAllUndefinedByType(ELoc::SIMU, iech)) return 0;
  }
  return 1;
}

/****************************************************************************/
/*!
 **  Mask the data sample in the case of cross-validation
 **
 ** \return  1 if the sample is masked; 0 otherwise
 **
 ** \param[in]  iech_in  Rank in the input Db structure
 ** \param[in]  iech_out Rank in the output Db structure
 ** \param[in]  eps      Tolerance
 **
 *****************************************************************************/
int ANeigh::_xvalid(int iech_in, int iech_out, double eps)
{
  if (! getFlagXvalid()) return 0;
  if (! getFlagKFold())
  {
    if (distance_inter(_dbin, _dbout, iech_in, iech_out, NULL) < eps) return 1;
  }
  else
  {
    if (! _dbin->hasLocVariable(ELoc::C)) return 0;
    if (_dbin->getLocVariable(ELoc::C,iech_in,0) ==
        _dbout->getLocVariable(ELoc::C,iech_out,0)) return 1;
  }
  return 0;
}

bool ANeigh::_isDimensionValid(int idim) const
{
  if (idim < 0 || idim >= (int) getNDim())
  {
    messerr("Error in 'idim'(%d). It should lie within [0,%d[",idim,getNDim());
    return false;
  }
  return true;
}

bool ANeigh::_deserialize(std::istream& is, bool /*verbose*/)
{
  int ndim = 0;

  bool ret = true;
  ret = ret && _recordRead<int>(is, "Space Dimension", ndim);
  if (ret) setNDim(ndim);
  return ret;
}

bool ANeigh::_serialize(std::ostream& os, bool /*verbose*/) const
{
  bool ret = true;
  ret = ret && _recordWrite<int>(os, "Space Dimension", getNDim());

  return ret;
}

void ANeigh::setBallSearch(bool status, int leaf_size)
{
  _useBallSearch = status;
  _ballLeafSize = leaf_size;
}