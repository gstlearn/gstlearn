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
#include "geoslib_old_f.h"

#include "Neigh/NeighWork.hpp"
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

NeighWork::NeighWork(const Db *dbin, const ANeighParam *neighparam)
    :
    _dbin(),
    _neighParam(),
    _flagInitialized(false),
    _flagIsUnchanged(false),
    _movingInd(),
    _movingIsect(),
    _movingNsect(),
    _movingX1(),
    _movingX2(),
    _movingDst(),
    _flagSimu(false),
    _dbout(nullptr),
    _iechOut(-1),
    _nbghMemo()
{
  initialize(dbin, neighparam);
}

NeighWork::NeighWork(const NeighWork &r)
    :
    _dbin(r._dbin),
    _neighParam(r._neighParam),
    _flagInitialized(r._flagInitialized),
    _flagIsUnchanged(r._flagIsUnchanged),
    _movingInd(r._movingInd),
    _movingIsect(r._movingIsect),
    _movingNsect(r._movingNsect),
    _movingX1(r._movingX1),
    _movingX2(r._movingX2),
    _movingDst(r._movingDst),
    _flagSimu(r._flagSimu),
    _dbout(r._dbout),
    _iechOut(r._iechOut),
    _nbghMemo(r._nbghMemo)
{
}

NeighWork& NeighWork::operator=(const NeighWork &r)
{
  if (this != &r)
  {
    _dbin = r._dbin;
    _neighParam = r._neighParam;
    _flagInitialized = r._flagInitialized;
    _flagIsUnchanged = r._flagIsUnchanged;
    _movingInd = r._movingInd;
    _movingIsect = r._movingIsect;
    _movingNsect = r._movingNsect;
    _movingX1 = r._movingX1;
    _movingX2 = r._movingX2;
    _movingDst = r._movingDst;
    _flagSimu = r._flagSimu;
    _dbout = r._dbout;
    _iechOut = r._iechOut;
    _nbghMemo = r._nbghMemo;
  }
  return *this;
}

NeighWork::~NeighWork()
{
}

/****************************************************************************/
/*!
 **  Initialize the neighborhood search
 **
 ** \param[in]  dbin          input Db structure
 ** \param[in]  neighparam    Description of the ANeighParam parameters
 **
 ** \remarks When the Neighborhood is performed in the case of Simulations
 ** \remarks checking for all variables being undefined is performed
 ** \remarks on ELoc::SIMU rather than on ELoc::Z
 **
 *****************************************************************************/
void NeighWork::initialize(const Db *dbin,
                           const ANeighParam *neighparam)
{
  if (neighparam == nullptr || dbin == nullptr) return;
  _neighParam = neighparam;
  _dbin = dbin;

  int nech  = _dbin->getSampleNumber();
  int ndim  = _dbin->getNDim();
  const NeighMoving* neighM = dynamic_cast<const NeighMoving*>(_neighParam);
  if (neighM != nullptr)
  {
    int nsect = neighM->getNSect();
    _movingInd = VectorInt(nech);
    _movingDst = VectorDouble(nech);
    _movingIsect = VectorInt(nsect);
    _movingNsect = VectorInt(nsect);
    _movingX1 = VectorDouble(ndim);
    _movingX2 = VectorDouble(ndim);
  }

  // Clear out the Memorization parameters

  _clearMemory();

  _flagInitialized = true;
}

/**
 * Clear all the local storage in the NeighWork structure
 */
void NeighWork::clear()
{
  /* Initialization */

  if (! _flagInitialized) return;

  // Clear the pointers

  _neighParam = nullptr;
  _dbin  = nullptr;

  /* Core deallocation */

  _movingInd.clear();
  _movingDst.clear();
  _movingIsect.clear();
  _movingNsect.clear();
  _movingX1.clear();
  _movingX2.clear();

  _nbghMemo.clear();

  _flagInitialized = false;
}

/****************************************************************************/
/*!
 **  Select the neighborhood
 **
 ** \return  Vector of sample ranks in neighborhood (empty when error)
 **
 ** \param[in]  dbout         output Db structure
 ** \param[in]  iech_out      Valid Rank of the sample in the output Db
 ** \param[in]  rankColCok    Vector of Colcok information (optional)
 ** \param[in]  verbose       Verbose option
 **
 *****************************************************************************/
VectorInt NeighWork::select(Db *dbout,
                            int iech_out,
                            const VectorInt& rankColCok,
                            bool verbose)
{
  if (! _flagInitialized) return VectorInt();
  if (! dbout->isSampleIndexValid(iech_out)) return VectorInt();
  int nech = _dbin->getSampleNumber();
  VectorInt ranks(nech, -1);

  // Optional title (only in verbose case)

  if (verbose)
    message(">>> Neighborhood search:\n");

  // Check if the current target coincides with the previous one
  // Then do not do anything (even in presence of colocation)
  if (_isSameTarget(dbout, iech_out, ranks, verbose)) return ranks;

  // Select the neighborhood samples as the target sample has changed
  bool doCompress = true;
  switch (_neighParam->getType().toEnum())
  {
    case ENeigh::E_IMAGE:
    case ENeigh::E_UNIQUE:
      if (! _isSameTargetUnique(dbout, iech_out, ranks, verbose))
        _unique(dbout, iech_out, ranks);
      else
        doCompress = false;
      break;

    case ENeigh::E_BENCH:
      if (! _isSameTargetBench(dbout, iech_out, ranks, verbose))
        _bench(dbout, iech_out, ranks);
      else
        doCompress = false;
      break;

    case ENeigh::E_MOVING:
      if (_moving(dbout, iech_out, ranks)) return VectorInt();
      break;
  }

  if (doCompress)
  {
    // In case of debug option, dump out neighborhood characteristics
    if (OptDbg::query(EDbg::NBGH)) _display(ranks);

    /* Compress the vector of returned sample ranks */
    int necr = 0;
    for (int iech = 0; iech < nech; iech++)
      if (ranks[iech] >= 0) ranks[necr++] = iech;
    ranks.resize(necr);

    // Set the flag telling if neighborhood has changed or not
    // and memorize the new set of ranks

    _checkUnchanged(dbout, iech_out, ranks);
  }

  // Update in case of Colocated option

  _updateColCok(rankColCok, ranks, iech_out);
  return ranks;
}

/****************************************************************************/
/*!
 **  Select the unique neighborhood (or Image Neighborhood)
 **
 ** \param[in]  dbout     output Db structure
 ** \param[in]  iech_out  rank of the output sample
 **
 ** \param[out]  ranks   Vector of samples elected in the Neighborhood
 **
 *****************************************************************************/
void NeighWork::_unique(Db *dbout, int iech_out, VectorInt& ranks)
{
  int nech = _dbin->getSampleNumber();

  /* Loop on samples */

  for (int iech = 0; iech < nech; iech++)
  {
    /* Discard the masked input sample */

    if (! _dbin->isActive(iech)) continue;

    /* Discard samples where all variables are undefined */

    if (_discardUndefined(iech)) continue;

    /* Discard the target sample for the cross-validation option */

    if (_neighParam->getFlagXvalid())
    {
      if (_xvalid(dbout, iech, iech_out)) continue;
    }
    ranks[iech] = 0;
  }
}

/****************************************************************************/
/*!
 **  Search for the bench neighborhood, according to the last
 **  coordinate
 **
 ** \param[in]  dbout     output Db structure
 ** \param[in]  iech_out  rank of the output sample
 **
 ** \param[out]  ranks    Vector of samples elected in the Neighborhood
 **
 *****************************************************************************/
void NeighWork::_bench(Db *dbout, int iech_out, VectorInt& ranks)
{
  int nech = _dbin->getSampleNumber();
  int idim_bench = _dbin->getNDim() - 1;
  double z0 = dbout->getCoordinate(iech_out, idim_bench);
  const NeighBench* neighB = dynamic_cast<const NeighBench*>(_neighParam);

  /* Loop on samples */

  for (int iech = 0; iech < nech; iech++)
  {
    /* Discard the masked input sample */

    if (! _dbin->isActive(iech)) continue;

    /* Discard samples where all variables are undefined */

    if (_discardUndefined(iech)) continue;

    /* Discard the target sample for the cross-validation option */

    if (_neighParam->getFlagXvalid())
    {
      if (_xvalid(dbout, iech, iech_out)) continue;
    }

    /* Discard sample located outside the bench */

    if (ABS(_dbin->getCoordinate(iech,idim_bench) - z0) <= neighB->getWidth())
      ranks[iech] = 0;
  }
}

/****************************************************************************/
/*!
 **  Moving neighborhood search
 **
 ** \return  Error return code
 ** \return  0 : No error
 ** \return  1 : The number of data is smaller than the minimum number of
 ** \return      data in Moving Neighborhood
 **
 ** \param[in]  dbout     Output Db structure
 ** \param[in]  iech_out  Rank of the target in the output Db structure
 ** \param[in]  eps       Tolerance
 **
 ** \param[out]  ranks    Vector of samples elected in the Neighborhood
 **
 *****************************************************************************/
int NeighWork::_moving(Db *dbout, int iech_out, VectorInt& ranks, double eps)
{
  const NeighMoving* neighM = dynamic_cast<const NeighMoving*>(_neighParam);
  int nech = _dbin->getSampleNumber();
  int isect = 0;
  if (nech < neighM->getNMini()) return 1;

  /* Loop on the data points */

  double distmax = 0.;
  int nsel = 0;

  for (int iech = 0; iech < nech; iech++)
  {

    /* Discard the masked input sample */

    if (! _dbin->isActive(iech)) continue;

    /* Discard samples where all variables are undefined */

    if (_discardUndefined(iech)) continue;

    /* Discard the target sample for the cross-validation option */

    if (_neighParam->getFlagXvalid())
    {
      if (_xvalid(dbout, iech, iech_out)) continue;
    }

    // In presence of Faults, check that sample 'iech' is still eligible

    if (neighM->hasFaults())
    {
      if (_hiddenByFault(dbout, iech, iech_out)) continue;
    }

    // Force selection if the sample belongs to the block

    double dist;
    if (neighM->getForceWithinBlock() && dbout->isGrid())
    {
      if (! _belongsToCell(dbout, iech, iech_out)) continue;
      dist = _movingDist(dbout, iech, iech_out);
    }
    else
    {
      /* Calculate the distance between data and target */

      dist = _movingDist(dbout, iech, iech_out);
      if (! FFFF(neighM->getRadius()) && dist > neighM->getRadius()) continue;
      if (dist > distmax) distmax = dist;

      /* Calculate the angular sector to which the sample belongs */

      if (neighM->getFlagSector())
        isect = _movingSectorDefine(_movingX1[0], _movingX1[1]);
    }

    /* The sample may be selected */

    _movingInd[nsel] = iech;
    _movingDst[nsel] = dist;
    ranks[iech] = isect;
    nsel++;
  }
  if (nsel < neighM->getNMini()) return 1;

  /* Slightly modify the distances in order to ensure the sorting results */
  /* In the case of equal distances                                       */

  for (int isel = 0; isel < nsel; isel++)
    _movingDst[isel] += distmax * isel * eps;

  /* Sort the selected samples according to the distance */

  ut_sort_double(0, nsel, _movingInd.data(), _movingDst.data());

  /* For each angular sector, select the first sample up to the maximum */

  if (! neighM->getForceWithinBlock() && neighM->getFlagSector() && neighM->getNSMax() > 0)
  {
    _movingSectorNsmax(nsel, ranks);
    if (nsel < neighM->getNMini()) return 1;
  }

    /* Select the first data samples (skipped if forcing all samples in block) */

  if (! neighM->getForceWithinBlock())
  {
    _movingSelect(nsel, ranks);
  }

  return 0;
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
bool NeighWork::_discardUndefined(int iech)
{
  if (_dbin->getVariableNumber() <= 0) return 0;

  if (! _flagSimu)
  {
    if (_dbin->isAllUndefined(iech)) return 0;
  }
  else
  {
    // In the case of simulations, the test is performed on the
    // simulation error for the first variable and first simulation
    if (! FFFF(_dbin->getSimvar(ELoc::SIMU, iech, 0, 0, 0, 1, 0))) return 0;
  }
  return 1;
}

/****************************************************************************/
/*!
 **  Mask the data sample in the case of cross-validation
 **
 ** \return  1 if the sample is masked; 0 otherwise
 **
 ** \param[in]  dbout    output Db structure
 ** \param[in]  iech_in  Rank in the input Db structure
 ** \param[in]  iech_out Rank in the output Db structure
 ** \param[in]  eps      Tolerance
 **
 *****************************************************************************/
int NeighWork::_xvalid(Db *dbout, int iech_in, int iech_out, double eps)
{
  if (! _neighParam->getFlagXvalid())
    return 0;
  else if (! _neighParam->getFlagKFold())
  {
    if (distance_inter(_dbin, dbout, iech_in, iech_out, NULL) < eps) return 1;
  }
  else
  {
    if (! _dbin->hasCode()) return 0;
    if (_dbin->getCode(iech_in) == dbout->getCode(iech_out)) return 1;
  }
  return 0;
}

bool NeighWork::_hiddenByFault(Db* dbout, int iech, int iech_out) const
{
  const NeighMoving* neighM = dynamic_cast<const NeighMoving*>(_neighParam);
  if (neighM == nullptr) return false;

  double xt1 = dbout->getCoordinate(iech_out, 0);
  double yt1 = dbout->getCoordinate(iech_out, 1);
  double xt2 = _dbin->getCoordinate(iech, 0);
  double yt2 = _dbin->getCoordinate(iech, 1);

  // Loop on the Fault polylines

  const Faults* faults = neighM->getFaults();
  if (faults == nullptr) return false;

  return faults->isSplitByFault(xt1, yt1, xt2, yt2);
}

/****************************************************************************/
/*!
 **  Calculates the distance between an input sample and the target
 **
 ** \return  Distance
 **
 ** \param[in]  dbout    output Db structure
 ** \param[in]  iech_in  Rank of the sample in the input Db structure
 ** \param[in]  iech_out Rank of the sample in the output Db structure
 **
 *****************************************************************************/
double NeighWork::_movingDist(Db *dbout, int iech_in, int iech_out)
{
  const NeighMoving* neighM = dynamic_cast<const NeighMoving*>(_neighParam);
  int ndim = _dbin->getNDim();

  /* Calculate the distance to the target */

  for (int idim = 0; idim < ndim; idim++)
    _movingX1[idim] = dbout->getCoordinate(iech_out, idim)
        - _dbin->getCoordinate(iech_in, idim);

  /* Anisotropic neighborhood */

  if (neighM->getFlagAniso())
  {

    /* Rotated anisotropy ellipsoid */

    if (neighM->getFlagRotation())
    {
      matrix_product(1, ndim, ndim, _movingX1.data(),
                     neighM->getAnisoRotMats().data(), _movingX2.data());
      _movingX1 = _movingX2;
    }
    for (int idim = 0; idim < ndim; idim++)
      _movingX1[idim] /= neighM->getAnisoCoeff(idim);
  }

  /* Calculate the distance */

  double dist;
  matrix_product(1, ndim, 1, _movingX1.data(), _movingX1.data(), &dist);
  dist = sqrt(dist);
  return dist;
}

bool NeighWork::_belongsToCell(Db* dbout, int iech, int iech_out)
{
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(dbout);
  if (dbgrid == nullptr) return false;

  // Get the coordinates of the sample
  VectorDouble coor = _dbin->getSampleCoordinates(iech);

  // Identify the dimensions of the cell
  VectorDouble dxsPerCell = dbgrid->getBlockExtensions(iech_out);

  // Check if the sample belongs to the cell
  return dbgrid->getGrid().sampleBelongsToCell(coor, iech_out, dxsPerCell);
}

/****************************************************************************/
/*!
 **  Returns the rank of the sector (using the first two coordinates)
 **
 ** \return  Sector index
 **
 ** \param[in]  dx    increment along X
 ** \param[in]  dy    increment along Y
 **
 *****************************************************************************/
int NeighWork::_movingSectorDefine(double dx, double dy)
{
  const NeighMoving* neighM = dynamic_cast<const NeighMoving*>(_neighParam);
  double angle;

  int isect = 0;
  if (neighM->getNSect() > 1)
  {
    if (dx == 0.)
    {
      if (dy >= 0.)
        angle = GV_PI / 2.;
      else
        angle = 1.5 * GV_PI;
    }
    else if (dx > 0.)
    {
      if (dy >= 0.)
        angle = atan(dy / dx);
      else
        angle = 2. * GV_PI - atan(-dy / dx);
    }
    else
    {
      if (dy > 0.)
        angle = GV_PI / 2. + atan(-dx / dy);
      else
        angle = GV_PI + atan(dy / dx);
    }
    isect = (int) (neighM->getNSect() * angle / (2. * GV_PI));
  }
  return (isect);
}

/****************************************************************************/
/*!
 **  For each angular sector, select the first samples until
 **  the maximum number of samples is reached
 **
 ** \param[in]  nsel   Number of selected samples
 **
 ** \param[out]  ranks Array of active data point ranks
 **
 *****************************************************************************/
void NeighWork::_movingSectorNsmax(int nsel, VectorInt& ranks)
{
  const NeighMoving* neighM = dynamic_cast<const NeighMoving*>(_neighParam);
  int n_end = 0;
  for (int isect = 0; isect < neighM->getNSect(); isect++)
  {
    int n_ang = 0;
    for (int i = 0; i < nsel; i++)
    {
      int j = _movingInd[i];
      if (ranks[j] != isect) continue;
      if (n_ang < neighM->getNSMax())
        n_ang++;
      else
        ranks[j] = -1;
    }
    n_end += n_ang;
  }
  return;
}

/****************************************************************************/
/*!
 **  Select the closest data per sector, until the maximum number
 **  neighbors is reached
 **
 ** \param[in]  nsel  Number of ellegible data points
 ** \param[in]  ranks Rank of the ellegible samples
 **
 ** \remark  The samples beyond the maximum number of neighbors have their
 ** \remark  rank turned into -1
 **
 *****************************************************************************/
void NeighWork::_movingSelect(int nsel, VectorInt& ranks)
{
  const NeighMoving* neighM = dynamic_cast<const NeighMoving*>(_neighParam);
  int number;

  if (neighM->getNMaxi() <= 0) return;
  for (int isect = 0; isect < neighM->getNSect(); isect++)
    _movingNsect[isect] = _movingIsect[isect] = 0;

  /* Count the number of samples per sector */

  number = 0;
  for (int i = 0; i < nsel; i++)
  {
    int j = _movingInd[i];
    int isect = ranks[j];
    if (isect < 0) continue;
    _movingNsect[isect]++;
    number++;
  }
  if (number < neighM->getNMaxi()) return;

  /* Find the rank of the admissible data per sector */

  number = 0;
  while (number < neighM->getNMaxi())
  {
    for (int isect = 0; isect < neighM->getNSect(); isect++)
    {
      if (_movingIsect[isect] >= _movingNsect[isect]) continue;
      _movingIsect[isect]++;
      number++;
      if (number >= neighM->getNMaxi()) break;
    }
  }

  /* Discard the data beyond the admissible rank per sector */

  for (int isect = 0; isect < neighM->getNSect(); isect++)
  {
    if (_movingIsect[isect] >= _movingNsect[isect]) continue;
    number = 0;
    for (int i = 0; i < nsel; i++)
    {
      int j = _movingInd[i];
      int jsect = ranks[j];
      if (jsect < 0) continue;
      if (isect != jsect) continue;
      number++;
      if (number > _movingIsect[isect]) ranks[j] = -1;
    }
  }
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
void NeighWork::_display(const VectorInt& ranks)
{
  String string;
  int ndim = _dbin->getNDim();
  int nech = _dbin->getSampleNumber();
  bool flag_ext = _dbin->getBlockExtensionNumber() > 0;

  /* Neighborhood data */

  mestitle(1, "Data selected in neighborhood");
  tab_prints(NULL, "Rank");
  tab_prints(NULL, "Sample");
  if (_dbin->hasCode()) tab_prints(NULL, "Code");
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
  if (_neighParam->getType() == ENeigh::MOVING)
    tab_prints(NULL, "Sector");
  message("\n");

  /* Loop on the sample points */

  int nsel = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    if (ranks[iech] < 0) continue;

    tab_printi(NULL, nsel + 1);
    tab_printi(NULL, iech + 1);
    if (_dbin->hasCode())
      tab_printi(NULL, static_cast<int>(_dbin->getCode(iech)));
    for (int idim = 0; idim < ndim; idim++)
      tab_printg(NULL, _dbin->getCoordinate(iech, idim));
    if (flag_ext)
    {
      for (int idim = 0; idim < ndim; idim++)
        tab_printg(NULL, _dbin->getBlockExtension(iech, idim));
    }
    if (_neighParam->getType() == ENeigh::MOVING)
      tab_printi(NULL, ranks[iech] + 1);
    message("\n");
    nsel++;
  }
  return;
}

void NeighWork::_checkUnchanged(const Db *dbout,
                                int iech_out,
                                const VectorInt &ranks)
{
  VectorInt rsorted = ranks;
  if (ranks.size() > 0)
  {
    std::sort(rsorted.begin(), rsorted.end());
  }

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

  _dbout = dbout;
  _iechOut = iech_out;
  _nbghMemo = rsorted;
}

void NeighWork::_clearMemory()
{
  _dbout = nullptr;
  _iechOut = -1;
  _nbghMemo = VectorInt();
}

/**
 * Checks if the current target matches the target previously treated
 * in the same procedure. If match is reached, then there is no need
 * to compute a new neighborhood: use the previous Vector of sample ranks.
 * Store the references of the new 'dbout' and 'iech_out' for next optimizations
 * @param dbout    Current 'Db' structure for output
 * @param iech_out Rank of the current target sample
 * @param ranks    Vector of selected samples
 * @param verbose  Verbose option
 * @return
 */
bool NeighWork::_isSameTarget(const Db* dbout,
                              int iech_out,
                              VectorInt& ranks,
                              bool verbose)
{
  // Check if the target remained unchanged
  bool flagSame = true;
  if (_dbout == nullptr || _iechOut < 0) flagSame = false;
  if (dbout != _dbout) flagSame = false;
  if (iech_out != _iechOut) flagSame = false;

  _resetFromMemory(flagSame, ranks, verbose);
  return flagSame;
}

bool NeighWork::_isSameTargetBench(const Db* dbout,
                                   int iech_out,
                                   VectorInt& ranks,
                                   bool verbose)
{
  // If no memorization is available, the match is false
  if (_dbout == nullptr || _iechOut < 0) return false;

  // Check if current target and previous target belong to the same bench

  bool flagSame = true;
  int ndim = dbout->getNDim();
  if (dbout->isGrid())
  {
    const DbGrid* dbgrid = dynamic_cast<const DbGrid*>(dbout);
    int nval = 1;
    for (int idim = 0; idim < ndim - 1; idim++)
      nval *= dbgrid->getNX(idim);
    if ((iech_out / nval) != (_iechOut / nval)) flagSame = false;
  }
  else
  {
    if (dbout->getCoordinate(iech_out, ndim - 1) != _dbout->getCoordinate(
        _iechOut, ndim - 1)) flagSame = false;
  }

  _resetFromMemory(flagSame, ranks, verbose);
  return flagSame;
}

bool NeighWork::_isSameTargetUnique(const Db* /*dbout*/,
                                    int /*iech_out*/,
                                    VectorInt& ranks,
                                    bool verbose)
{
  // If no memorization is available, the match is false
  if (_dbout == nullptr || _iechOut < 0) return false;
  bool flagSame = true;
  _resetFromMemory(flagSame, ranks, verbose);
  return flagSame;
}

void NeighWork::_resetFromMemory(bool flagSame, VectorInt& ranks, bool verbose)
{
  if (flagSame)
  {
    // If target is unchanged, upload the previously stored vector of Sample ranks

    ranks = _nbghMemo;
    _flagIsUnchanged = true;
    if (verbose)
      message(">>> Search is bypassed as already calculated\n");
  }
  else
  {
    if (verbose)
      message(">>> Search must probably be performed\n");
  }
}

/**
 * Update the set of selected samples in case of colocated option
 * This is done only if:
 * - the colocation option is ON (vector of colocated variable is defined)
 * - at least one of the colocated variables at the target is valid
 * - the target does not coincide with a sample already selected
 * If the colocation option is validated, an additional member is added to 'ranks':
 * it value is conventionally set to -1.
 * @param rankColCok Vector of Colocated Variables
 * @param ranks      Vector of samples already selected
 * @param iech_out   Rank of the targt site (in dbout)
 * @return
 */
void NeighWork::_updateColCok(const VectorInt &rankColCok,
                              VectorInt &ranks,
                              int iech_out)
{
  if (rankColCok.empty()) return;
  int nvarin = (int) rankColCok.size();

  /* Do not add the target if no variable is defined */
  bool found = false;
  for (int ivar = 0; ivar < nvarin && !found; ivar++)
  {
    int jvar = rankColCok[ivar];
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

/****************************************************************************/
/*!
 **  Returns the main Neighborhood Parameters for a given target as a vector:
 ** \li                    0 : Number of active samples
 ** \li                    1 : Minimum distance
 ** \li                    2 : Maximum distance
 ** \li                    3 : Number of non-empty sectors
 ** \li                    4 : Number of consecutive empty sectors
 **
 ** \param[in]  dbout         output Db structure
 ** \param[in]  iech_out      Valid Rank of the sample in the output Db
 ** \param[in]  rankColCok    Vector of ColCok information (optional)
 **
 *****************************************************************************/
VectorDouble NeighWork::summary(Db *dbout,
                                int iech_out,
                                const VectorInt& rankColCok)
{
  VectorDouble tab(5,0.);

  const NeighMoving* neighM = dynamic_cast<const NeighMoving*>(_neighParam);
  if (neighM == nullptr) return tab;

  /* Number of selected samples */

  VectorInt nbgh_ranks = select(dbout, iech_out, rankColCok);
  int nsel = (int) nbgh_ranks.size();
  tab[0] = (double) nsel;

  /* Maximum distance */

  double dmax = TEST;
  for (int iech = 0; iech < nsel; iech++)
  {
    double dist = _movingDst[iech];
    if (FFFF(dmax) || dist > dmax) dmax = dist;
  }
  tab[1] = dmax;

  /* Minimum distance */

  double dmin = TEST;
  for (int iech = 0; iech < nsel; iech++)
  {
    double dist = _movingDst[iech];
    if (FFFF(dmin) || dist < dmin) dmin = dist;
  }
  tab[2] = dmin;

  /* Number of sectors containing neighborhood information */

  int number = 0;
  for (int isect = 0; isect < neighM->getNSect(); isect++)
  {
    if (_movingNsect[isect] > 0) number++;
  }
  tab[3] = (double) number;

  /* Number of consecutive empty sectors */

  number = 0;
  int n_empty = 0;
  for (int isect = 0; isect < neighM->getNSect(); isect++)
  {
    if (_movingNsect[isect] > 0)
      n_empty = 0;
    else
    {
      n_empty++;
      if (n_empty > number) number = n_empty;
    }
  }
  if (neighM->getNSect() > 0)
  {
    if(_movingNsect[0] > 0)
      n_empty = 0;
    else
    {
      n_empty++;
      if (n_empty > number) number = n_empty;
    }
  }
  tab[4] = (double) number;

  return tab;
}
