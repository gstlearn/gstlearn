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
#include "Neigh/NeighWork.hpp"
#include "Neigh/Neigh.hpp"
#include "Db/Db.hpp"
#include "Basic/Vector.hpp"
#include "Basic/AException.hpp"
#include "geoslib_f.h"
#include "geoslib_old_f.h"

#include <math.h>

NeighWork::NeighWork(const Neigh* neigh,
                     const Db* dbin,
                     bool flag_var_nocheck,
                     bool flag_simu)
    : _neigh(),
      _dbin(),
      _flagInitialized(false),
      _nbghInd(),
      _nbghIsect(),
      _nbghNsect(),
      _nbghX1(),
      _nbghX2(),
      _nbghDst(),
      _flagVariableNoCheck(false),
      _flagSimu(false)
{
  initialize(neigh, dbin, flag_var_nocheck, flag_simu);
}

NeighWork::NeighWork(const NeighWork& r)
    : _neigh(r._neigh),
      _dbin(r._dbin),
      _flagInitialized(r._flagInitialized),
      _nbghInd(r._nbghInd),
      _nbghIsect(r._nbghIsect),
      _nbghNsect(r._nbghNsect),
      _nbghX1(r._nbghX1),
      _nbghX2(r._nbghX2),
      _nbghDst(r._nbghDst),
      _flagVariableNoCheck(r._flagVariableNoCheck),
      _flagSimu(r._flagSimu)
{
}

NeighWork& NeighWork::operator=(const NeighWork& r)
{
  if (this != &r)
  {
    _neigh = r._neigh;
    _dbin = r._dbin;
    _flagInitialized = r._flagInitialized;
    _nbghInd = r._nbghInd;
    _nbghIsect = r._nbghIsect;
    _nbghNsect = r._nbghNsect;
    _nbghX1 = r._nbghX1;
    _nbghX2 = r._nbghX2;
    _nbghDst = r._nbghDst;
    _flagVariableNoCheck = r._flagVariableNoCheck;
    _flagSimu = r._flagSimu;
   }
  return *this;
}

NeighWork::~NeighWork()
{
}

/****************************************************************************/
/*!
 **  Initialilze the neighborhood search
 **
 ** \param[in]  neigh         Description of the Neigh parameters
 ** \param[in]  dbin          input Db structure
 ** \param[in]  flag_simu     1 if used for Simulation
 ** \param[in]  flag_var_nocheck 1 if no check on variable definition
 **
 ** \remarks When the Neighborhood is performed in the case of Simulations
 ** \remarks checking for all variables being undefined is performed
 ** \remarks on ELoc::SIMU rather than on ELoc::Z
 **
 *****************************************************************************/
void NeighWork::initialize(const Neigh* neigh, const Db* dbin,
                           bool flag_var_nocheck, bool flag_simu)
{
  if (neigh == nullptr || dbin == nullptr) return;
  _neigh = neigh;
  _dbin = dbin;
  _flagVariableNoCheck = flag_var_nocheck;
  _flagSimu = flag_simu;

  int nech  = _dbin->getSampleNumber();
  int ndim  = _dbin->getNDim();
  int nsect = _neigh->getNSect();

  _nbghInd = VectorInt(nech);
  _nbghDst = VectorDouble(nech);
  _nbghIsect = VectorInt(nsect);
  _nbghNsect = VectorInt(nsect);
  _nbghX1 = VectorDouble(ndim);
  _nbghX2 = VectorDouble(ndim);

  _flagInitialized = true;
}

void NeighWork::clear()
{
  /* Initialization */

  if (! isInitialized()) return;

  // Unset the pointers

  _neigh = nullptr;
  _dbin  = nullptr;

  /* Core deallocation */

  _nbghInd.clear();
  _nbghDst.clear();
  _nbghIsect.clear();
  _nbghNsect.clear();
  _nbghX1.clear();
  _nbghX2.clear();

  _flagInitialized = false;
}

/****************************************************************************/
/*!
 **  Select the neighborhood
 **
 ** \return  Vector of sample ranks in neighborhood (empty when error)
 **
 ** \param[in]  dbout         output Db structure
 ** \param[in]  iech_out      rank of the sample in the output Db
 **
 *****************************************************************************/
VectorInt NeighWork::select(Db *dbout,
                            int iech_out)
{
  VectorInt ranks;
  if (! isInitialized()) return ranks;
  int nech = _dbin->getSampleNumber();

  ranks.resize(nech, -1);

  /* Select the active data points */

  switch (_neigh->getType().toEnum())
  {
    case ENeigh::E_IMAGE:
    case ENeigh::E_UNIQUE:
      _unique(dbout, iech_out, ranks);
      break;

    case ENeigh::E_BENCH:
      _bench(dbout, iech_out, ranks);
      break;

    case ENeigh::E_MOVING:
      if (_moving(dbout, iech_out, ranks)) return VectorInt();
      break;
  }

  /* Print the Neighborhood search result */

  if (debug_query("nbgh")) _display(ranks);

  /* Compress the vector of returned sample ranks */

  int necr = 0;
  for (int iech = 0; iech < nech; iech++)
    if (ranks[iech] >= 0) ranks[necr++] = iech;
  ranks.resize(necr);

  return ranks;
}

/****************************************************************************/
/*!
 **  Select the unique neighborhood (or Image Neighborhood)
 **  Pay attention to the cross-validation flag
 **
 ** \param[in]  dbout     output Db structure
 ** \param[in]  iech_out  rank of the output sample
 **
 ** \param[out]  ranks    Neighborhood selection array
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

    if (_neigh->getFlagXvalid() != 0)
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
 ** \param[out]  ranks     Neighborhood selection array
 **
 *****************************************************************************/
void NeighWork::_bench(Db *dbout, int iech_out, VectorInt& ranks)
{
  int nech = _dbin->getSampleNumber();
  int idim_bench = _dbin->getNDim() - 1;
  double z0 = dbout->getCoordinate(iech_out, idim_bench);

  /* Loop on samples */

  for (int iech = 0; iech < nech; iech++)
  {
    /* Discard the masked input sample */

    if (! _dbin->isActive(iech)) continue;

    /* Discard samples where all variables are undefined */

    if (_discardUndefined(iech)) continue;

    /* Discard the target sample for the cross-validation option */

    if (_neigh->getFlagXvalid() != 0)
    {
      if (_xvalid(dbout, iech, iech_out)) continue;
    }

    /* Discard sample located outside the bench */

    if (ABS(_dbin->getCoordinate(iech,idim_bench) - z0) <= _neigh->getWidth())
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
 ** \param[out]  ranks    Array of active data point ranks
 **
 *****************************************************************************/
int NeighWork::_moving(Db *dbout, int iech_out, VectorInt& ranks, double eps)
{
  int nech = _dbin->getSampleNumber();
  int isect = 0;
  if (nech < _neigh->getNMini()) return 1;

  /* Loop on the data points */

  double distmax = 0.;
  int nsel = 0;
  for (int iech = nsel = 0; iech < nech; iech++)
  {

    /* Discard the masked input sample */

    if (! _dbin->isActive(iech)) continue;

    /* Discard samples where all variables are undefined */

    if (_discardUndefined(iech)) continue;

    /* Discard the target sample for the cross-validation option */

    if (_neigh->getFlagXvalid() != 0)
    {
      if (_xvalid(dbout, iech, iech_out)) continue;
    }

    /* Calculate the distance between data and target */

    double dist = _movingDist(dbout, iech, iech_out);
    if (! FFFF(_neigh->getRadius()) && dist > _neigh->getRadius()) continue;
    if (dist > distmax) distmax = dist;

    /* Calculate the angular sector to which the sample belongs */

    if (_neigh->getFlagSector())
      isect = _movingSectorDefine(_nbghX1[0], _nbghX1[1]);

    /* The sample may be selected */

    _nbghInd[nsel] = iech;
    _nbghDst[nsel] = dist;
    ranks[iech] = isect;
    nsel++;
  }
  if (nsel < _neigh->getNMini()) return 1;

  /* Slightly modify the distances in order to ensure the sorting results */
  /* In the case of equal distances                                       */

  for (int isel = 0; isel < nsel; isel++)
    _nbghDst[isel] += distmax * isel * eps;

  /* Sort the selected samples according to the distance */

  ut_sort_double(0, nsel, _nbghInd.data(), _nbghDst.data());

  /* For each angular sector, select the first sample up to the maximum */

  if (_neigh->getFlagSector() && _neigh->getNSMax() > 0)
  {
    _movingSectorNsmax(nsel, ranks);
    if (nsel < _neigh->getNMini()) return 1;
  }

  /* Select the first data samples */

  _movingSelect(nsel, ranks);

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
  int nvar = _dbin->getVariableNumber();

  if (_flagVariableNoCheck)
  {
    if (_dbin->getVariableNumber() <= 0)
      return 0;
    else
    {
      for (int ivar = 0; ivar < nvar; ivar++)
        if (! FFFF(_dbin->getVariable(iech, ivar))) return 0;
      return 1;
    }
  }

  if (! _flagSimu)
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      if (!FFFF(_dbin->getVariable(iech, ivar))) return 0;
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
  if (_neigh->getFlagXvalid() == 0)
    return 0;
  else if (_neigh->getFlagXvalid() > 0)
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
  int ndim = _dbin->getNDim();

  /* Calculate the distance to the target */

  for (int idim = 0; idim < ndim; idim++)
    _nbghX1[idim] = dbout->getCoordinate(iech_out, idim)
        - _dbin->getCoordinate(iech_in, idim);

  /* Anisotropic neighborhood */

  if (_neigh->getFlagAniso())
  {

    /* Rotated anisotropy ellipsoid */

    if (_neigh->getFlagRotation())
    {
      matrix_product(1, ndim, ndim, _nbghX1.data(),
                     _neigh->getAnisoRotMats().data(), _nbghX2.data());
      _nbghX1 = _nbghX2;
    }
    for (int idim = 0; idim < ndim; idim++)
      _nbghX1[idim] /= _neigh->getAnisoCoeff(idim);
  }

  /* Calculate the distance */

  double dist;
  matrix_product(1, ndim, 1, _nbghX1.data(), _nbghX1.data(), &dist);
  dist = sqrt(dist);
  return dist;
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
  double angle;

  int isect = 0;
  if (_neigh->getNSect() > 1)
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
    isect = (int) (_neigh->getNSect() * angle / (2. * GV_PI));
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
  int n_end = 0;
  for (int isect = 0; isect < _neigh->getNSect(); isect++)
  {
    int n_ang = 0;
    for (int i = 0; i < nsel; i++)
    {
      int j = _nbghInd[i];
      if (ranks[j] != isect) continue;
      if (n_ang < _neigh->getNSMax())
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
  int number;

  if (_neigh->getNMaxi() <= 0) return;
  for (int isect = 0; isect < _neigh->getNSect(); isect++)
    _nbghNsect[isect] = _nbghIsect[isect] = 0;

  /* Count the number of samples per sector */

  number = 0;
  for (int i = 0; i < nsel; i++)
  {
    int j = _nbghInd[i];
    int isect = ranks[j];
    if (isect < 0) continue;
    _nbghNsect[isect]++;
    number++;
  }
  if (number < _neigh->getNMaxi()) return;

  /* Find the rank of the admissible data per sector */

  number = 0;
  while (number < _neigh->getNMaxi())
  {
    for (int isect = 0; isect < _neigh->getNSect(); isect++)
    {
      if (_nbghIsect[isect] >= _nbghNsect[isect]) continue;
      _nbghIsect[isect]++;
      number++;
      if (number >= _neigh->getNMaxi()) break;
    }
  }

  /* Discard the data beyond the admissible rank per sector */

  for (int isect = 0; isect < _neigh->getNSect(); isect++)
  {
    if (_nbghIsect[isect] >= _nbghNsect[isect]) continue;
    number = 0;
    for (int i = 0; i < nsel; i++)
    {
      int j = _nbghInd[i];
      int jsect = ranks[j];
      if (jsect < 0) continue;
      if (isect != jsect) continue;
      number++;
      if (number > _nbghIsect[isect]) ranks[j] = -1;
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
  bool flag_ext = is_flag_data_disc_defined();

  /* Neighborhood data */

  mestitle(1, "Data selected in neighborhood");
  tab_prints(NULL, 1, EJustify::RIGHT, "Rank");
  tab_prints(NULL, 1, EJustify::RIGHT, "Sample");
  if (_dbin->hasCode()) tab_prints(NULL, 1, EJustify::RIGHT, "Code");
  for (int idim = 0; idim < ndim; idim++)
  {
    string = getLocatorName(ELoc::X, idim);
    tab_prints(NULL, 1, EJustify::RIGHT, string.c_str());
  }
  if (flag_ext)
  {
    for (int idim = 0; idim < ndim; idim++)
    {
      string = getLocatorName(ELoc::BLEX, idim);
      tab_prints(NULL, 1, EJustify::RIGHT, string.c_str());
    }
  }
  if (_neigh->getType() == ENeigh::MOVING)
    tab_prints(NULL, 1, EJustify::RIGHT, "Sector");
  message("\n");

  /* Loop on the sample points */

  int nsel = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    if (ranks[iech] < 0) continue;

    tab_printi(NULL, 1, EJustify::RIGHT, nsel + 1);
    tab_printi(NULL, 1, EJustify::RIGHT, iech + 1);
    if (_dbin->hasCode())
      tab_printi(NULL, 1, EJustify::RIGHT,
                 static_cast<int>(_dbin->getCode(iech)));
    for (int idim = 0; idim < ndim; idim++)
      tab_printg(NULL, 1, EJustify::RIGHT, _dbin->getCoordinate(iech, idim));
    if (flag_ext)
    {
      for (int idim = 0; idim < ndim; idim++)
        tab_printg(NULL, 1, EJustify::RIGHT,
                   _dbin->getBlockExtension(iech, idim));
    }
    if (_neigh->getType() == ENeigh::MOVING)
      tab_printi(NULL, 1, EJustify::RIGHT, ranks[iech] + 1);
    message("\n");
    nsel++;
  }
  return;
}
