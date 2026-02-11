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

#include "Anamorphosis/AAnam.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Model/Model.hpp"
#include "Polygon/Polygons.hpp"
#include "Stats/Classical.hpp"
#include "Variogram/VCloud.hpp"
#include "Variogram/Vario.hpp"
#include "Variogram/VarioParam.hpp"

namespace gstlrn
{
static Polygons* POLYGON = nullptr;
static VectorDouble IDS;

VCloud::VCloud(DbGrid* dbcloud, const VarioParam* varioparam)
  : AVario()
  , _dbcloud(dbcloud)
  , _varioparam(varioparam)
  , _IPTR(-1)
{
}

VCloud::VCloud(const VCloud& r)
  : AVario(r)
  , _dbcloud(r._dbcloud)
  , _varioparam(r._varioparam)
  , _IPTR(r._IPTR)
{
}

VCloud& VCloud::operator=(const VCloud& r)
{
  if (this != &r)
  {
    AVario::operator=(r);
    _dbcloud    = r._dbcloud;
    _varioparam = r._varioparam;
    _IPTR       = r._IPTR;
  }
  return *this;
}

VCloud::~VCloud()
{
}

double VCloud::_getIVAR(const Db* db, Id iech, Id ivar) const
{
  return db->getZVariable(iech, ivar);
}

/****************************************************************************/
/*!
 **  Internal function for setting a VCloud value
 **
 ** \param[in]  iech1       Rank of the first sample
 ** \param[in]  iech2       Rank of the second sample
 ** \param[in]  nvar        Number of variables
 ** \param[in]  idir        Rank of the direction
 ** \param[in]  ilag        Rank of the variogram lag
 ** \param[in]  ivar        Index of the first variable
 ** \param[in]  jvar        Index of the second variable
 ** \param[in]  orient      Orientation
 ** \param[in]  ww          Weight
 ** \param[in]  w1          Weight of the first sample
 ** \param[in]  w2          Weight of the second sample
 ** \param[in]  z1          Value of first variable
 ** \param[in]  z2          Value of second variable
 ** \param[in]  dist        Distance
 ** \param[in]  value       Variogram value
 **
 *****************************************************************************/
void VCloud::_setAVarioResult(Id iech1,
                              Id iech2,
                              Id nvar,
                              Id idir,
                              Id ilag,
                              Id ivar,
                              Id jvar,
                              Id orient,
                              double ww,
                              double w1,
                              double w2,
                              double z1,
                              double z2,
                              double dist,
                              double value)
{
  DECLARE_UNUSED(nvar);
  DECLARE_UNUSED(idir);
  DECLARE_UNUSED(ilag);
  DECLARE_UNUSED(ivar);
  DECLARE_UNUSED(jvar);
  DECLARE_UNUSED(orient);
  DECLARE_UNUSED(ww);
  DECLARE_UNUSED(w1);
  DECLARE_UNUSED(w2);
  DECLARE_UNUSED(z1);
  DECLARE_UNUSED(z2);

  Id igrid = _update_discretization_grid(dist, value);
  if (igrid < 0) return;

  if (POLYGON == nullptr)
  {
    // Store in the output grid
    _dbcloud->updArray(igrid, _IPTR, EOperator::ADD, 1.);
  }
  else
  {
    VectorInt indg(2);
    VectorDouble coor(2);
    _dbcloud->rankToIndice(igrid, indg);
    _dbcloud->indicesToCoordinateInPlace(indg, coor);
    if (!POLYGON->inside(coor, false)) return;
    {
      IDS[iech1] += 1.;
      IDS[iech2] += 1.;
    }
  }

  _storage(iech1, iech2, dist, value);
}

/****************************************************************************/
/*!
 **  Evaluate the experimental variogram cloud on irregular data
 **  This method creates one variable per direction of 'dirparam'
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db descriptor
 ** \param[in]  calculType  Calculation type
 ** \param[in]  flag_ergodic Ergodic flag
 ** \param[in]  namconv      Naming convention
 **
 *****************************************************************************/
Id VCloud::compute(Db* db,
                   const ECalcVario& calculType,
                   bool flag_ergodic,
                   const NamingConvention& namconv)
{
  if (db == nullptr) return (1);

  /* Preliminary checks */

  if (db->getNDim() != _varioparam->getNDim())
  {
    messerr("Inconsistent parameters:");
    messerr("Data Base: NDIM=%d", db->getNDim());
    messerr("Variogram: NDIM=%d", _varioparam->getNDim());
    return (1);
  }
  if (!db->isNVarComparedTo(1)) return 1;
  if (_dbcloud->getNDim() != 2)
  {
    messerr("The output Db for storing the variogram cloud must be 2-D");
    return (1);
  }

  /* Allocate new variables */

  setCalcul(calculType);
  setErgodic(flag_ergodic);
  if (getFlagAsym())
  {
    messerr("VCloud calculation only available for symmetrical variograms");
    return 1;
  }
  Id ndir = _varioparam->getNDir();
  Id iptr = _dbcloud->addColumnsByConstant(ndir, 0.);
  if (iptr < 0) return (1);

  /* Loop on the directions to evaluate */

  for (Id idir = 0; idir < ndir; idir++)
  {
    _IPTR = iptr + idir;
    _variogram_cloud(db, idir);
    _final_discretization_grid();
  }

  // Naming of the newly created variables

  namconv.setNamesAndLocators(db, VectorString(), ELoc::Z, -1, _dbcloud, iptr,
                              String(), ndir, false);

  return (0);
}

/****************************************************************************/
/*!
 **  Replace zero values by TEST values
 **
 *****************************************************************************/
void VCloud::_final_discretization_grid()
{
  for (Id iech = 0, nech = _dbcloud->getNSample(); iech < nech; iech++)
  {
    double value = _dbcloud->getArray(iech, _IPTR);
    if (value != 0.) continue;
    _dbcloud->setArray(iech, _IPTR, TEST);
  }
}

/****************************************************************************/
/*!
 **  Evaluate the variogram cloud
 **
 ** \param[in]  db      Db descriptor
 ** \param[in]  idir    Rank of the Direction
 **
 *****************************************************************************/
void VCloud::_variogram_cloud(Db* db, Id idir)
{
  double dist;
  SpaceTarget T1(_varioparam->getSpace());
  SpaceTarget T2(_varioparam->getSpace());

  // Creating a local Vario structure (to constitute the BiTargetCheck list
  Vario* vario = Vario::create(*_varioparam);
  vario->setDb(db);
  if (vario->prepare()) return;

  // Local variables to speed up calculations
  bool hasSel = db->hasLocVariable(ELoc::SEL);
  Id nech     = db->getNSample();
  Id nvar     = db->getNLoc(ELoc::Z);

  /* Loop on the first point */

  for (Id iech = 0; iech < nech - 1; iech++)
  {
    if (hasSel && !db->isActive(iech)) continue;
    db->getSampleAsSTInPlace(iech, T1);

    Id ideb = (_varioparam->isDateUsed(db)) ? 0 : iech + 1;
    for (Id jech = ideb; jech < nech; jech++)
    {
      if (hasSel && !db->isActive(jech)) continue;
      db->getSampleAsSTInPlace(jech, T2);

      // Reject the point as soon as one BiTargetChecker is not correct
      if (!vario->keepPair(idir, T1, T2, &dist)) continue;

      (this->*_evaluate)(db, nvar, iech, jech, idir, 0, dist, false);
    }
  }

  delete vario;
}

/****************************************************************************/
/*!
 **  Add one pick to the discretization grid
 **
 ** \return  Index of the grid cell (or -1)
 **
 ** \param[in]  x     Coordinate along the first axis
 ** \param[in]  y     Coordinate along the first axis
 **
 *****************************************************************************/
Id VCloud::_update_discretization_grid(double x, double y)
{
  Id ndim = _dbcloud->getNDim();
  VectorInt indg(ndim, 0);

  Id ix = static_cast<Id>(floor((x - _dbcloud->getX0(0)) / _dbcloud->getDX(0) + 0.5));
  Id iy = static_cast<Id>(floor((y - _dbcloud->getX0(1)) / _dbcloud->getDX(1) + 0.5));
  if (ix < 0 || ix >= _dbcloud->getNX(0)) return (-1);
  if (iy < 0 || iy >= _dbcloud->getNX(1)) return (-1);
  indg[0] = ix;
  indg[1] = iy;
  return _dbcloud->indiceToRank(indg);
}

DbGrid* vcloudGrid(const Db* db, double lagmax, double varmax, Id lagnb, Id varnb)
{
  if (FFFF(lagmax)) lagmax = db->getExtensionDiagonal();
  if (FFFF(varmax)) varmax = 3. * db->getVariance(db->getNameByLocator(ELoc::Z));

  // Create a grid as a support for the variogram cloud calculations

  VectorInt nx(2);
  nx[0] = lagnb;
  nx[1] = varnb;
  VectorDouble dx(2);
  dx[0] = lagmax / static_cast<double>(lagnb);
  dx[1] = varmax / static_cast<double>(varnb);
  VectorDouble x0(2);
  x0[0]          = 0.;
  x0[1]          = 0.;
  DbGrid* dbgrid = DbGrid::create(nx, dx, x0);
  return dbgrid;
}

/****************************************************************************/
/*!
 **  Evaluate the experimental variogram cloud
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db descriptor
 ** \param[in]  varioparam   VarioParam structure
 ** \param[in]  calculType   Calculation type
 ** \param[in]  flag_ergodic Ergodic flag
 ** \param[in]  lagmax       Maximum distance
 ** \param[in]  varmax       Maximum Variance value (see remarks)
 ** \param[in]  lagnb        Number of discretization steps along distance axis
 ** \param[in]  varnb        Number of discretization steps along variance axis
 ** \param[in]  namconv      Naming convention
 **
 ** \remarks If 'lagmax' is not provided, it is set to the diagonal of the
 ** area covered by the active samples within 'db'.
 ** \remarks If 'varmax' is not defined, it is set to 3 times the experimental
 ** variance of the first variable (Z_locator)
 **
 *****************************************************************************/
DbGrid* db_vcloud(Db* db,
                  const VarioParam* varioparam,
                  const ECalcVario& calculType,
                  bool flag_ergodic,
                  double lagmax,
                  double varmax,
                  Id lagnb,
                  Id varnb,
                  const NamingConvention& namconv)
{
  // Create a grid as a support for the variogram cloud calculations

  DbGrid* dbgrid = vcloudGrid(db, lagmax, varmax, lagnb, varnb);

  // Initialize the variogram cloud structure

  VCloud vcloud(dbgrid, varioparam);

  // Calling the variogram cloud calculation function

  if (vcloud.compute(db, calculType, flag_ergodic, namconv))
  {
    delete dbgrid;
    dbgrid = nullptr;
  }
  return dbgrid;
}

/****************************************************************************/
/*!
 **  Check the samples which are involved in the pairs which are located
 **  within the polygon
 **
 ** \param[in]  db      Db descriptor
 ** \param[in]  polygon Polygons structure
 ** \param[in]  idir    Rank of the direction of itnerest
 **
 *****************************************************************************/
Id VCloud::selectFromPolygon(Db* db, Polygons* polygon, Id idir)
{
  POLYGON = polygon;
  Id nech = db->getNSample();
  IDS.resize(nech, 0.);

  _variogram_cloud(db, idir);

  /* Printout the scores: they are ranked by decreasing number */

  mestitle(0, "Samples in variogram cloud (by decreasing order of occurence)");
  VectorInt rank = VH::sequence(nech);
  ut_sort_double(0, nech, rank.data(), IDS.data());

  for (Id iech = 0; iech < nech; iech++)
  {
    Id jech = nech - iech - 1;
    if (IDS[jech] <= 0.) break;
    message("Sample #%3d: %d occurence(s)\n", rank[jech] + 1, static_cast<Id>(IDS[jech]));
  }

  POLYGON = nullptr;
  IDS.clear();
  return 0;
}

String VCloud::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  // Print the calculation type

  sstr << _elemString(strfmt) << std::endl;

  // TODO: to be completed

  return sstr.str();
}

DbGrid* vcloudCalculate(Db* db,
                        const ECalcVario& calculType,
                        bool flag_ergodic,
                        Id nlag,
                        double dlag,
                        Id ndir,
                        const VectorDouble& angles,
                        double toldis,
                        double tolang,
                        Id lagnb,
                        Id varnb)
{
  VarioParam* varioparam = nullptr;
  if (ndir > 1)
    varioparam = VarioParam::createMultiple(ndir, nlag, dlag, toldis);
  else if (!angles.empty())
    varioparam = VarioParam::createSeveral2D(angles, nlag, dlag, toldis, tolang);
  else
    varioparam = VarioParam::createOmniDirection(nlag, dlag, toldis);

  auto* dbgrid = db_vcloud(db, varioparam, calculType, flag_ergodic, TEST, TEST, lagnb, varnb);

  return dbgrid;
}

/**
 * @brief Dump the storage of variogram pair calculations
 *
 * @param mode See remarks
 * @param distmax  Maximum distance for pairs to be considered (mode=1 or 2)
 * @param varmin   Minimum variogram value for pairs to be considered (mode=2)
 * @param npairmax Maximum number of pairs to be displayed (mode=1)
 * @param ndatamax Maximum number of data to be displayed (mode=2)
 *
 * @remarks The argument 'mode' can be:
 * - 0 : Full dump of all the pairs
 * - 1 : Printout of the pairs (up to 'nparimax') with largest variogram values for a distance lower then 'distmax'
 * - 2 : Printout of the sample mostly involved in large varioagrma value (> 'varmin') for a distance lower then 'distmax'
 */
void VCloud::dumpStorage(Id mode, double distmax, double varmin, Id npairmax, Id ndatamax)
{
  if (!_isStorage()) return;

  mestitle(0, "Storage of variogram pair calculations");
  Id npairs = static_cast<Id>(_getStorageSize());
  message("Number of pairs = %d\n", npairs);

  if (mode == 0)
  {
    // Provide the full dump of all the pairs
    for (Id ipair = 0; ipair < npairs; ++ipair)
    {
      const auto& record = _getStorage(ipair);
      message("Pair (%3d, %3d) - Dist: %lf - Value: %lf\n",
              static_cast<Id>(record[0]), static_cast<Id>(record[1]), record[2], record[3]);
    }
  }
  else if (mode == 1)
  {
    // Provide only the summary of the pairs

    if (!FFFF(distmax))
      message("- Pairs are considered up to a Maximum Distance = %lf\n", distmax);
    if (!IFFFF(npairmax))
      message("- Only the %d largest variogram values are displayed\n", npairmax);

    VectorDouble values;
    VectorInt ipairs;
    for (Id ipair = 0; ipair < npairs; ++ipair)
    {
      double distloc = _getStorage(ipair)[2];
      double valloc  = _getStorage(ipair)[3];
      if (!FFFF(distmax) && distloc > distmax) continue;
      values.push_back(valloc);
      ipairs.push_back(ipair);
    }
    VectorInt order = VH::orderRanks(values, false);

    // Dump the pairs in descending order of the value
    Id nstat = order.size();
    if (!IFFFF(npairmax)) nstat = MIN(npairmax, nstat);
    for (Id irank = 0; irank < nstat; ++irank)
    {
      Id jrank           = order[irank];
      const auto& record = _getStorage(ipairs[jrank]);
      message("Pair (%3d, %3d) - Dist: %lf - Value: %lf\n",
              static_cast<Id>(record[0]), static_cast<Id>(record[1]), record[2], record[3]);
    }
  }
  else if (mode == 2)
  {
    // Provide only the summary of the pairs

    if (!FFFF(distmax))
      message("- Pairs are considered up to a Maximum Distance = %lf\n", distmax);
    if (varmin > 0.)
      message("- Pairs are considered with a Minimum Variogram Value = %lf\n", varmin);

    // Find the largest sample number;
    Id nbdata = -1;
    for (Id ipair = 0; ipair < npairs; ++ipair)
    {
      double iech1 = static_cast<Id>(_getStorage(ipair)[0]);
      if (iech1 > nbdata) nbdata = iech1;
      double iech2 = static_cast<Id>(_getStorage(ipair)[1]);
      if (iech2 > nbdata) nbdata = iech2;
    }
    VectorInt count(nbdata + 1, 0.);

    for (Id ipair = 0; ipair < npairs; ++ipair)
    {
      double distloc = _getStorage(ipair)[2];
      double valloc  = _getStorage(ipair)[3];
      if (!FFFF(distmax) && distloc > distmax) continue;
      if (varmin > 0. && valloc < varmin) continue;
      double iech1 = static_cast<Id>(_getStorage(ipair)[0]);
      count[iech1]++;
      double iech2 = static_cast<Id>(_getStorage(ipair)[1]);
      count[iech2]++;
    }
    VectorInt order = VH::orderRanks(count, false);

    Id nstat = count.size();
    if (!IFFFF(ndatamax)) nstat = MIN(ndatamax, nstat);
    for (Id irank = 0; irank < nstat; ++irank)
    {
      Id jrank = order[irank];
      message("Sample %3d - Number of calls: %d\n", jrank, count[jrank]);
    }
  }
}

} // namespace gstlrn