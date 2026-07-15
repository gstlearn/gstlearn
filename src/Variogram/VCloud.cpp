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
  VCloud::VCloud(DbGrid* dbcloud, const VarioParam* varioparam)
    : AVario()
    , _dbcloud(dbcloud)
    , _varioparam(varioparam)
    , _IPTR(-1)
    , _samples()
    , _flagSelection(false)
    , _IDS()
    , _polygon(nullptr)
    , _distmax(TEST)
    , _varmin(TEST)
    , _countmax(ITEST)
  {
  }

  VCloud::VCloud(const VCloud& r)
    : AVario(r)
    , _dbcloud(r._dbcloud)
    , _varioparam(r._varioparam)
    , _IPTR(r._IPTR)
    , _samples(r._samples)
    , _flagSelection(r._flagSelection)
    , _IDS(r._IDS)
    , _polygon(r._polygon)
    , _distmax(r._distmax)
    , _varmin(r._varmin)
    , _countmax(r._countmax)
  {
  }

  VCloud& VCloud::operator=(const VCloud& r)
  {
    if (this != &r)
    {
      AVario::operator=(r);
      _dbcloud = r._dbcloud;
      _varioparam = r._varioparam;
      _IPTR = r._IPTR;
      _samples = r._samples;
      _flagSelection = r._flagSelection;
      _IDS = r._IDS;
      _polygon = r._polygon;
      _distmax = r._distmax;
      _varmin = r._varmin;
      _countmax = r._countmax;
    }
    return *this;
  }

  VCloud::~VCloud() {}

  double VCloud::_getIVAR(const Db* db, Id iech, Id ivar) const
  {
    return db->getZVariable(iech, ivar);
  }

  void VCloud::setPolygon(const Polygons* polygon)
  {
    if (polygon == nullptr) return;
    _polygon = polygon;
    _flagSelection = true;
  }

  void VCloud::setIntervals(double distmax, double varmin, Id countmax)
  {
    if (isNA(distmax) && isNA(varmin)) return;
    _flagSelection = true;
    _distmax = distmax;
    _varmin = varmin;
    _countmax = countmax;
  }

  void VCloud::setSamples(const VectorInt& samples)
  {
    _samples = samples;
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
   **IDS
   *****************************************************************************/
  void VCloud::_setAVarioResult(
    Id iech1,
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

    Id igrid = _getDiscretizedCellRank(dist, value);
    if (igrid < 0) return;

    // Update the output grid
    _dbcloud->updArray(igrid, _IPTR, EOperator::ADD, 1.);

    if (!_flagSelection) return;

    if (_flagByPolygon())
    {
      VectorInt indg(2);
      VectorDouble coor(2);
      _dbcloud->rankToIndice(igrid, indg);
      _dbcloud->indicesToCoordinateInPlace(indg, coor);
      if (!_polygon->inside(coor, false)) return;
    }
    else if (_flagByIntervals())
    {
      if (!isNA(_distmax) && dist > _distmax) return;
      if (!isNA(_varmin) && value < _varmin) return;
    }

    _IDS[iech1] += 1;
    _IDS[iech2] += 1;

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
  Id VCloud::compute(
    Db* db,
    const ECalcVario& calculType,
    bool flag_ergodic,
    const NamingConvention& namconv)
  {
    if (db == nullptr) return (1);

    // Preliminary checks
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

    // Allocate new variables
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

    // Loop on the directions
    Id nech = db->getNSample();
    if (_flagSelection) _IDS.resize(nech, 0.);
    for (Id idir = 0; idir < ndir; idir++)
    {
      _IPTR = iptr + idir;

      if (_flagSelection) _IDS.fill(0.);

      if (_flagBySamples())
        _variogramCloudBySamples(db, idir);
      else
        _variogramCloud(db, idir);

      _finalDiscretizationOnGrid();

      if (_flagSelection)
      {

        // Printout the scores: they are ranked by decreasing number
        mestitle(
          0,
          "Samples for Direction %d (by decreasing order of "
          "occurence)",
          idir + 1);
        VectorInt ranks = VH::orderRanks(_IDS, false);
        for (Id iech = 0; iech < nech; iech++)
        {
          if (_countmax > 0 && iech >= _countmax) break;
          Id jech = ranks[iech];
          if (_IDS[jech] <= 0.) break;
          message(
            " Sample #%3d/%3d: %d occurence(s)\n", jech + 1, nech, _IDS[jech]);
        }
        message("\n");
      }
    }

    // Naming of the newly created variables
    auto names = db->getNamesByLocator(ELoc::Z);
    namconv.setOutput(names, 0, _dbcloud, iptr, String(), ndir, false);

    return (0);
  }

  /****************************************************************************/
  /*!
   **  Replace zero values by TEST values
   **
   *****************************************************************************/
  void VCloud::_finalDiscretizationOnGrid()
  {
    for (Id iech = 0, nech = _dbcloud->getNSample(); iech < nech; iech++)
    {
      if (_dbcloud->getArray(iech, _IPTR) != 0.) continue;
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
  void VCloud::_variogramCloud(Db* db, Id idir)
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
    Id nech = db->getNSample();
    Id nvar = db->getNLoc(ELoc::Z);

    // Loop on the first point
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

        // Add the contribution of the pair to the Variogram Cloud
        (this->*_evaluate)(db, nvar, iech, jech, idir, 0, dist, false);
      }
    }

    delete vario;
  }

  /****************************************************************************/
  /*!
   **  Evaluate the variogram cloud
   **
   ** \param[in]  db      Db descriptor
   ** \param[in]  idir    Rank of the Direction
   **
   *****************************************************************************/
  void VCloud::_variogramCloudBySamples(Db* db, Id idir)
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
    Id nech = db->getNSample();
    Id nvar = db->getNLoc(ELoc::Z);
    Id nsamples = static_cast<Id>(_samples.size());
    Id iech;

    // Loop on the first point
    for (Id iiech = 0; iiech < nsamples; iiech++)
    {
      iech = _samples[iiech];
      if (hasSel && !db->isActive(iech)) continue;
      db->getSampleAsSTInPlace(iech, T1);

      Id ideb = (_varioparam->isDateUsed(db)) ? 0 : iech + 1;
      for (Id jech = ideb; jech < nech; jech++)
      {
        if (hasSel && !db->isActive(jech)) continue;
        db->getSampleAsSTInPlace(jech, T2);

        // Reject the point as soon as one BiTargetChecker is not correct
        if (!vario->keepPair(idir, T1, T2, &dist)) continue;

        // Add the contribution of the pair to the Variogram Cloud
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
  Id VCloud::_getDiscretizedCellRank(double x, double y)
  {
    Id ndim = _dbcloud->getNDim();
    VectorInt indg(ndim, 0);

    Id ix = static_cast<Id>(
      floor((x - _dbcloud->getX0(0)) / _dbcloud->getDX(0) + 0.5));
    if (ix < 0 || ix >= _dbcloud->getNX(0)) return (-1);
    Id iy = static_cast<Id>(
      floor((y - _dbcloud->getX0(1)) / _dbcloud->getDX(1) + 0.5));
    if (iy < 0 || iy >= _dbcloud->getNX(1)) return (-1);
    indg[0] = ix;
    indg[1] = iy;
    return _dbcloud->indiceToRank(indg);
  }

  String VCloud::toString(const AStringFormat* strfmt) const
  {
    std::stringstream sstr;

    // Print the calculation type

    sstr << _elemString(strfmt) << std::endl;

    // TODO: to be completed

    return sstr.str();
  }

} // namespace gstlrn
