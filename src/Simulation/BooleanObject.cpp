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
#include "Simulation/BooleanObject.hpp"
#include "Basic/Law.hpp"
#include "Boolean/AShape.hpp"
#include "Boolean/ModelBoolean.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Simulation/SimuBooleanParam.hpp"

#include <cmath>

namespace gstlrn
{
  BooleanObject::BooleanObject(const AShape* ashape, Id ndim)
    : AStringable()
    , _mode(0)
    , _ndim(ndim)
    , _token(ashape)
    , _center({0., 0., 0.})
    , _extension(
        {
          0.,
          0.,
          0.,
        })
    , _orientation(0.)
    , _values({0., 0., 0.})
    , _box()
  {
  }

  BooleanObject::BooleanObject(const BooleanObject& r)
    : AStringable(r)
    , _mode(r._mode)
    , _ndim(r._ndim)
    , _token(r._token)
    , _center(r._center)
    , _extension(r._extension)
    , _orientation(r._orientation)
    , _values(r._values)
    , _box(r._box)
  {
  }

  BooleanObject& BooleanObject::operator=(const BooleanObject& r)
  {
    if (this != &r)
    {
      AStringable::operator=(r);
      _mode = r._mode;
      _ndim = r._ndim;
      _token = r._token;
      _center = r._center;
      _extension = r._extension;
      _orientation = r._orientation;
      _values = r._values;
      _box = r._box;
    }
    return *this;
  }

  BooleanObject::~BooleanObject() {}

  String BooleanObject::toString(const AStringFormat* /*strfmt*/) const
  {
    std::stringstream sstr;

    if (_mode == 0)
      sstr << "Tentative Object" << std::endl;
    else if (_mode == 1)
      sstr << "Primary Object" << std::endl;
    else if (_mode == 2)
      sstr << "Secondary Object" << std::endl;
    else
      sstr << "Object with unknown mode = " << _mode << std::endl;
    sstr << "- Type        = " << _token->getType().getDescr() << std::endl;
    constvect center_ndim(_center.begin(), _center.begin() + _ndim);
    sstr << "- Center      = " << toStrVectorVec(String(), center_ndim);
    constvect extension_ndim(_extension.begin(), _extension.begin() + _ndim);
    sstr << "- Extension   = " << toStrVectorVec(String(), extension_ndim);
    sstr << "- Orientation = " << _orientation << std::endl;

    return sstr.str();
  }

  void BooleanObject::setCenter(const VectorDouble& center)
  {
    for (Id idim = 0; idim < _ndim; idim++) _center[idim] = center[idim];
  }

  VectorDouble BooleanObject::getValues() const
  {
    VectorDouble tab;
    tab.push_back(static_cast<double>(_mode));
    tab.push_back(static_cast<double>(_token->getType().getValue()));
    for (const auto c: _center) tab.push_back(c);
    for (const auto e: _extension) tab.push_back(e);
    tab.push_back(_orientation);
    return tab;
  }

  void BooleanObject::_defineBoundingBox(double eps)
  {
    double dx, dy, dz;

    if (ABS(_orientation) < eps)
    {
      dx = _extension[0];
      dy = _extension[1];
      dz = _extension[2];
    }
    else
    {
      double angle = _orientation * GV_PI / 180.;
      double sint = ABS(sin(angle));
      double cost = ABS(cos(angle));
      dx = cost * _extension[0] + sint * _extension[1];
      dy = sint * _extension[0] + cost * _extension[1];
      dz = _extension[2];
    }

    _box[0][0] = _center[0] - dx / 2;
    _box[0][1] = _center[0] + dx / 2;
    _box[1][0] = _center[1] - dy / 2;
    _box[1][1] = _center[1] + dy / 2;

    if (!_token->getFlagCutZ())
    {
      _box[2][0] = _center[2] - dz / 2;
      _box[2][1] = _center[2] + dz / 2;
    }
    else
    {
      _box[2][0] = _center[2] - dz;
      _box[2][1] = _center[2];
    }
  }

  void BooleanObject::_generateCoordinatesInPlace(
    const DbGrid* dbout,
    const SimuBooleanParam& boolparam,
    VectorDouble& coor)
  {
    Id ndim = static_cast<Id>(coor.size());
    for (Id idim = 0; idim < ndim; idim++)
    {
      double origin = dbout->getX0(idim);
      origin -= boolparam.getDilate(idim);
      double field = dbout->getExtend(idim);
      field += 2. * boolparam.getDilate(idim);
      coor[idim] = origin + field * law_uniform(0., 1.);
    }
  }

  /*****************************************************************************/
  /*!
   **  Function used to generate the geometry of an object
   **
   *****************************************************************************/
  BooleanObject* BooleanObject::generate(
    const DbGrid* dbout,
    const VectorDouble& cdgrain,
    const ModelBoolean* tokens,
    const SimuBooleanParam& boolparam,
    double eps)
  {
    Id ndim = dbout->getNDim();

    // Define the (primary) location of the object
    Id iter = 0;
    VectorDouble coor(ndim);
    if (!cdgrain.empty())
    {
      coor = cdgrain;
    }
    else
    {
      do
      {
        iter++;
        if (iter > boolparam.getMaxiter()) return nullptr;
        _generateCoordinatesInPlace(dbout, boolparam, coor);
      } while (_invalidTokenFromIntensity(dbout, tokens, coor));
    }

    // Generate an object of the correct Token type
    BooleanObject* object = tokens->generateObject(ndim);

    /* Operate the linkage */
    object->_extensionLinkage();

    // Store the coordinates of the object center
    double valrand;
    if (!cdgrain.empty())
    {
      do
      {
        iter++;
        if (iter > boolparam.getMaxiter())
        {
          delete object;
          return nullptr;
        }
        for (Id idim = 0; idim < ndim; idim++)
        {
          if (idim < 2)
          {
            valrand = law_uniform(0., 1.) - 0.5;
          }
          else
          {
            if (!object->getToken()->getFlagCutZ())
              valrand = law_uniform(0., 1.) - 0.5;
            else
              valrand = law_uniform(0., 1.);
          }
          object->setCenter(
            idim, coor[idim] + object->getExtension(idim) * valrand);
        }
      } while (!object->_isInObject(cdgrain));
    }
    else
      object->setCenter(coor);

    /* Determine the inclusive box */

    object->_defineBoundingBox(eps);

    return object;
  }

  /*****************************************************************************/
  /*!
   **  Function to link the geometries of an object
   **
   *****************************************************************************/
  void BooleanObject::_extensionLinkage()
  {
    if (_token->getFactorX2Y() > 0.)
      _extension[1] = _extension[0] * _token->getFactorX2Y();
    if (_token->getFactorX2Z() > 0.)
      _extension[2] = _extension[0] * _token->getFactorX2Z();
    if (_token->getFactorY2Z() > 0.)
      _extension[2] = _extension[1] * _token->getFactorY2Z();
  }

  /*****************************************************************************/
  /*!
   **  Check if an object may be generated according to the value
   **  of the Intensity
   **  This Intensity can be local or not (if Flag_stat)
   **
   ** \return  Error return code:
   ** \return  FALSE the token is created
   ** \return  TRUE  the token may not be created
   **
   *****************************************************************************/
  bool BooleanObject::_invalidTokenFromIntensity(
    const DbGrid* dbout,
    const ModelBoolean* tokens,
    const VectorDouble& coor,
    double eps)
  {
    double theta;
    if (tokens->isFlagStat())
    {
      theta = tokens->getThetaCst();
    }
    else
    {
      Id iech = dbout->coordinateToRank(coor, false, eps);
      theta = dbout->getLocVariable(ELoc::P, iech, 0);
    }
    return (law_uniform(0., 1.) > theta);
  }

  /*****************************************************************************/
  /*!
   **  Check if the pixel (coor) belongs to the grain.
   **
   ** \return  True if the pixel is in the object, False otherwise
   **
   ** \param[in]  coor     location of the pixel
   **
   *****************************************************************************/
  bool BooleanObject::_isInObject(const VectorDouble& coor)
  {
    VectorDouble incr(_ndim);
    for (Id idim = 0; idim < _ndim; idim++)
      incr[idim] = coor[idim] - _center[idim];

    if (_orientation)
    {
      double angle = _orientation * GV_PI / 180.;
      double sint = sin(angle);
      double cost = cos(angle);
      double dxr = incr[0] * cost + incr[1] * sint;
      double dyr = incr[0] * cost - incr[1] * sint;
      incr[0] = dxr;
      incr[1] = dyr;
    }

    /* Check if the grain is outside the box */

    if (ABS(incr[0]) > _extension[0] / 2.) return false;
    if (ABS(incr[1]) > _extension[1] / 2.) return false;
    if (_ndim > 2)
    {
      if (!_token->getFlagCutZ())
      {
        if (ABS(incr[2]) > _extension[2] / 2.) return false;
      }
      else
      {
        if (incr[2] > 0) return false;
        if (ABS(incr[2]) > _extension[2]) return false;
      }
    }

    /* Check the pixel according to the grain definition */

    bool answer = _token->belongObject(incr, this);

    return answer;
  }

  /*****************************************************************************/
  /*!
   **  Check if the pixel (coor) belongs to the object bounding box
   **
   **  \return True if the pixel is in the box and False otherwise
   **
   ** \param[in]  coor     location of the pixel
   **
   *****************************************************************************/
  bool BooleanObject::_isInBoundingBox(const VectorDouble& coor)

  {
    for (Id idim = 0; idim < _ndim; idim++)
    {
      if (coor[idim] < _box[idim][0]) return false;
      if (coor[idim] > _box[idim][1]) return false;
    }
    return true;
  }

  bool BooleanObject::isPore(const Db* db, Id iech)
  {
    return (db->getZVariable(iech, 0) == 0);
  }

  bool BooleanObject::isGrain(const Db* db, Id iech)
  {
    return (db->getZVariable(iech, 0) != 0);
  }

  Id BooleanObject::_getCoverageAtSample(const Db* db, Id iptr_cover, Id iech)
  {
    return static_cast<Id>(db->getArray(iech, iptr_cover));
  }

  void BooleanObject::_updateCoverageAtSample(
    Db* db,
    Id iptr_cover,
    Id iech,
    Id ival)
  {
    db->setArray(iech, iptr_cover, db->getArray(iech, iptr_cover) + ival);
  }

  /*****************************************************************************/
  /*!
   **  Check if the current object is compatible with the constraining pores
   **
   ** \return True if it is compatible; False otherwise
   **
   ** \param[in]  db  Constraining data set
   **
   *****************************************************************************/
  bool BooleanObject::isCompatiblePore(const Db* db)
  {
    if (db == nullptr) return true;
    for (Id iech = 0; iech < db->getNSample(); iech++)
    {
      // Discard if the data is maked off or corresponds to a Grain
      if (!db->isActive(iech)) continue;
      if (!isPore(db, iech)) continue;

      // The data is a Pore
      VectorDouble coor = db->getSampleCoordinates(iech);

      // Check if the data is in the bounding box of the object
      if (!_isInBoundingBox(coor)) continue;

      // Check if the data is in the object
      if (!_isInObject(coor)) continue;

      // The current object is refused because it contains a constraining pore
      return false;
    }
    return true;
  }

  /*****************************************************************************/
  /*!
   **  Check if an object can be added with regards to the constraining grains
   **
   ** \param[in]  db  Constraining data set
   **
   *****************************************************************************/
  bool BooleanObject::isCompatibleGrainAdd(const Db* db)
  {
    if (db == nullptr) return true;

    for (Id iech = 0; iech < db->getNSample(); iech++)
    {
      if (!db->isActive(iech)) continue;
      if (!isGrain(db, iech)) continue;
      VectorDouble coor = db->getSampleCoordinates(iech);
      if (!_isInBoundingBox(coor)) continue;
      if (!_isInObject(coor)) continue;
    }
    return true;
  }

  /*****************************************************************************/
  /*!
   **  Check if an object can be deleted with regards to the constraining grains
   **
   ** \param[in]  db         Constraining data set
   ** \param[in]  iptr_cover UIUD for coverage variable
   **
   *****************************************************************************/
  bool BooleanObject::isCompatibleGrainDelete(const Db* db, Id iptr_cover)
  {
    if (db == nullptr) return true;
    for (Id iech = 0; iech < db->getNSample(); iech++)
    {
      if (!db->isActive(iech)) continue;
      if (!isGrain(db, iech)) continue;
      VectorDouble coor = db->getSampleCoordinates(iech);
      if (!_isInBoundingBox(coor)) continue;
      if (_getCoverageAtSample(db, iptr_cover, iech) > 1) continue;
      if (_isInObject(coor)) return false;
    }
    return true;
  }

  /*****************************************************************************/
  /*!
   **  Update the covering value of each constraining grain after a
   **  deletion or an addition operation
   **
   ** \return  Count of grains not covered after the operation
   **
   ** \param[in]  db       Db structure
   ** \param[in]  iptr_cover UID for the covering variable
   ** \param[in]  val      type of the operation to be tested
   **                      1 for addition; -1 for deletion
   **
   *****************************************************************************/
  Id BooleanObject::coverageUpdate(Db* db, Id iptr_cover, Id val)
  {
    if (db == nullptr) return 0;
    Id not_covered = 0;
    for (Id iech = 0; iech < db->getNSample(); iech++)
    {
      if (!db->isActive(iech)) continue;
      if (!isGrain(db, iech)) continue;
      VectorDouble coor = db->getSampleCoordinates(iech);
      if (_isInBoundingBox(coor))
      {
        if (_isInObject(coor))
        {
          _updateCoverageAtSample(db, iptr_cover, iech, val);
        }
      }
      if (_getCoverageAtSample(db, iptr_cover, iech) <= 0) not_covered++;
    }
    return not_covered;
  }

  void BooleanObject::projectToGrid(
    DbGrid* dbout,
    Id iptr_simu,
    Id iptr_rank,
    Id facies,
    Id rank)
  {
    Id ix0, ix1, iy0, iy1, iz0, iz1;
    VectorDouble coor(_ndim);
    VectorInt indice(_ndim);

    /* Look for the nodes in the box of influence of the object */

    if (_ndim >= 1)
    {
      ix0 =
        static_cast<Id>((_box[0][0] - dbout->getX0(0)) / dbout->getDX(0) - 1);
      ix0 = MAX(ix0, 0);
      ix1 =
        static_cast<Id>((_box[0][1] - dbout->getX0(0)) / dbout->getDX(0) + 1);
      ix1 = MIN(ix1, dbout->getNX(0) - 1);
    }
    else
    {
      ix0 = 0;
      ix1 = 0;
    }
    if (_ndim >= 2)
    {
      iy0 =
        static_cast<Id>((_box[1][0] - dbout->getX0(1)) / dbout->getDX(1) - 1);
      iy0 = MAX(iy0, 0);
      iy1 =
        static_cast<Id>((_box[1][1] - dbout->getX0(1)) / dbout->getDX(1) + 1);
      iy1 = MIN(iy1, dbout->getNX(1) - 1);
    }
    else
    {
      iy0 = 0;
      iy1 = 0;
    }
    if (_ndim >= 3)
    {
      iz0 =
        static_cast<Id>((_box[2][0] - dbout->getX0(2)) / dbout->getDX(2) - 1);
      iz0 = MAX(iz0, 0);
      iz1 =
        static_cast<Id>((_box[2][1] - dbout->getX0(2)) / dbout->getDX(2) + 1);
      iz1 = MIN(iz1, dbout->getNX(2) - 1);
    }
    else
    {
      iz0 = 0;
      iz1 = 0;
    }

    /* Check the pixels within the box */

    for (Id ix = ix0; ix <= ix1; ix++)
      for (Id iy = iy0; iy <= iy1; iy++)
        for (Id iz = iz0; iz <= iz1; iz++)
        {
          if (_ndim >= 1) coor[0] = dbout->getX0(0) + ix * dbout->getDX(0);
          if (_ndim >= 2) coor[1] = dbout->getX0(1) + iy * dbout->getDX(1);
          if (_ndim >= 3) coor[2] = dbout->getX0(2) + iz * dbout->getDX(2);

          if (!_isInObject(coor)) continue;

          if (_ndim >= 1) indice[0] = ix;
          if (_ndim >= 2) indice[1] = iy;
          if (_ndim >= 3) indice[2] = iz;

          Id iad = dbout->indiceToRank(indice);

          /* Bypass writing if the cell is masked off */

          if (!dbout->isActive(iad)) continue;

          /* Set the values */

          if (iptr_simu >= 0)
          {
            dbout->setArray(iad, iptr_simu, facies);
          }
          if (iptr_rank >= 0)
          {
            double value = dbout->getArray(iad, iptr_rank);
            if (FFFF(value) || value == 0)
              dbout->setArray(iad, iptr_rank, static_cast<double>(rank));
          }
        }
  }
} // namespace gstlrn
