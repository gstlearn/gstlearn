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
#include "Boolean/AShape.hpp"
#include "Boolean/ModelBoolean.hpp"
#include "Simulation/SimuBooleanParam.hpp"
#include "Simulation/BooleanObject.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Law.hpp"

#include <math.h>

BooleanObject::BooleanObject(const AShape* ashape)
    : AStringable(),
      _mode(0),
      _token(ashape),
      _center({0.,0.,0.}),
      _extension({0.,0.,0.,}),
      _orientation(0.),
      _values({0.,0.,0.}),
      _box({{0.,0.},{0.,0.},{0.,0.}})
{
}

BooleanObject::BooleanObject(const BooleanObject &r)
    : AStringable(r),
      _mode(r._mode),
      _token(r._token),
      _center(r._center),
      _extension(r._extension),
      _orientation(r._orientation),
      _values(r._values),
      _box(r._box)
{
}

BooleanObject& BooleanObject::operator=(const BooleanObject &r)
{
  if (this != &r)
  {
    AStringable::operator =(r);
    _mode = r._mode;
    _token = r._token;
    _center = r._center;
    _extension = r._extension;
    _orientation = r._orientation;
    _values = r._values;
    _box = r._box;
  }
  return *this;
}

BooleanObject::~BooleanObject()
{
}

String BooleanObject::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  if (_mode == 1)
    sstr << "Primary Object" << std::endl;
  else
    sstr << "Secondary Object" << std::endl;
  sstr << "- Type        = " << _token->getType().getDescr() << std::endl;
  sstr << "- Center      = " << VH::toString(_center);
  sstr << "- Extension   = " << VH::toString(_extension);
  sstr << "- Orientation = " << _orientation << std::endl;

  return sstr.str();
}

void BooleanObject::setCenter(const VectorDouble& center)
{
  _center.resize(3,0.);
  for (int idim = 0; idim < (int) center.size(); idim++)
    _center[idim] = center[idim];
}

VectorDouble BooleanObject::getValues() const
{
  VectorDouble tab;
  tab.push_back((double) _mode);
  tab.push_back((double) _token->getType().getValue());
  tab.insert(tab.end(), _center.begin(), _center.end());
  tab.insert(tab.end(), _extension.begin(), _extension.end());
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

  if (_token->getFlagCutZ())
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

void BooleanObject::_drawCoordinate(const DbGrid *dbout,
                             const SimuBooleanParam& boolparam,
                             VectorDouble& coor)
{
  int ndim = dbout->getNDim();
  for (int idim = 0; idim < ndim; idim++)
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
BooleanObject* BooleanObject::generate(const DbGrid* dbout,
                         const VectorDouble& cdgrain,
                         const ModelBoolean* tokens,
                         const SimuBooleanParam& boolparam,
                         double eps)
{
  int ndim = dbout->getNDim();

  // Define the (primary) location of the object

  int iter = 0;
  VectorDouble coor(ndim);
  if (! cdgrain.empty())
  {
    coor = cdgrain;
  }
  else
  {
    do
    {
      iter++;
      if (iter > boolparam.getMaxiter()) return nullptr;
      _drawCoordinate(dbout, boolparam, coor);
    }
    while (_checkIntensity(dbout, tokens, coor));
  }

  // Generate an object of the correct Token type

  BooleanObject* object = tokens->generateObject(ndim);

  /* Operate the linkage */

  object->_extensionLinkage();

  /* Store the coordinates of the object center */

  double valrand;
  if (! cdgrain.empty())
  {
    do
    {
      iter++;
      if (iter > boolparam.getMaxiter())
      {
        delete object;
        return nullptr;
      }
      for (int idim = 0; idim < ndim; idim++)
      {
        if (idim < 2)
        {
          valrand = law_uniform(0., 1.) - 0.5;
        }
        else
        {
          if (object->getToken()->getFlagCutZ())
            valrand = law_uniform(0., 1.) - 0.5;
          else
            valrand = law_uniform(0., 1.);
        }
        object->setCenter(idim, coor[idim] +
                          object->getExtension(idim) * valrand);
      }
    }
    while (! object->_checkObject(cdgrain, ndim));
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
  ** \return  0 the token is created
  ** \return  1 the token may not be created
  **
  *****************************************************************************/
bool BooleanObject::_checkIntensity(const DbGrid* dbout,
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
     int iech = dbout->coordinateToRank(coor, false, eps);
     theta = dbout->getLocVariable(ELoc::P,iech, 0);
   }
   return (law_uniform(0., 1.) > theta);
 }

/*****************************************************************************/
/*!
 **  Check if the pixel (coor) belongs to the grain.
 **
 ** \return  true if the pixel is in the grain, 0false otherwise
 **
 ** \param[in]  coor     location of the pixel
 ** \param[in]  ndim     Space dimension
 **
 *****************************************************************************/
bool BooleanObject::_checkObject(const VectorDouble& coor, int ndim)
{
  VectorDouble incr(ndim);
  for (int idim = 0; idim < ndim; idim++)
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

  if (ndim > 2)
  {
    if (_token->getFlagCutZ())
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
 ** \param[in]  coor     location of the pixel
 ** \param[in]  ndim     Space dimension
 **
 *****************************************************************************/
bool BooleanObject::_checkBoundingBox(const VectorDouble& coor, int ndim)

{
  for (int idim = 0; idim < ndim; idim++)
  {
    if (coor[idim] < _box[idim][0]) return false;
    if (coor[idim] > _box[idim][1]) return false;
  }
  return true;
}

bool BooleanObject::_isPore(const Db* db, int iech)
{
  return (db->getLocVariable(ELoc::Z,iech, 0) == 0);
}

bool BooleanObject::_isGrain(const Db* db, int iech)
{
  return (db->getLocVariable(ELoc::Z,iech, 0) != 0);
}

int BooleanObject::_getCoverageAtSample(const Db* db, int iech)
{
  return (int) db->getLocVariable(ELoc::Z,iech,1);
}

void BooleanObject::_updateCoverageAtSample(Db* db, int iech, int ival)
{
  db->setLocVariable(ELoc::Z,iech, 1, db->getLocVariable(ELoc::Z,iech, 1) + ival);
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
  int ndim = db->getNDim();
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (! db->isActive(iech)) continue;
    if (! _isPore(db, iech)) continue;
    VectorDouble coor = db->getSampleCoordinates(iech);
    if (! _checkBoundingBox(coor, ndim)) continue;
    if (_checkObject(coor, ndim)) return true;
  }
  return false;
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

  int ndim = db->getNDim();
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (! db->isActive(iech)) continue;
    if (! _isGrain(db,iech)) continue;
    VectorDouble coor = db->getSampleCoordinates(iech);
    if (! _checkBoundingBox(coor, ndim)) continue;
    if (! _checkObject(coor, ndim)) continue;
  }
  return true;
}

/*****************************************************************************/
/*!
 **  Check if an object can be deleted with regards to the constraining grains
 **
 ** \param[in]  db       Constraining data set
 **
 *****************************************************************************/
bool BooleanObject::isCompatibleGrainDelete(const Db* db)
{
  if (db == nullptr) return true;
  int ndim = db->getNDim();

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (! db->isActive(iech)) continue;
    if (! _isGrain(db, iech)) continue;
    VectorDouble coor = db->getSampleCoordinates(iech);
    if (! _checkBoundingBox(coor, ndim)) continue;
    if (_getCoverageAtSample(db, iech) > 1) continue;
    if (_checkObject(coor, ndim)) return false;
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
 ** \param[in]  val      type of the operation to be tested
 **                      1 for addition; -1 for deletion
 **
 *****************************************************************************/
int BooleanObject::coverageUpdate(Db* db, int val)
{
  if (db == nullptr) return 0;
  int ndim = db->getNDim();
  int not_covered = 0;
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (! db->isActive(iech)) continue;
    if (! _isGrain(db, iech)) continue;
    VectorDouble coor = db->getSampleCoordinates(iech);
    if (_checkBoundingBox(coor, ndim))
    {
      if (_checkObject(coor, ndim))
      {
        if (val < 0)
        {
          // Deletion
          _updateCoverageAtSample(db, iech, -1);
        }
        else
        {
          // Addition
          _updateCoverageAtSample(db, iech, +1);
        }
      }
    }
    if (_getCoverageAtSample(db, iech) <= 0) not_covered++;
  }
  return not_covered;
}

void BooleanObject::projectToGrid(DbGrid* dbout,
                           int iptr_simu,
                           int iptr_rank,
                           int facies,
                           int rank)
{
  int ix0, ix1, iy0, iy1, iz0, iz1;
  int ndim = dbout->getNDim();
  VectorDouble coor(ndim);
  VectorInt indice(ndim);

  /* Look for the nodes in the box of influence of the object */

  if (ndim >= 1)
  {
    ix0 = (int) ((_box[0][0] - dbout->getX0(0)) / dbout->getDX(0) - 1);
    ix0 = MAX(ix0, 0);
    ix1 = (int) ((_box[0][1] - dbout->getX0(0)) / dbout->getDX(0) + 1);
    ix1 = MIN(ix1, dbout->getNX(0) - 1);
  }
  else
  {
    ix0 = 0;
    ix1 = 0;
  }
  if (ndim >= 2)
  {
    iy0 = (int) ((_box[1][0] - dbout->getX0(1)) / dbout->getDX(1) - 1);
    iy0 = MAX(iy0, 0);
    iy1 = (int) ((_box[1][1] - dbout->getX0(1)) / dbout->getDX(1) + 1);
    iy1 = MIN(iy1, dbout->getNX(1) - 1);
  }
  else
  {
    iy0 = 0;
    iy1 = 0;
  }
  if (ndim >= 3)
  {
    iz0 = (int) ((_box[2][0] - dbout->getX0(2)) / dbout->getDX(2) - 1);
    iz0 = MAX(iz0, 0);
    iz1 = (int) ((_box[2][1] - dbout->getX0(2)) / dbout->getDX(2) + 1);
    iz1 = MIN(iz1, dbout->getNX(2) - 1);
  }
  else
  {
    iz0 = 0;
    iz1 = 0;
  }

  /* Check the pixels within the box */

  for (int ix = ix0; ix <= ix1; ix++)
    for (int iy = iy0; iy <= iy1; iy++)
      for (int iz = iz0; iz <= iz1; iz++)
      {
        if (ndim >= 1)
          coor[0] = dbout->getX0(0) + ix * dbout->getDX(0);
        if (ndim >= 2)
          coor[1] = dbout->getX0(1) + iy * dbout->getDX(1);
        if (ndim >= 3)
          coor[2] = dbout->getX0(2) + iz * dbout->getDX(2);

        if (! _checkObject(coor, ndim)) continue;

        if (ndim >= 1) indice[0] = ix;
        if (ndim >= 2) indice[1] = iy;
        if (ndim >= 3) indice[2] = iz;

        int iad = dbout->indiceToRank(indice);

        /* Bypass writing if the cell is masked off */

        if (! dbout->isActive(iad)) continue;

        /* Set the values */

        if (iptr_simu >= 0)
        {
          dbout->setArray(iad, iptr_simu, facies);
        }
        if (iptr_rank >= 0)
        {
          if (FFFF(dbout->getArray(iad, iptr_rank)))
            dbout->setArray(iad, iptr_rank, (double) (rank + 1));
        }
      }
}
