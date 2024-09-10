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
#include "Covariances/ANoStatCov.hpp"
#include "geoslib_define.h"
#include "Covariances/ACov.hpp"

#include <math.h>

ANoStatCov::ANoStatCov()
    : AStringable(),
      _items(),
      _amesh(nullptr),
      _dbin(nullptr),
      _dbout(nullptr)
{
}

ANoStatCov::ANoStatCov(const VectorString& codes)
  : AStringable(),
    _items(),
    _amesh(nullptr),
    _dbin(nullptr),
    _dbout(nullptr)
{
  (void) addNoStatElems(codes);
}

ANoStatCov::ANoStatCov(const ANoStatCov &m)
  : AStringable(m),
    _items(m._items),
    _amesh(m._amesh),
    _dbin(m._dbin),
    _dbout(m._dbout)
{
}

ANoStatCov& ANoStatCov::operator= (const ANoStatCov &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _items = m._items;
    _amesh = m._amesh;
    _dbin  = m._dbin;
    _dbout = m._dbout;
  }
  return *this;
}

ANoStatCov::~ANoStatCov()
{
}

String ANoStatCov::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  if (getNoStatElemNumber() <= 0) return sstr.str();

  sstr << toTitle(1, "Non-Stationary Parameters");
  for (int i = 0; i < (int) getNoStatElemNumber(); i++)
    sstr << std::to_string(i+1) << " - " << _items[i].toString(strfmt);
  return sstr.str();
}

/**
 * Look if a Non-stationary parameter is defined
 * @return
 */
bool ANoStatCov::isDefinedByCov() const
{
  if (_items.empty()) return false;
  return !_items.empty();
}

/**
 * Look if a Non-stationary parameter is defined
 * @param type Rank of Target Type (or EConsElem::UNKNOWN for any)
 * @return
 */
bool ANoStatCov::isDefinedByType(const EConsElem& type) const
{
  if (_items.empty()) return false;
  for (int ipar = 0; ipar < (int) getNoStatElemNumber(); ipar++)
  {
      if (matchType(ipar,type)) return true;
  }
  return false;
}

/**
 * Look if a Non-stationary parameter is defined
 * @param type Rank of Target Type (or EConsElem::UNKNOWN for any)
 * @return
 */
bool ANoStatCov::isDefinedByCovType(const EConsElem& type) const
{
  if (_items.empty()) return false;
  for (int ipar = 0; ipar < (int) getNoStatElemNumber(); ipar++)
  {
    if (matchType(ipar,type)) return true;
  }
  return false;
}

bool ANoStatCov::isDefined(const EConsElem& type,
                           int iv1,
                           int iv2) const
{
  if (_items.empty()) return false;
  for (int ipar = 0; ipar < (int) getNoStatElemNumber(); ipar++)
  {
    if (matchType(ipar,type) &&
        matchIV1 (ipar,iv1)  &&
        matchIV2 (ipar,iv2)) return true;
  }
  return false;
}

/**
 * Look if a Non-stationary parameter for Variance is defined
 * for any pairs of variables
 * @return
 */
bool ANoStatCov::isDefinedForVariance() const
{
  if (_items.empty()) return false;
  for (int ipar = 0; ipar < (int) getNoStatElemNumber(); ipar++)
  {
    if (getType(ipar) == EConsElem::SILL) return true;
  }
  return false;
}

/**
 * Look if a Non-stationary parameter for Anisotropy is defined
 * either by Tensor or by (angle/range/scale)
 * @return
 */
bool ANoStatCov::isDefinedforAnisotropy() const
{
  if (_items.empty()) return false;
  for (int ipar = 0; ipar < (int) getNoStatElemNumber(); ipar++)
  {
    if (getType(ipar) == EConsElem::ANGLE ||
        getType(ipar) == EConsElem::RANGE ||
        getType(ipar) == EConsElem::SCALE ||
        getType(ipar) == EConsElem::TENSOR) return true;
  }
  return false;
}

/**
 * Look if a Non-stationary parameter for Anisotropy is defined
 * only by angle / range / scale
 * @return
 */
bool ANoStatCov::isDefinedforRotation() const
{
  if (_items.empty()) return false;
  for (int ipar = 0; ipar < (int) getNoStatElemNumber(); ipar++)
  {
    if (getType(ipar) == EConsElem::ANGLE ||
        getType(ipar) == EConsElem::RANGE ||
        getType(ipar) == EConsElem::SCALE) return true;
  }
  return false;
}

/**
 * Return the rank for a Non-stationary parameter
 * @param type Rank of Target Type (or EConsElem::UNKNOWN for any)
 * @param iv1  Rank of the first additional element designation (or -1 for any)
 * @param iv2  Rank of the second additional element designation (or -1 for any)
 * @return -1 if no match is found
 */
int ANoStatCov::getRank(const EConsElem& type, int iv1, int iv2) const
{
  if (_items.empty()) return -1;
  for (int ipar = 0; ipar < (int) getNoStatElemNumber(); ipar++)
  {
    if (matchType(ipar,type) &&
        matchIV1 (ipar,iv1)  &&
        matchIV2 (ipar,iv2)) return ipar;
  }
  return -1;
}

/**
 * Add an element to the pile of Non-stationary parameters
 * @param type Type of the non-stationary parameter
 * @param iv1  Rank of the first designation
 * @param iv2  Rank of the second designation
 */
int ANoStatCov::addNoStatElem(const EConsElem& type, int iv1, int iv2)
{
  int nelem = static_cast<int> (_items.size());
  _items.resize(nelem+1);
  _items[nelem].init(type, iv1, iv2);
  if (! _checkConsistency())
  {
    deleteNoStatElem(nelem);
    return 1;
  }
  return 0;
}

int ANoStatCov::addNoStatElemByItem(const ParamId& item)
{
  return addNoStatElem( item.getType(),
                       item.getIV1(), item.getIV2());
}

int ANoStatCov::addNoStatElems(const VectorString &codes)
{
  int iv1, iv2;
  EConsElem type;

  int nerr = 0;
  for (int i = 0; i < (int) codes.size(); i++)
  {
    if (_understandCode(codes[i],&type, &iv1, &iv2)) continue;
    nerr += addNoStatElem( type, iv1, iv2);
  }
  return (nerr > 0);
}

void ANoStatCov::deleteNoStatElem(int ipar)
{
  _items.erase(_items.begin() + ipar);
}

void ANoStatCov::deleteAllNoStatElem()
{
  _items.clear();
}

/**
 * This function tries to attach the current Non-Stationary environment
 * to the cova, checking that references to Covariance ranks are valid
 * @param cova ACov structure
 * @return Error return code
 */
int ANoStatCov::attachCov(const ACov* cova)
{
  if (cova == nullptr)
  {
    messerr("CovAniso must be defined beforehand");
    return 1;
  }
  if (cova->getNDim() > 3)
  {
    messerr("Non stationary cova is restricted to Space Dimension <= 3");
    return 1;
  }

  // Patch. It may happen that the CovAniso already contains the parameters of the ANostat
  // which are better defined than in the Current ANostat structure
  // In this case copy the ANostat parameters from CovAniso to Current
  // TODO: Remove this part of code
  _updateFromCova(cova);

  for (int ipar=0; ipar<(int) getNoStatElemNumber(); ipar++)
  {
    EConsElem type = getType(ipar);

    // Check that the Non-stationary parameter is valid with respect
    // to the CovAniso definition

    if (type == EConsElem::PARAM)
    {
      messerr("The current methodology does not handle constraint on third parameter");
      return 1;
    }
  }
  return 0;
}

/** Decoding a NoStat Element code which should follow the following syntax
 *
 * @param code Input code
 * @param type Type of parameter ("R","A","P","V","S","T","C","I","H")
 *             G : Rank of the GRF
 *             M : Rank of the Structure
 *             R : Range
 *             A : Angle
 *             P : Third parameter
 *             V : Sill
 *             S : Scale
 *             T : Tapering range
 *             C : Velocity (advection)
 *             I : Rotation angle for Sphere
 *             H : Anisotropy matrix terms
 * @param iv1  Rank of the first variable
 * @param iv2  Rank of the second variable
 * @return 0 if a valid constraint has been defined; 1 otherwise
 */
int ANoStatCov::_understandCode(const String& code,
                                EConsElem *type,
                                int *iv1,
                                int *iv2)
{
  *iv1 = *iv2 = 0;
  *type = EConsElem::UNKNOWN;
  VectorString keywords = separateKeywords(code);

  bool flagV1  = false;
  bool flagV2  = false;
  int size = static_cast<int> (keywords.size());
  int lec = 0;

  // Decoding the operator keyword ("R","A","P","V","S","T","C","I","H")
  if (lec < size)
  {
    // TODO better decoding which would not be dependent on the order in ENUM_CONS
    VectorString list = {"R","A","P","V","S","T","C","I","H"};
    int itype = getRankInList(list, keywords[lec], false) + 1;
    (*type) = EConsElem::fromValue(itype);
    if ((*type) == EConsElem::UNKNOWN) return 1;
    flagV1 = true;
    lec++;
  }

  // Decoding the rank of the operator (conditional)
  if (lec < size && flagV1)
  {
    *iv1 = toInteger(keywords[lec]);
    if (IFFFF(*iv1)) return 1;
    (*iv1)--;
    lec++;
  }

  // Decoding the separator (keyword "-")
  if (lec < size)
  {
    if (matchKeyword(keywords[lec],"-",false))
    {
      flagV2 = true;
      lec++;
    }
  }

  // Decoding the rank of the second variable (conditional)
  if (lec < size)
  {
    if (flagV2)
    {
      *iv2 = toInteger(keywords[lec]);
      if (IFFFF(*iv2)) return 1;
      (*iv2)--;
      lec++;
    }
    else
    {
      messerr("Wrong character ('%s') found in code %s", keywords[lec].c_str(), code.c_str());
      return 1;
    }
  }
  return 0;
}

/**
 * The CovAniso may contain some non-stationary parameters already initialized
 * that the current structure is not aware of. They are simply duplicated here
 * This temporary patch is applied in the following conditions:
 * 1) the number of stationary elements in current structure is 0
 * 2) the number of stationary elements in CovAniso is positive
 * @param cova Input ACov
 */
void ANoStatCov::_updateFromCova(const ACov* cova)
{
  int nelemFromCovAniso = cova->getNoStatElemNumber();
  int nelem = getNoStatElemNumber();
  if (nelem > 0) return;
  if (nelemFromCovAniso <= 0) return;

  for (int ipar = 0; ipar < nelemFromCovAniso; ipar++)
  {
    ParamId item = cova->getCovParamId(ipar);
    (void) addNoStatElemByItem(item);
  }
}

/**
 * Get the information from the storage in Dbin and/or Dbout
 * @param ipar  Rank of the non-stationary parameter
 * @param icas1 Type of first Db: 1 for Input; 2 for Output
 * @param iech1 Rank of the first sample (in Dbin)
 * @param icas2 Type of first Db: 1 for Input; 2 for Output
 * @param iech2 Rank of the second sample (in Dbout)
 * @param val1  Returned value at first sample
 * @param val2  Returned value at the second sample
 *
 * @return true if the two values are different
 */
bool ANoStatCov::getInfoFromDb(int ipar,
                            int icas1,
                            int iech1,
                            int icas2,
                            int iech2,
                            double *val1,
                            double *val2) const
{
  *val1 = getValueByParam(ipar, icas1, iech1);
  *val2 = getValueByParam(ipar, icas2, iech2);

  if (FFFF(*val1) && FFFF(*val2)) return false;

  if (FFFF(*val1)) *val1 = *val2;
  if (FFFF(*val2)) *val2 = *val1;

  return *val1 != *val2;
}

int ANoStatCov::attachToMesh(const AMesh* mesh, bool center, bool /*verbose*/) const
{
  DECLARE_UNUSED_(center);
  _setAmesh(mesh);
  return 0;
}

int ANoStatCov::attachToDb(Db* db, int icas, bool /*verbose*/) const
{
  if (icas == 1)
    _setDbin(db);
  else
    _setDbout(db);
  return 0;
}

void ANoStatCov::detachFromMesh() const
{
  _setAmesh(nullptr);
}

void ANoStatCov::detachFromDb(Db* /*db*/, int icas) const
{
  if (icas == 1)
    _setDbin(nullptr);
  else
    _setDbout(nullptr);
}

/**
 * This function is meant to check the consistency of the different
 * non-stationary parameters, i.e:
 * - HH is incompatible with (angle / range/ scale)
 * @return
 */
bool ANoStatCov::_checkConsistency() const
{
  bool flagHH  = false;
  bool flagRot = false;
  int nitem = getNoStatElemNumber();

  for (int ipar = 0; ipar < nitem; ipar++)
  {
    flagHH = flagHH   ||
        matchType(ipar, EConsElem::TENSOR);
    flagRot = flagRot ||
        matchType(ipar, EConsElem::ANGLE) ||
        matchType(ipar, EConsElem::RANGE) ||
        matchType(ipar, EConsElem::SCALE);
  }
  if (flagHH && flagRot)
  {
    messerr("Error in the definition of CovAniso Non-stationarity");
    messerr("You cannot mix the following two parameterizations:");
    messerr("- in Tensor using HH");
    messerr("- in rotation matrix using Angle / [Scale | Range]");
    return false;
  }
  return true;
}

/**
 * This (temporary) function checks the validity between arguments 'icas' and 'rank'
 * @param icas  Source definition:
 *              0 : from Meshing (rank: absolute rank to be converted into relative)
 *              1 : from Dbin
 *              2 : from Dbout
 * @param rank  Rank of the target
 * @return
 */
bool ANoStatCov::_isValid(int icas, int rank) const
{
  switch (icas)
  {
    case 0:
      if (_amesh == nullptr)
      {
        messerr("Checking the validity of the argument");
        messerr("Meshing: This requires '_amesh' to be defined beforehand");
        return false;
      }
      if (rank < 0 || rank >= _amesh->getNMeshes())
      {
        messerr("Check the validity of the argument");
        messerr("Meshing: 'rank' (%d) should be smaller than number of meshes (%d)",
                rank, _amesh->getNMeshes());
        return false;
      }
      break;

    case 1:
      if (_dbin == nullptr)
      {
        messerr("Checking the validity of the argument");
        messerr("Dbin: This requires '_dbin' to be defined beforehand");
        return false;
      }
      if (rank < 0 || rank >= _dbin->getSampleNumber(0))
      {
        messerr("Check the validity of the argument");
        messerr("Dbin: 'rank' (%d) should be smaller than number of samples (%d)",
                rank, _dbin->getSampleNumber(0));
        return false;
      }
      break;

    case 2:
      const Db* dbout = _getDbout();
      if (dbout == nullptr)
      {
        messerr("Checking the validity of the argument");
        messerr("Dbout: This requires '_dbout' to be defined beforehand");
        return false;
      }
      if (rank < 0 || rank >= dbout->getSampleNumber(0))
      {
        messerr("Check the validity of the argument");
        messerr("Dbout: 'rank' (%d) should be smaller than number of samples (%d)",
                rank, dbout->getSampleNumber(0));
        return false;
      }
      break;
  }
  return true;
}

/*****************************************************************************/
/*!
 **  Derive the non-stationary information(s) from the Output db (if Grid)
 **  to the Input Db
 **
 ** \return  Error return code
 **
 ** \param[in]  mode        1 for allocation; -1 for deallocation
 ** \param[in]  dbin        Descriptor of the input Db
 ** \param[in]  dbout       Descriptor of the output Db
 **
 *****************************************************************************/
int ANoStatCov::manageInfo(int mode, Db *dbin, Db *dbout) const
{

  /* Dispatch */

  if (mode > 0)
  {
    if (dbin != nullptr)
    {
      // Attach the Input Db
      if (attachToDb(dbin, 1)) return 1;
    }

    if (dbout != nullptr)
    {
      // Attach the Output Db
      if (attachToDb(dbout, 2)) return 1;
    }
  }
  else
  {
    if (dbin != nullptr)
    {
      // Detach the Input Db
      detachFromDb(dbin, 1);
    }

    if (dbout != nullptr)
    {
      // Detach the output Db
      detachFromDb(dbout, 2);
    }
  }
  return (0);
}

void ANoStatCov::checkCode(const String& code)
{
  int iv1, iv2;
  EConsElem type;

  if (ANoStatCov::_understandCode(code, &type, &iv1, &iv2)) return;

  message("Code '%s' decodes as follows (1-based)\n", code.c_str());
  message("TYPE = %s\n", type.getDescr().c_str());
  message("IV1  = %d\n", iv1+1);
  message("IV2  = %d\n", iv2+1);
}

const Db* ANoStatCov::_getDbout() const
{
  if (_dbout == nullptr && _dbin != nullptr) return _dbin;
  return _dbout;
}
