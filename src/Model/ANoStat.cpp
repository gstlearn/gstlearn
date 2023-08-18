/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Model/ANoStat.hpp"
#include "Model/Model.hpp"
#include "Covariances/CovAniso.hpp"
#include "Basic/String.hpp"
#include "Basic/Utilities.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"

#include <math.h>

ANoStat::ANoStat()
    : AStringable(),
      _items(),
      _amesh(nullptr),
      _dbin(nullptr),
      _dbout(nullptr)
{
}

ANoStat::ANoStat(const VectorString& codes)
  : AStringable(),
    _items(),
    _amesh(nullptr),
    _dbin(nullptr),
    _dbout(nullptr)
{
  (void) addNoStatElems(codes);
}

ANoStat::ANoStat(const ANoStat &m)
  : AStringable(m),
    _items(m._items),
    _amesh(m._amesh),
    _dbin(m._dbin),
    _dbout(m._dbout)
{
}

ANoStat& ANoStat::operator= (const ANoStat &m)
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

ANoStat::~ANoStat()
{
}

String ANoStat::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  if (getNoStatElemNumber() <= 0) return sstr.str();

  sstr << toTitle(1, "Non-Stationary Parameters");
  for (int i = 0; i < (int) getNoStatElemNumber(); i++)
    sstr << _items[i].toString(strfmt);
  return sstr.str();
}

/**
 * Look if a Non-stationary parameter is defined
 * @param icov Rank of Target Covariance (or -1 for any)
 * @param igrf Rank of Target GRF (or -1 for any)
 * @return
 */
bool ANoStat::isDefinedByCov(int icov, int igrf) const
{
  if (_items.empty()) return false;
  for (int ipar = 0; ipar < (int) getNoStatElemNumber(); ipar++)
  {
    if (matchIGrf(ipar,igrf) &&
        matchICov(ipar,icov)) return true;
  }
  return false;
}

/**
 * Look if a Non-stationary parameter is defined
 * @param type Rank of Target Type (or EConsElem::UNKNOWN for any)
 * @param igrf Rank of Target GRF (or -1 for any)
 * @return
 */
bool ANoStat::isDefinedByType(const EConsElem& type, int igrf) const
{
  if (_items.empty()) return false;
  for (int ipar = 0; ipar < (int) getNoStatElemNumber(); ipar++)
  {
    if (matchIGrf(ipar,igrf) &&
        matchType(ipar,type)) return true;
  }
  return false;
}

/**
 * Look if a Non-stationary parameter is defined
 * @param type Rank of Target Type (or EConsElem::UNKNOWN for any)
 * @param icov Rank of the Target Covariance (or -1 for any)
 * @param igrf Rank of Target GRF (or -1 for any)
 * @return
 */
bool ANoStat::isDefinedByCovType(const EConsElem& type, int icov, int igrf) const
{
  if (_items.empty()) return false;
  for (int ipar = 0; ipar < (int) getNoStatElemNumber(); ipar++)
  {
    if (matchIGrf(ipar,igrf)  &&
        matchICov(ipar,icov)  &&
        matchType(ipar,type)) return true;
  }
  return false;
}

bool ANoStat::isDefined(const EConsElem& type,
                        int icov,
                        int iv1,
                        int iv2,
                        int igrf) const
{
  if (_items.empty()) return false;
  for (int ipar = 0; ipar < (int) getNoStatElemNumber(); ipar++)
  {
    if (matchIGrf(ipar,igrf) &&
        matchICov(ipar,icov) &&
        matchType(ipar,type) &&
        matchIV1 (ipar,iv1)  &&
        matchIV2 (ipar,iv2)) return true;
  }
  return false;
}

/**
 * Look if a Non-stationary parameter for Anisotropy is defined
 * either by Tensor or by (angle/range/scale)
 * @param icov Rank of the Target Covariance (or -1 for any)
 * @param igrf Rank of Target GRF (or -1 for any)
 * @return
 */
bool ANoStat::isDefinedforAnisotropy(int icov, int igrf) const
{
  if (_items.empty()) return false;
  for (int ipar = 0; ipar < (int) getNoStatElemNumber(); ipar++)
  {
    if (! matchIGrf(ipar,igrf)) continue;
    if (! matchICov(ipar,icov)) continue;
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
 * @param icov Rank of the Target Covariance (or -1 for any)
 * @param igrf Rank of Target GRF (or -1 for any)
 * @return
 */
bool ANoStat::isDefinedforRotation(int icov, int igrf) const
{
  if (_items.empty()) return false;
  for (int ipar = 0; ipar < (int) getNoStatElemNumber(); ipar++)
  {
    if (! matchIGrf(ipar,igrf)) continue;
    if (! matchICov(ipar,icov)) continue;
    if (getType(ipar) == EConsElem::ANGLE ||
        getType(ipar) == EConsElem::RANGE ||
        getType(ipar) == EConsElem::SCALE) return true;
  }
  return false;
}

/**
 * Return the rank for a Non-stationary parameter
 * @param type Rank of Target Type (or EConsElem::UNKNOWN for any)
 * @param icov Rank of the Target Covariance (or -1 for any)
 * @param iv1  Rank of the first additional element designation (or -1 for any)
 * @param iv2  Rank of the second additional element designation (or -1 for any)
 * @param igrf Rank of Target GRF (or -1 for any)
 * @return -1 if no match is found
 */
int ANoStat::getRank(const EConsElem& type, int icov, int iv1, int iv2, int igrf) const
{
  if (_items.empty()) return -1;
  for (int ipar = 0; ipar < (int) getNoStatElemNumber(); ipar++)
  {
    if (matchIGrf(ipar,igrf) &&
        matchICov(ipar,icov) &&
        matchType(ipar,type) &&
        matchIV1 (ipar,iv1)  &&
        matchIV2 (ipar,iv2)) return ipar;
  }
  return -1;
}

/**
 * Add an element to the pile of Non-stationary parameters
 * @param igrf Rank of the Target GRF
 * @param icov Rank of the target covariance
 * @param type Type of the non-stationary parameter
 * @param iv1  Rank of the first designation
 * @param iv2  Rank of the second designation
 */
int ANoStat::addNoStatElem(int igrf, int icov, const EConsElem& type, int iv1, int iv2)
{
  int nelem = static_cast<int> (_items.size());
  _items.resize(nelem+1);
  _items[nelem].init(igrf, icov, type, iv1, iv2);
  if (! _checkConsistency())
  {
    deleteNoStatElem(nelem);
    return 1;
  }
  return 0;
}

int ANoStat::addNoStatElemByItem(const CovParamId& item)
{
  return addNoStatElem(item.getIGrf(), item.getICov(), item.getType(),
                       item.getIV1(), item.getIV2());
}

int ANoStat::addNoStatElems(const VectorString &codes)
{
  int igrf, icov, iv1, iv2;
  EConsElem type;

  int nerr = 0;
  for (int i = 0; i < (int) codes.size(); i++)
  {
    if (_understandCode(codes[i], &igrf, &icov, &type, &iv1, &iv2)) continue;
    nerr += addNoStatElem(igrf, icov, type, iv1, iv2);
  }
  return (nerr > 0);
}

void ANoStat::deleteNoStatElem(int ipar)
{
  _items.erase(_items.begin() + ipar);
}

void ANoStat::deleteAllNoStatElem()
{
  _items.clear();
}

/**
 * This function tries to attach the current Non-Stationary environment
 * to the model, checking that references to Covariance ranks are valid
 * @param model Model structure
 * @return Error return code
 */
int ANoStat::attachModel(const Model* model)
{
  if (model == nullptr)
  {
    messerr("Model must be defined beforehand");
    return 1;
  }
  if (model->getDimensionNumber() > 3)
  {
    messerr("Non stationary model is restricted to Space Dimension <= 3");
    return 1;
  }

  // Patch. It may happen that the Model already contains the parameters of the ANostat
  // which are better defined than in the Current ANostat structure
  // In this case copy the ANostat parameters from Model to Current
  // TODO: Remove this part of code
  _updateFromModel(model);

  for (int ipar=0; ipar<(int) getNoStatElemNumber(); ipar++)
  {
    int icov = getICov(ipar);
    EConsElem type = getType(ipar);

    // Check that the Non-stationary parameter is valid with respect
    // to the Model definition

    if (icov < 0 || icov >= model->getCovaNumber())
    {
      messerr("Invalid Covariance rank (%d) for the Non-Stationary Parameter (%d)",
              icov,ipar);
      return 1;
    }
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
 * @param igrf Rank of the GRF
 * @param icov Rank of the Covariance
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
int ANoStat::_understandCode(const String& code,
                             int *igrf,
                             int *icov,
                             EConsElem *type,
                             int *iv1,
                             int *iv2)
{
  *igrf = *icov = *iv1 = *iv2 = 0;
  *type = EConsElem::UNKNOWN;
  VectorString keywords = separateKeywords(code);

  bool flagGRF = false;
  bool flagSTR = false;
  bool flagV1  = false;
  bool flagV2  = false;
  int size = static_cast<int> (keywords.size());
  int lec = 0;

  // Decoding the GRF keyword (keyword "G")
  if (lec < size)
  {
    if (matchKeyword(keywords[lec],"G",false))
    {
      flagGRF = true;
      lec++;
    }
  }

  // Decoding the rank of the GRF (conditional)
  if (lec < size && flagGRF)
  {
    *igrf = toInteger(keywords[lec]);
    if (IFFFF(*igrf)) return 1;
    (*igrf)--;
    lec++;
  }

  // Decoding the structure keyword (keyword "M")
  if (lec < size)
  {
    if (matchKeyword(keywords[lec],"M",false))
    {
      flagSTR = true;
      lec++;
    }
  }

  // Decoding the rank of the Structure (conditional)
  if (lec < size && flagSTR)
  {
    *icov = toInteger(keywords[lec]);
    if (IFFFF(*icov)) return 1;
    (*icov)--;
    lec++;
  }

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
  if (lec < size && flagV2)
  {
    *iv2 = toInteger(keywords[lec]);
    if (IFFFF(*iv2)) return 1;
    (*iv2)--;
    lec++;
  }
  return 0;
}

/**
 * The Model may contain some non-stationary parameters already initialized
 * that the current structure is not aware of. They are simply duplicated here
 * This temporary patch is applied in the following conditions:
 * 1) the number of stationary elements in current structure is 0
 * 2) the number of stationary elements in Model is positive
 * @param model Input Model
 */
void ANoStat::_updateFromModel(const Model* model)
{
  int nelemFromModel = model->getNoStatElemNumber();
  int nelem = getNoStatElemNumber();
  if (nelem > 0) return;
  if (nelemFromModel <= 0) return;

  for (int ipar = 0; ipar < nelemFromModel; ipar++)
  {
    CovParamId item = model->getCovParamId(ipar);
    (void) addNoStatElemByItem(item);
  }
}

/**
 * Update the Model according to the Non-stationary parameters
 * @param model Model to be patched
 * @param icas1 Type of first Db: 1 for Input; 2 for Output
 * @param iech1 Rank of the target within Db1 (or -1)
 * @param icas2 Type of first Db: 1 for Input; 2 for Output
 * @param iech2 Rank of the target within Dbout (or -2)
 */
void ANoStat::updateModel(Model* model,
                          int icas1,
                          int iech1,
                          int icas2,
                          int iech2) const
{
  double val1, val2;

  // If no non-stationary parameter is defined, simply skip
  if (! model->isNoStat()) return;

  // Loop on the elements that can be updated one-by-one

  for (int ipar = 0; ipar < getNoStatElemNumber(); ipar++)
  {
    int icov = getICov(ipar);
    EConsElem type = getType(ipar);

    if (type == EConsElem::SILL)
    {
      _getInfoFromDb(ipar, icas1, iech1, icas2, iech2, &val1, &val2);
      int iv1  = getIV1(ipar);
      int iv2  = getIV2(ipar);
      model->setSill(icov, iv1, iv2, sqrt(val1 * val2));
    }
    else if (type == EConsElem::PARAM)
    {
      _getInfoFromDb(ipar, icas1, iech1, icas2, iech2, &val1, &val2);
      model->getCova(icov)->setParam(0.5 * (val1 + val2));
    }
  }

  // Loop on the other parameters (Anisotropy) that must be processed globally

  for (int icov = 0; icov < model->getCovaNumber(); icov++)
  {
    if (! isDefinedforAnisotropy(icov)) continue;
    CovAniso* cova = model->getCova(icov);

    VectorDouble angle0(cova->getAnisoAngles());
    VectorDouble angle1(angle0);
    VectorDouble angle2(angle0);

    VectorDouble scale0(cova->getScales());
    VectorDouble scale1(scale0);
    VectorDouble scale2(scale0);

    VectorDouble range0(cova->getRanges());
    VectorDouble range1(range0);
    VectorDouble range2(range0);

    // Define the angles (for all space dimensions)
    bool flagRot = false;
    if (isDefined(EConsElem::ANGLE, icov))
    {
      flagRot = true;
      for (int idim = 0; idim < model->getDimensionNumber(); idim++)
      {
        if (isDefined(EConsElem::ANGLE, icov, idim, 0))
        {
          int ipar = getRank(EConsElem::ANGLE, icov, idim);
          if (ipar < 0) continue;
          _getInfoFromDb(ipar, icas1, iech1, icas2, iech2, &angle1[idim], &angle2[idim]);
        }
      }
    }

    // Define the Theoretical ranges (for all space dimensions)

    bool flagScale = false;
    if (isDefined(EConsElem::SCALE, icov))
    {
      flagScale = true;
      for (int idim = 0; idim < model->getDimensionNumber(); idim++)
      {
        if (isDefined(EConsElem::SCALE, icov, idim, 0))
        {
          int ipar = getRank(EConsElem::SCALE, icov, idim);
          if (ipar < 0) continue;
          _getInfoFromDb(ipar, icas1, iech1, icas2, iech2, &scale1[idim], &scale2[idim]);
        }
      }
    }

    // Define the Practical ranges (for all space dimensions)

    bool flagRange = false;
    if (isDefined(EConsElem::RANGE, icov))
    {
      flagRange = true;
      for (int idim = 0; idim < model->getDimensionNumber(); idim++)
      {
        if (isDefined(EConsElem::RANGE, icov, idim))
        {
          int ipar = getRank(EConsElem::RANGE, icov, idim);
          if (ipar < 0) continue;
          _getInfoFromDb(ipar, icas1, iech1, icas2, iech2, &range1[idim], &range2[idim]);
        }
      }
    }

    // Create and exploit the Tensor for First location
    if (flagRot || flagRange || flagScale)
    {
      if (flagRot)   cova->setAnisoAngles(angle1);
      if (flagRange) cova->setRanges(range1);
      if (flagScale) cova->setScales(scale1);
      const MatrixSquareGeneral direct1 = cova->getAniso().getTensorDirect();

      // Create and exploit the Tensor for the second location
      if (flagRot)   cova->setAnisoAngles(angle2);
      if (flagRange) cova->setRanges(range2);
      if (flagScale) cova->setScales(scale2);

      // Build the new Tensor (as average of tensors at end-points)
      Tensor tensor = cova->getAniso();
      MatrixSquareGeneral direct = tensor.getTensorDirect();
      direct.linearCombination(0.5, 0.5, direct1);
      tensor.setTensorDirect(direct);
      cova->setAniso(tensor);
    }
  }
}

/**
 * Update the Model according to the Non-stationary parameters
 * @param model Model to be patched
 * @param imesh Rank of the target mesh
 */
void ANoStat::updateModelByMesh(Model* model, int imesh) const
{
  // If no non-stationary parameter is defined, simply skip
  if (! model->isNoStat()) return;

  // Loop on the elements that can be updated one-by-one

  for (int ipar = 0; ipar < model->getNoStatElemNumber(); ipar++)
  {
    int icov = getICov(ipar);
    EConsElem type = getType(ipar);

    if (type == EConsElem::SILL)
    {
      double sill = getValueByParam(ipar, 0, imesh);
      int iv1  = getIV1(ipar);
      int iv2  = getIV2(ipar);
      model->setSill(icov, iv1, iv2, sill);
    }
  }

  // Loop on the other parameters (Anisotropy) that must be processed globally

  for (int icov = 0; icov < model->getCovaNumber(); icov++)
  {
    if (! isDefinedforAnisotropy(icov)) continue;
    CovAniso* cova = model->getCova(icov);

    VectorDouble angle(cova->getAnisoAngles());
    VectorDouble scale(cova->getScales());
    VectorDouble range(cova->getRanges());

    // Define the angles (for all space dimensions)
    bool flagRot = false;
    if (isDefined(EConsElem::ANGLE, icov))
    {
      flagRot = true;
      for (int idim = 0; idim < model->getDimensionNumber(); idim++)
      {
        if (isDefined(EConsElem::ANGLE, icov, idim, 0))
        {
          int ipar = getRank(EConsElem::ANGLE, icov, idim);
          if (ipar < 0) continue;
          angle[idim] = getValueByParam(ipar, 0, imesh);
        }
      }
    }

    // Define the Theoretical ranges (for all space dimensions)

    bool flagScale = false;
    if (isDefined(EConsElem::SCALE, icov))
    {
      flagScale = true;
      for (int idim = 0; idim < model->getDimensionNumber(); idim++)
      {
        if (isDefined(EConsElem::SCALE, icov, idim))
        {
          int ipar = getRank(EConsElem::SCALE, icov, idim);
          if (ipar < 0) continue;
          scale[idim] = getValueByParam(ipar, 0, imesh);
        }
      }
    }

    // Define the Practical ranges (for all space dimensions)

    bool flagRange = false;
    if (isDefined(EConsElem::RANGE, icov))
    {
      flagRange = true;
      for (int idim = 0; idim < model->getDimensionNumber(); idim++)
      {
        if (isDefined(EConsElem::RANGE, icov, idim))
        {
          int ipar = getRank(EConsElem::RANGE, icov, idim);
          if (ipar < 0) continue;
          range[idim] = getValueByParam(ipar, 0, imesh);
        }
      }
    }

    // Exploit the Anisotropy
    if (flagRot)   cova->setAnisoAngles(angle);
    if (flagRange) cova->setRanges(range);
    if (flagScale) cova->setScales(scale);
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
 */
void ANoStat::_getInfoFromDb(int ipar,
                             int icas1,
                             int iech1,
                             int icas2,
                             int iech2,
                             double *val1,
                             double *val2) const
{
  *val1 = getValueByParam(ipar, icas1, iech1);
  *val2 = getValueByParam(ipar, icas2, iech2);

  if (FFFF(*val1) && FFFF(*val2)) return;

  if (! FFFF(*val1)) *val2 = *val1;
  if (! FFFF(*val2)) *val1 = *val2;
}

int ANoStat::attachToMesh(const AMesh* mesh, bool /*verbose*/) const
{
  _setAmesh(mesh);
  return 0;
}

int ANoStat::attachToDb(Db* db, int icas, bool /*verbose*/) const
{
  if (icas == 1)
    _setDbin(db);
  else
    _setDbout(db);
  return 0;
}

void ANoStat::detachFromMesh() const
{
  _setAmesh(nullptr);
}

void ANoStat::detachFromDb(Db* /*db*/, int icas) const
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
bool ANoStat::_checkConsistency() const
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
    messerr("Error in the definition of Model Non-stationarity");
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
bool ANoStat::_isValid(int icas, int rank) const
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
      if (_dbout == nullptr)
      {
        messerr("Checking the validity of the argument");
        messerr("Dbout: This requires '_dbout' to be defined beforehand");
        return false;
      }
      if (rank < 0 || rank >= _dbout->getSampleNumber(0))
      {
        messerr("Check the validity of the argument");
        messerr("Dbout: 'rank' (%d) should be smaller than number of samples (%d)",
                rank, _dbout->getSampleNumber(0));
        return false;
      }
      break;
  }
  return true;
}
