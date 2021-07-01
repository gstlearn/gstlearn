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
#include "Model/ANoStat.hpp"

#include "Model/Model.hpp"
#include "Basic/String.hpp"
#include "Basic/Utilities.hpp"
#include "geoslib_f.h"

ANoStat::ANoStat()
  : _items()
{
}

ANoStat::ANoStat(const VectorString& codes)
  : _items()
{
  addNoStatElems(codes);
}

ANoStat::ANoStat(const ANoStat &m)
  : _items(m._items)
{
}

ANoStat& ANoStat::operator= (const ANoStat &m)
{
  _items       = m._items;
  return *this;
}

ANoStat::~ANoStat()
{
}

String ANoStat::toString(int level) const
{
  std::stringstream sstr;
  if (getNoStatElemNumber() <= 0) return sstr.str();

  sstr << toTitle(1, "Non-Stationary Parameters");
  for (int i = 0; i < (int) getNoStatElemNumber(); i++)
    sstr << _items[i].toString(level);
  return sstr.str();
}

/**
 * Look if a Non-stationary parameter is defined
 * @param igrf Rank of Target GRF (or -1 for any)
 * @param icov Rank of Target Covariance (or -1 for any)
 * @return
 */
bool ANoStat::isDefinedByCov(int igrf, int icov)
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
 * @param igrf Rank of Target GRF (or -1 for any)
 * @param type Rank of Target Type (or CONS_UNKNOWN for any)
 * @return
 */
bool ANoStat::isDefinedByType(int igrf, ENUM_CONS type)
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
 * @param igrf Rank of Target GRF (or -1 for any)
 * @param icov Rank of the Target Covariance (or -1 for any)
 * @param type Rank of Target Type (or CONS_UNKNOWN for any)
 * @return
 */
bool ANoStat::isDefinedByCovType(int igrf, int icov, ENUM_CONS type)
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

bool ANoStat::isDefined(int igrf, int icov, ENUM_CONS type, int iv1, int iv2)
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
 * @param igrf Rank of Target GRF (or -1 for any)
 * @param icov Rank of the Target Covariance (or -1 for any)
 * @return
 */
bool ANoStat::isDefinedforAnisotropy(int igrf, int icov)
{
  if (_items.empty()) return false;
  for (int ipar = 0; ipar < (int) getNoStatElemNumber(); ipar++)
  {
    if (! matchIGrf(ipar,igrf)) continue;
    if (! matchICov(ipar,icov)) continue;
    if (getType(ipar) == CONS_ANGLE ||
        getType(ipar) == CONS_RANGE ||
        getType(ipar) == CONS_SCALE) return true;
  }
  return false;
}

/**
 * Return the rank for a Non-stationary parameter
 * @param igrf Rank of Target GRF (or -1 for any)
 * @param icov Rank of the Target Covariance (or -1 for any)
 * @param type Rank of Target Type (or CONS_UNKNOWN for any)
 * @param iv1  Rank of the first additional element designation (or -1 for any)
 * @param iv2  Rank of the second additional element designation (or -1 for any)
 * @return -1 if no match is found
 */
int ANoStat::getRank(int igrf, int icov, ENUM_CONS type, int iv1, int iv2) const
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
void ANoStat::addNoStatElem(int igrf, int icov, ENUM_CONS type, int iv1, int iv2)
{
  int nelem = _items.size();
  _items.resize(nelem+1);
  _items[nelem].init(CONS_TYPE_DEFAULT, igrf, icov, type, iv1, iv2, TEST);
}

void ANoStat::addNoStatElem(const ConsItem& item)
{
  addNoStatElem(item.getIGrf(),item.getICov(),item.getType(),item.getIV1(),item.getIV2());
}

void ANoStat::addNoStatElems(const VectorString &codes)
{
  int igrf, icov, iv1, iv2;
  ENUM_CONS type;

  for (int i = 0; i < (int) codes.size(); i++)
  {
    if (_understandCode(codes[i], &igrf, &icov, &type, &iv1, &iv2)) continue;
    addNoStatElem(igrf, icov, type, iv1, iv2);
  }
}

/**
 * This function tries to attach the current Non-Stationary environment
 * to the model, checking that references to Covariance ranks are valid
 * @param model Model structure
 * @return Error return code
 */
const int ANoStat::attachModel(const Model* model)
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

  // Patch. It may happen that the Model already contain the parameters of the ANostat
  // which are better defined than in the Current ANostat structure
  // In this case copy the ANostat parameters from Model to Current
  // TODO: Remove this part of code
  _updateFromModel(model);

  for (int ipar=0; ipar<(int) getNoStatElemNumber(); ipar++)
  {
    int icov = getICov(ipar);
    ENUM_CONS type = getType(ipar);

    // Check that the Non-stationary parameter is valid with respect
    // to the Model definition

    if (icov < 0 || icov >= model->getCovaNumber())
    {
      messerr("Invalid Covariance rank (%d) for the Non-Stationary Parameter (%d)",
              icov,ipar);
      return 1;
    }
    if (type == CONS_PARAM)
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
 * @param type Type of parameter ("R","A","P","V","S","T","C","I")
 * @param iv1  Rank of the first variable
 * @param iv2  Rank of the second variable
 * @return 0 if a valid constraint has been defined; 1 otherwise
 */
int ANoStat::_understandCode(const String& code,
                             int *igrf,
                             int *icov,
                             ENUM_CONS *type,
                             int *iv1,
                             int *iv2)
{
  *igrf = *icov = *iv1 = *iv2 = 0;
  *type = CONS_UNKNOWN;
  VectorString keywords = separateKeywords(code);

  bool flagGRF = false;
  bool flagSTR = false;
  bool flagV1  = false;
  bool flagV2  = false;
  int size = keywords.size();
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
    *igrf = toInt(keywords[lec]);
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
    *icov = toInt(keywords[lec]);
    if (IFFFF(*icov)) return 1;
    (*icov)--;
    lec++;
  }

  // Decoding the operator keyword ("R","A","P","V","S","T","C","I")
  if (lec < size)
  {
    // TODO better decoding which would not be dependent on the order in ENUM_CONS
    VectorString list = {"R","A","P","V","S","T","C","I"};
    (*type) = (ENUM_CONS) (getRankInList(list, keywords[lec], false) + 1);
    if ((*type) < 0) return 1;
    flagV1 = true;
    lec++;
  }

  // Decoding the rank of the operator (conditional)
  if (lec < size && flagV1)
  {
    *iv1 = toInt(keywords[lec]);
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
    *iv2 = toInt(keywords[lec]);
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
    ConsItem item = model->getConsItem(ipar);
    addNoStatElem(item);
  }
}
