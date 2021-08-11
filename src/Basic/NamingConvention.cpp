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
#include "Basic/NamingConvention.hpp"
#include "Basic/Vector.hpp"
#include "Basic/String.hpp"
#include "Db/Db.hpp"
#include <string>

NamingConvention::NamingConvention(String radix,
                                   ENUM_LOCS locatorOutType,
                                   String delim,
                                   bool flagvariter,
                                   bool flagclean)
    : _radix(radix),
      _delim(delim),
      _locatorOutType(locatorOutType),
      _flagVarIter(flagvariter),
      _flagClean(flagclean)
{

}

NamingConvention::NamingConvention(const NamingConvention &m)
    : _radix(m._radix),
      _delim(m._delim),
      _locatorOutType(m._locatorOutType),
      _flagVarIter(m._flagVarIter),
      _flagClean(m._flagClean)
{

}

NamingConvention& NamingConvention::operator=(const NamingConvention &m)
{
  if (this != &m)
  {
    _radix = m._radix;
    _locatorOutType = m._locatorOutType;
    _delim = m._delim;
    _flagVarIter = m._flagVarIter;
    _flagClean = m._flagClean;
  }
  return *this;
}

NamingConvention::~NamingConvention()
{

}

/**
 * Naming from a set of variables identified by their names
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param suffix Optional suffix
 * @param nitems Number of items
 * @param flagLocate True if the variable must be assigned the locator
 */
void NamingConvention::setNamesAndLocators(Db* dbout,
                                           int iattout_start,
                                           const String& suffix,
                                           int nitems,
                                           bool flagLocate) const
{
  _setNames(dbout, iattout_start, VectorString(), suffix, nitems);
  setLocators(dbout, iattout_start, 1, nitems, flagLocate);
}

/**
 * Naming from a set of variables identified by their names
 * @param names Vector of variable names in Input Db
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param suffix Optional suffix
 * @param nitems Number of items
 * @param flagLocate True if the variable must be assigned the locator
 */
void NamingConvention::setNamesAndLocators(const VectorString& names,
                                           Db* dbout,
                                           int iattout_start,
                                           const String& suffix,
                                           int nitems,
                                           bool flagLocate) const
{
  if (iattout_start <= 0) return;
  int nvar = static_cast<int> (names.size());
  if (nvar <= 0) return;
  _setNames(dbout, iattout_start, names, suffix, nitems);
  setLocators(dbout, iattout_start, nvar, nitems, flagLocate);
}

/**
 * Naming from one variable identified by its names
 * @param namin variable name in Input Db
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param suffix Optional suffix
 * @param nitems Number of items
 * @param flagLocate True if the variable must be assigned the locator
 */
void NamingConvention::setNamesAndLocators(const String& namin,
                                           Db* dbout,
                                           int iattout_start,
                                           const String& suffix,
                                           int nitems,
                                           bool flagLocate) const
{
  if (iattout_start <= 0) return;
  VectorString names;
  names.push_back(namin);
  _setNames(dbout, iattout_start, names, suffix, nitems);
  setLocators(dbout, iattout_start, 1, nitems, flagLocate);
}

/**
 * Naming from a set of variables identified by their locatorType
 * @param dbin  Pointer to the input Db (kept for symmetry)
 * @param locatorInType Locator Tyoe of the variables in Input Db
 * @param nvar Number of items belonging to the locatorType
 *             (if -1, all the items available for this locator are used)
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param suffix Optional suffix
 * @param nitems Number of items
 * @param flagLocate True if the variable must be assigned the locator
 */
void NamingConvention::setNamesAndLocators(Db *dbin,
                                           ENUM_LOCS locatorInType,
                                           int nvar,
                                           Db* dbout,
                                           int iattout_start,
                                           const String& suffix,
                                           int nitems,
                                           bool flagLocate) const
{
  if (iattout_start <= 0) return;
  VectorString names;
  if (dbin != nullptr && locatorInType == LOC_UNKNOWN)
  {
    names = dbin->getNames(locatorInType);
    if (nvar <= 0) nvar = static_cast<int> (names.size());
    if (nvar < (int) names.size()) names.resize(nvar);
  }
  else
  {
    if (nvar < 0) nvar = 1;
  }
  _setNames(dbout, iattout_start, names, suffix, nitems);
  setLocators(dbout, iattout_start, nvar, nitems, flagLocate);
}

/**
 * Naming from a set of variables identified by their attribute indices
 * @param dbin  Pointer to the input Db (kept for symmetry)
 * @param iatts Vector of attribute indices of the variables in Input Db
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param suffix Optional suffix
 * @param nitems Number of items
 * @param flagLocate True if the variable must be assigned the locator
 */
void NamingConvention::setNamesAndLocators(Db *dbin,
                                           const VectorInt& iatts,
                                           Db* dbout,
                                           int iattout_start,
                                           const String& suffix,
                                           int nitems,
                                           bool flagLocate) const
{
  if (iattout_start <= 0) return;
  if (dbin == nullptr) return;
  int nvar = static_cast<int> (iatts.size());
  if (nvar <= 0) return;

  VectorString names;
  for (int ivar = 0; ivar < nvar; ivar++)
    names.push_back(dbin->getName(iatts[ivar]));
  _setNames(dbout, iattout_start, names, suffix, nitems);
  setLocators(dbout, iattout_start, nvar, nitems, flagLocate);
 }

/**
 * Naming from a set of variables identified by their attribute indices
 * @param dbin  Pointer to the input Db (kept for symmetry)
 * @param iatt  Attribute index of the variables in Input Db
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param suffix Optional suffix
 * @param nitems Number of items
 * @param flagLocate True if the variable must be assigned the locator
 */
void NamingConvention::setNamesAndLocators(Db *dbin,
                                           int iatt,
                                           Db* dbout,
                                           int iattout_start,
                                           const String& suffix,
                                           int nitems,
                                           bool flagLocate) const
{
  if (iattout_start <= 0) return;
  if (dbin == nullptr) return;

  VectorString names;
  names.push_back(dbin->getName(iatt));
  _setNames(dbout, iattout_start, names, suffix, nitems);
  setLocators(dbout, iattout_start, 1, nitems, flagLocate);
 }

void NamingConvention::setLocators(Db *dbout,
                                   int iattout_start,
                                   int nvar,
                                   int nitems,
                                   bool flagLocate) const
{
  if (! flagLocate) return;
  if (_locatorOutType == LOC_UNKNOWN) return;

  // Erase already existing locators of the same Type
  if (_flagClean)
    dbout->clearLocators(_locatorOutType);

  // Set the locator for all variables
  for (int ecr = 0; ecr < nvar * nitems; ecr++)
    dbout->setLocatorByAttribute(iattout_start + ecr, _locatorOutType, ecr + 1);
}

/**
 * Defines the names given to the output variables
 *
 * @param dbout   Pointer to the output Db structure
 * @param iattout_start Rank of the first variable to be named
 * @param names Vector of Names (dimension: nvar) or String()
 * @param suffix Suffix provided to construct the names
 * @param nitems Number of items to be renamed
 *
 */
void NamingConvention::_setNames(Db *dbout,
                                 int iattout_start,
                                 VectorString names,
                                 const String& suffix,
                                 int nitems) const
{
  int ecr = 0;
  int nvar = (names.empty()) ? 1 : static_cast<int> (names.size());

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    String local = (names.empty() || (int) names.size() != nvar) ? String() : names[ivar];

    for (int item = 0; item < nitems; item++)
    {
      String locnum = (nitems <= 1) ? String() : std::to_string(item+1);
      String name = concatenateStrings(_delim,_radix, local, suffix, locnum);
      dbout->setName(iattout_start + ecr, name);
      ecr++;
    }
  }
}
