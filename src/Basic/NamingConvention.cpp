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

NamingConvention::NamingConvention(String prefix,
                                   bool flag_varname,
                                   bool flag_qualifier,
                                   const ELoc &locatorOutType,
                                   String delim,
                                   bool cleanSameLocator)
    : _prefix(prefix),
      _delim(delim),
      _flagVarname(flag_varname),
      _flagQualifier(flag_qualifier),
      _locatorOutType(locatorOutType),
      _cleanSameLocator(cleanSameLocator)
{
}

NamingConvention::NamingConvention(const NamingConvention &m)
    : _prefix(m._prefix),
      _delim(m._delim),
      _flagVarname(m._flagVarname),
      _flagQualifier(m._flagQualifier),
      _locatorOutType(m._locatorOutType),
      _cleanSameLocator(m._cleanSameLocator)
{
}

NamingConvention& NamingConvention::operator=(const NamingConvention &m)
{
  if (this != &m)
  {
    _prefix = m._prefix;
    _locatorOutType = m._locatorOutType;
    _flagVarname = m._flagVarname;
    _flagQualifier = m._flagQualifier;
    _delim = m._delim;
    _cleanSameLocator = m._cleanSameLocator;
  }
  return *this;
}

NamingConvention::~NamingConvention()
{
}

/**
 * Naming a set of variables of 'dbout' identified by their names
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param qualifier Optional qualifier
 * @param nitems Number of items
 * @param flagSetLocator True if the variable must be assigned the locator
 */
void NamingConvention::setNamesAndLocators(Db* dbout,
                                           int iattout_start,
                                           const String& qualifier,
                                           int nitems,
                                           bool flagSetLocator) const
{
  _setNames(dbout, iattout_start, VectorString(), qualifier, nitems);

  if (flagSetLocator)
    setLocators(dbout, iattout_start, 1, nitems);
}

/**
 * Naming from a set of variables identified by their names
 * @param names Vector of variable names in Input Db
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param qualifier Optional qualifier
 * @param nitems Number of items
 * @param flagSetLocator True if the variable must be assigned the locator
 */
void NamingConvention::setNamesAndLocators(const VectorString& names,
                                           Db* dbout,
                                           int iattout_start,
                                           const String& qualifier,
                                           int nitems,
                                           bool flagSetLocator) const
{
  if (iattout_start < 0) return;
  int nvar = static_cast<int> (names.size());
  if (nvar <= 0) return;

  _setNames(dbout, iattout_start, names, qualifier, nitems);

  if (flagSetLocator)
    setLocators(dbout, iattout_start, nvar, nitems);
}

/**
 * Naming given the set of output variable names
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param names Vector of output variable names
 * @param flagSetLocator True if the variable must be assigned the locator
 */
void NamingConvention::setNamesAndLocators(Db* dbout,
                                           int iattout_start,
                                           const VectorString& names,
                                           bool flagSetLocator) const
{
  if (iattout_start < 0) return;
  int nvar = static_cast<int> (names.size());

  for (int ivar = 0; ivar < nvar; ivar++)
    dbout->setNameByUID(iattout_start+ivar, names[ivar]);

  if (flagSetLocator)
    setLocators(dbout, iattout_start, nvar, 1);
}

/**
 * Naming from one variable identified by its names
 * @param namin variable name in Input Db
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param qualifier Optional qualifier
 * @param nitems Number of items
 * @param flagSetLocator True if the variable must be assigned the locator
 */
void NamingConvention::setNamesAndLocators(const String& namin,
                                           Db* dbout,
                                           int iattout_start,
                                           const String& qualifier,
                                           int nitems,
                                           bool flagSetLocator) const
{
  if (iattout_start < 0) return;
  VectorString names;
  names.push_back(namin);

  _setNames(dbout, iattout_start, names, qualifier, nitems);

  if (flagSetLocator)
    setLocators(dbout, iattout_start, 1, nitems);
}

/**
 * Naming from a set of variables identified by their locatorType
 * @param dbin  Pointer to the input Db (kept for symmetry)
 * @param locatorInType Locator Tyoe of the variables in Input Db
 * @param nvar Number of items belonging to the locatorType
 *             (if -1, all the items available for this locator are used)
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param qualifier Optional qualifier
 * @param nitems Number of items
 * @param flagSetLocator True if the variable must be assigned the locator
 */
void NamingConvention::setNamesAndLocators(const Db *dbin,
                                           const ELoc& locatorInType,
                                           int nvar,
                                           Db* dbout,
                                           int iattout_start,
                                           const String& qualifier,
                                           int nitems,
                                           bool flagSetLocator) const
{
  if (iattout_start < 0) return;
  VectorString names;
  if (dbin != nullptr && locatorInType != ELoc::UNKNOWN)
  {
    names = dbin->getNamesByLocator(locatorInType);
    if (nvar <= 0) nvar = static_cast<int>(names.size());
  }
  else
  {
    if (nvar < 0) nvar = 1;
  }
  if (nvar != static_cast<int>(names.size())) names.resize(nvar);

  _setNames(dbout, iattout_start, names, qualifier, nitems);

  if (flagSetLocator)
    setLocators(dbout, iattout_start, nvar, nitems);
}

/**
 * Naming from a set of variables identified by their attribute indices
 * @param dbin  Pointer to the input Db (kept for symmetry)
 * @param iatts Vector of attribute indices of the variables in Input Db
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param qualifier Optional qualifier
 * @param nitems Number of items
 * @param flagSetLocator True if the variable must be assigned the locator
 */
void NamingConvention::setNamesAndLocators(const Db *dbin,
                                           const VectorInt& iatts,
                                           Db* dbout,
                                           int iattout_start,
                                           const String& qualifier,
                                           int nitems,
                                           bool flagSetLocator) const
{
  if (iattout_start < 0) return;
  if (dbin == nullptr) return;
  int nvar = static_cast<int>(iatts.size());
  if (nvar <= 0) return;

  VectorString names;
  for (int ivar = 0; ivar < nvar; ivar++)
    names.push_back(dbin->getNameByUID(iatts[ivar]));

  _setNames(dbout, iattout_start, names, qualifier, nitems);

  if (flagSetLocator)
    setLocators(dbout, iattout_start, nvar, nitems);
 }

/**
 * Naming from a set of variables identified by their attribute indices
 * @param dbin  Pointer to the input Db (kept for symmetry)
 * @param iatt  Attribute index of the variables in Input Db
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param qualifier Optional qualifier
 * @param nitems Number of items
 * @param flagSetLocator True if the variable must be assigned the locator
 */
void NamingConvention::setNamesAndLocators(const Db *dbin,
                                           int iatt,
                                           Db* dbout,
                                           int iattout_start,
                                           const String& qualifier,
                                           int nitems,
                                           bool flagSetLocator) const
{
  if (iattout_start < 0) return;
  if (dbin == nullptr) return;

  VectorString names;
  names.push_back(dbin->getNameByUID(iatt));

  _setNames(dbout, iattout_start, names, qualifier, nitems);

  if (flagSetLocator)
    setLocators(dbout, iattout_start, 1, nitems);
 }

void NamingConvention::setLocators(Db *dbout,
                                   int iattout_start,
                                   int nvar,
                                   int nitems) const
{
  if (_locatorOutType == ELoc::UNKNOWN) return;

  // Erase already existing locators of the same Type
  if (_cleanSameLocator)
    dbout->clearLocators(_locatorOutType);

  // Set the locator for all variables
  for (int ecr = 0; ecr < nvar * nitems; ecr++)
    dbout->setLocatorByUID(iattout_start + ecr, _locatorOutType, ecr);
}

/**
 * Defines the names of the output variables. These variables are located
 * in 'dbout'; they have consecutive UIDs, starting from 'iattout_start'
 *
 * @param dbout   Pointer to the output Db structure
 * @param iattout_start Rank of the first variable to be named
 * @param names Vector of Names (dimension: nvar) or String()
 * @param qualifier Optional qualifier
 * @param nitems Number of items to be renamed
 */
void NamingConvention::_setNames(Db *dbout,
                                 int iattout_start,
                                 const VectorString& names,
                                 const String& qualifier,
                                 int nitems) const
{
  int ecr = 0;
  int nvar = (names.empty()) ? 1 : static_cast<int>(names.size());

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    // Variable 'local' defined for each variable, is:
    // - extracted from the array 'names' (if defined)
    // - generated as the rank of the variable (if several)
    String loc_varname;
    if (_flagVarname)
    {
      if (static_cast<int>(names.size()) == nvar) loc_varname = names[ivar];
      if (loc_varname.empty() && nvar > 1) loc_varname = std::to_string(ivar+1);
    }

    for (int item = 0; item < nitems; item++)
    {
      String loc_qualifier;
      String loc_number;
      if (_flagQualifier)
      {
        loc_qualifier = qualifier;
        if (nitems > 1) loc_number = std::to_string(item+1);
      }
      String name = concatenateStrings(_delim, _prefix,
                                       loc_varname, loc_qualifier, loc_number);
      dbout->setNameByUID(iattout_start + ecr, name);
      ecr++;
    }
  }
}
