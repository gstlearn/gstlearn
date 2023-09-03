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
#include "Basic/NamingConvention.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/String.hpp"
#include "Db/Db.hpp"

#include <string>

NamingConvention::NamingConvention(String prefix,
                                   bool flag_varname,
                                   bool flag_qualifier,
                                   bool flag_locator,
                                   const ELoc& locatorOutType,
                                   String delim,
                                   bool cleanSameLocator)
    : _prefix(prefix),
      _delim(delim),
      _flagVarname(flag_varname),
      _flagQualifier(flag_qualifier),
      _flagLocator(flag_locator),
      _locatorOutType(locatorOutType),
      _cleanSameLocator(cleanSameLocator)
{
}

NamingConvention::NamingConvention(const NamingConvention &m)
    : _prefix(m._prefix),
      _delim(m._delim),
      _flagVarname(m._flagVarname),
      _flagQualifier(m._flagQualifier),
      _flagLocator(m._flagLocator),
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
    _flagLocator = m._flagLocator;
    _delim = m._delim;
    _cleanSameLocator = m._cleanSameLocator;
  }
  return *this;
}

NamingConvention::~NamingConvention()
{
}

/**
 * Construct an item of the Naming Convention Class
 * @param prefix Name given to the prefix
 * @param flag_varname When TRUE, the 'varname' is included in the output names
 * @param flag_qualifier When TRUE, the 'qualifier' is included in the output names
 * @param flag_locator When TRUE, the output variables receive a 'locator'
 * @param locatorOutType Type of locator assigned to the output variables
 * @param delim Symbol used as a delimitor separating the different parts of the output names
 * @param cleanSameLocator When TRUE and if 'flag_locator' is TRUE, all variables assigned to the same locator are cancelled beforehand
 * @return
 */
NamingConvention* NamingConvention::create(String prefix,
                                           bool flag_varname,
                                           bool flag_qualifier,
                                           bool flag_locator,
                                           const ELoc &locatorOutType,
                                           String delim,
                                           bool cleanSameLocator)
{
  return new NamingConvention(prefix, flag_varname, flag_qualifier,
    flag_locator, locatorOutType, delim, cleanSameLocator);
}

/**
 * Newly created variables are named as follows:
 *
 * 'prefix'.'qualifier'.'item_rank'
 *
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param qualifier Optional qualifier
 * @param nitems Number of items
 * @param flagSetLocator True if the variable must be assigned the locator
 * @param locatorShift Shift to be applied to the locator currently defined
 */
void NamingConvention::setNamesAndLocators(Db* dbout,
                                           int iattout_start,
                                           const String& qualifier,
                                           int nitems,
                                           bool flagSetLocator,
                                           int locatorShift) const
{
  _setNames(dbout, iattout_start, VectorString(), qualifier, nitems);

  if (flagSetLocator)
    setLocators(dbout, iattout_start, 1, nitems, locatorShift);
}

/**
 * Newly created variables are named as follows:
 *
 * 'prefix'.'names[i]'.qualifier'.'item_rank'
 *
 * @param names Vector of variable names
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param qualifier Optional qualifier
 * @param nitems Number of items
 * @param flagSetLocator True if the variable must be assigned the locator
 * @param locatorShift Shift to be applied to the locator currently defined
 */
void NamingConvention::setNamesAndLocators(const VectorString& names,
                                           Db* dbout,
                                           int iattout_start,
                                           const String& qualifier,
                                           int nitems,
                                           bool flagSetLocator,
                                           int locatorShift) const
{
  if (iattout_start < 0) return;
  int nvar = static_cast<int> (names.size());
  if (nvar <= 0) return;

  _setNames(dbout, iattout_start, names, qualifier, nitems);

  if (flagSetLocator)
    setLocators(dbout, iattout_start, nvar, nitems, locatorShift);
}

/**
 * Newly created variables are named as follows:
 *
 * "names[i]"
 *
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param names Vector of output variable names
 * @param flagSetLocator True if the variable must be assigned the locator
 * @param locatorShift Shift to be applied to the locator currently defined
 */
void NamingConvention::setNamesAndLocators(Db* dbout,
                                           int iattout_start,
                                           const VectorString& names,
                                           bool flagSetLocator,
                                           int locatorShift) const
{
  if (iattout_start < 0) return;
  int nvar = static_cast<int> (names.size());

  for (int ivar = 0; ivar < nvar; ivar++)
    dbout->setNameByUID(iattout_start+ivar, names[ivar]);

  if (flagSetLocator)
    setLocators(dbout, iattout_start, nvar, 1, locatorShift);
}

/**
 * Newly created variables are named as follow:
 *
 * 'prefix'.'namin'.'qualifier'.'item_rank'
 *
 * @param namin variable name in Input Db
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param qualifier Optional qualifier
 * @param nitems Number of items
 * @param flagSetLocator True if the variable must be assigned the locator
 * @param locatorShift Shift to be applied to the locator currently defined
 */
void NamingConvention::setNamesAndLocators(const String& namin,
                                           Db* dbout,
                                           int iattout_start,
                                           const String& qualifier,
                                           int nitems,
                                           bool flagSetLocator,
                                           int locatorShift) const
{
  if (iattout_start < 0) return;
  VectorString names;
  names.push_back(namin);

  _setNames(dbout, iattout_start, names, qualifier, nitems);

  if (flagSetLocator)
    setLocators(dbout, iattout_start, 1, nitems, locatorShift);
}

/**
 * Newly created variables are named as follow:
 *
 * 'prefix'.'v_Loc'.'qualifier'.'item_rank'
 * where 'v_Loc' stands for the name of the variable(s) with locator 'Loc' in 'dbin'
 *
 * @param dbin  Pointer to the input Db (kept for symmetry)
 * @param locatorInType Locator Type of the variables in Input Db
 * @param nvar Number of items belonging to the locatorType
 *             (if -1, all the items available for this locator are used)
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param qualifier Optional qualifier
 * @param nitems Number of items
 * @param flagSetLocator True if the variable must be assigned the locator
 * @param locatorShift Shift to be applied to the locator currently defined
 */
void NamingConvention::setNamesAndLocators(const Db *dbin,
                                           const ELoc& locatorInType,
                                           int nvar,
                                           Db* dbout,
                                           int iattout_start,
                                           const String& qualifier,
                                           int nitems,
                                           bool flagSetLocator,
                                           int locatorShift) const
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
    setLocators(dbout, iattout_start, nvar, nitems, locatorShift);
}

/**
 * Newly created variables are named as follow:
 *
 * 'prefix'.'v[i]'.'qualifier'.'item_rank'
 * where v[i] is the variable with rank 'i' within 'dbin'
 *
 * @param dbin  Pointer to the input Db (kept for symmetry)
 * @param iatts Vector of attribute indices of the variables in Input Db
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param qualifier Optional qualifier
 * @param nitems Number of items
 * @param flagSetLocator True if the variable must be assigned the locator
 * @param locatorShift Shift to be applied to the locator currently defined
 */
void NamingConvention::setNamesAndLocators(const Db *dbin,
                                           const VectorInt& iatts,
                                           Db* dbout,
                                           int iattout_start,
                                           const String& qualifier,
                                           int nitems,
                                           bool flagSetLocator,
                                           int locatorShift) const
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
    setLocators(dbout, iattout_start, nvar, nitems, locatorShift);
 }

/**
 * Newly created variables are named as follow:
 *
 * 'prefix'.'v[iatt]'.'qualifier'.'item_rank'
 *
 * @param dbin  Pointer to the input Db (kept for symmetry)
 * @param iatt  Attribute index of the variables in Input Db
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param qualifier Optional qualifier
 * @param nitems Number of items
 * @param flagSetLocator True if the variable must be assigned the locator
 * @param locatorShift Shift to be applied to the locator currently defined
 */
void NamingConvention::setNamesAndLocators(const Db *dbin,
                                           int iatt,
                                           Db* dbout,
                                           int iattout_start,
                                           const String& qualifier,
                                           int nitems,
                                           bool flagSetLocator,
                                           int locatorShift) const
{
  if (iattout_start < 0) return;
  if (dbin == nullptr) return;

  VectorString names;
  names.push_back(dbin->getNameByUID(iatt));

  _setNames(dbout, iattout_start, names, qualifier, nitems);

  if (flagSetLocator)
    setLocators(dbout, iattout_start, 1, nitems, locatorShift);
 }

void NamingConvention::setLocators(Db *dbout,
                                   int iattout_start,
                                   int nvar,
                                   int nitems,
                                   int locatorShift) const
{
  if (! _flagLocator || _locatorOutType == ELoc::UNKNOWN) return;

  // Erase already existing locators of the same Type
  // (this is only done if you are not precisely adding higher order version for given locator)
  if (_cleanSameLocator && locatorShift == 0)
    dbout->clearLocators(_locatorOutType);

  // Set the locator for all variables
  for (int ecr = 0; ecr < nvar * nitems; ecr++)
    dbout->setLocatorByUID(iattout_start + ecr, _locatorOutType, ecr + locatorShift);
}

/**
 * Defines the names of the output variables. These variables are located
 * in 'dbout'; they have consecutive UIDs, starting from 'iattout_start'
 *
 * @param dbout   Pointer to the output Db structure
 * @param iattout_start Rank of the first variable to be named
 * @param names Vector of Names (dimension: nvar)
 * @param qualifier Optional qualifier
 * @param nitems Number of items to be renamed
 */
void NamingConvention::_setNames(Db *dbout,
                                 int iattout_start,
                                 const VectorString& names,
                                 const String& qualifier,
                                 int nitems) const
{
  VectorString outnames = _createNames(names, qualifier, nitems);

  int ecr = 0;
  int nvar = (names.empty()) ? 1 : static_cast<int>(names.size());
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    for (int item = 0; item < nitems; item++)
    {
      dbout->setNameByUID(iattout_start + ecr, outnames[ecr]);
      ecr++;
    }
  }
}

/**
 * Defines the names of the output variables.
 *
 * @param names Vector of Names (dimension: nvar)
 * @param qualifier Optional qualifier
 * @param nitems Number of items to be renamed
 */
VectorString NamingConvention::_createNames(const VectorString &names,
                                            const String &qualifier,
                                            int nitems) const
{
  VectorString outnames;

  int nvar = (names.empty()) ? 1 : static_cast<int>(names.size());

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    // Variable 'local' defined for each variable, is:
    // - extracted from the array 'names' (if defined)
    // - generated as the rank of the variable (if several)
    String loc_varname;
    String loc_number;
    if (_flagVarname)
    {
      if (static_cast<int>(names.size()) == nvar) loc_varname = names[ivar];
      if (loc_varname.empty() && nvar > 1) loc_varname = std::to_string(ivar+1);
    }
    else
    {
      // Build the rank from the variable number (possibly overwritten by item number)
      if (nvar > 1) loc_number = std::to_string(ivar+1);
    }

    for (int item = 0; item < nitems; item++)
    {
      String loc_qualifier;

      if (_flagQualifier)
      {
        loc_qualifier = qualifier;
        if (nitems > 1) loc_number = std::to_string(item+1);
      }

      // Compose the variable name
      String name = concatenateStrings(_delim, _prefix,
                                       loc_varname, loc_qualifier, loc_number);

      outnames.push_back(name);
    }
  }
  return outnames;
}
