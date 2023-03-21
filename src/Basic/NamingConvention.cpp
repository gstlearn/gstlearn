/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
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
 * Naming a set of variables of 'dbout' identified by their names
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
 * Naming from a set of variables identified by their names
 * @param names Vector of variable names in Input Db
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
 * Naming given the set of output variable names
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
 * Naming from one variable identified by its names
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
 * Naming from a set of variables identified by their attribute indices
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
 * Naming from a set of variables identified by their attribute indices
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
  VectorString outnames = createNames(names, qualifier, nitems);

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
VectorString NamingConvention::createNames(const VectorString& names,
                                           const String& qualifier,
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

      // Compose the variable name
      String name = concatenateStrings(_delim, _prefix,
                                       loc_varname, loc_qualifier, loc_number);

      outnames.push_back(name);
    }
  }
  return outnames;
}
