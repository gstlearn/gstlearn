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
                                   ENUM_LOCS locatorType,
                                   String delim,
                                   bool flagvariter,
                                   bool flagclean)
    : _radix(radix),
      _delim(delim),
      _locatorType(locatorType),
      _flagVarIter(flagvariter),
      _flagClean(flagclean)
{

}

NamingConvention::NamingConvention(const NamingConvention &m)
    : _radix(m._radix),
      _delim(m._delim),
      _locatorType(m._locatorType),
      _flagVarIter(m._flagVarIter),
      _flagClean(m._flagClean)
{

}

NamingConvention& NamingConvention::operator=(const NamingConvention &m)
{
  if (this != &m)
  {
    _radix = m._radix;
    _locatorType = m._locatorType;
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
 * @param flagLocator true if the Locators must be set
 */
void NamingConvention::setNamesAndLocators(Db* dbout,
                                           int iattout_start,
                                           const String& suffix,
                                           int nitems,
                                           bool flagLocator) const
{
  _setNames(dbout,iattout_start,VectorString(),suffix,nitems);

  if (flagLocator)
    setLocators(dbout,iattout_start,1,nitems);
}

/**
 * Naming from a set of variables identified by their names
 * @param names Vector of variable names in Input Db
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param suffix Optional suffix
 * @param nitems Number of items
 * @param flagLocator true if the Locators must be set
 */
void NamingConvention::setNamesAndLocators(const VectorString& names,
                                           Db* dbout,
                                           int iattout_start,
                                           const String& suffix,
                                           int nitems,
                                           bool flagLocator) const
{
  if (iattout_start <= 0) return;
  int nvar = names.size();
  if (nvar <= 0) return;
  _setNames(dbout,iattout_start,names,suffix,nitems);
  if (flagLocator)
    setLocators(dbout,iattout_start,nvar,nitems);
}

/**
 * Naming from one variable identified by its names
 * @param namin variable name in Input Db
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param suffix Optional suffix
 * @param nitems Number of items
 * @param flagLocator true if the Locators must be set
 */
void NamingConvention::setNamesAndLocators(const String& namin,
                                           Db* dbout,
                                           int iattout_start,
                                           const String& suffix,
                                           int nitems,
                                           bool flagLocator) const
{
  if (iattout_start <= 0) return;
  VectorString names;
  names.push_back(namin);
  _setNames(dbout,iattout_start,names,suffix,nitems);
  if (flagLocator)
    setLocators(dbout,iattout_start,1,nitems);
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
 * @param flagLocator true if the Locators must be set
 */
void NamingConvention::setNamesAndLocators(Db *dbin,
                                           ENUM_LOCS locatorInType,
                                           int nvar,
                                           Db* dbout,
                                           int iattout_start,
                                           const String& suffix,
                                           int nitems,
                                           bool flagLocator) const
{
  if (iattout_start <= 0) return;
  VectorString names;
  if (dbin != nullptr)
  {
    names = dbin->getNames(locatorInType);
    if (nvar <= 0) nvar = names.size();
    if (nvar < (int) names.size()) names.resize(nvar);
  }

  _setNames(dbout,iattout_start,names,suffix,nitems);
  if (flagLocator)
    setLocators(dbout,iattout_start,nvar,nitems);
}

/**
 * Naming from a set of variables identified by their attribute indices
 * @param dbin  Pointer to the input Db (kept for symmetry)
 * @param iatts Vector of attribute indices of the variables in Input Db
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param suffix Optional suffix
 * @param nitems Number of items
 * @param flagLocator true if the Locators must be set
 */
void NamingConvention::setNamesAndLocators(Db *dbin,
                                           const VectorInt& iatts,
                                           Db* dbout,
                                           int iattout_start,
                                           const String& suffix,
                                           int nitems,
                                           bool flagLocator) const
{
  if (iattout_start <= 0) return;
  if (dbin == nullptr) return;
  int nvar = iatts.size();
  if (nvar <= 0) return;

  VectorString names;
  for (int ivar = 0; ivar < nvar; ivar++)
    names.push_back(dbin->getName(iatts[ivar]));
  _setNames(dbout,iattout_start,names,suffix,nitems);
  if (flagLocator)
    setLocators(dbout,iattout_start,nvar,nitems);
 }

/**
 * Naming from a set of variables identified by their attribute indices
 * @param dbin  Pointer to the input Db (kept for symmetry)
 * @param iatt  Attribute index of the variablesin Input Db
 * @param dbout Pointer to the output Db
 * @param iattout_start Starting attribute index
 * @param suffix Optional suffix
 * @param nitems Number of items
 * @param flagLocator true if the Locators must be set
 */
void NamingConvention::setNamesAndLocators(Db *dbin,
                                           int iatt,
                                           Db* dbout,
                                           int iattout_start,
                                           const String& suffix,
                                           int nitems,
                                           bool flagLocator) const
{
  if (iattout_start <= 0) return;
  if (dbin == nullptr) return;

  VectorString names;
  names.push_back(dbin->getName(iatt));
  _setNames(dbout,iattout_start,names,suffix,nitems);
  if (flagLocator)
    setLocators(dbout,iattout_start,1,nitems);
 }

void NamingConvention::setLocators(Db *dbout,
                                    int iattout_start,
                                    int nvar,
                                    int nitems) const
{
  if (_locatorType < 0) return;

  if (_flagClean)
  {
    // Erase already existing locators of the same Type
    dbout->clearLocators(_locatorType);
  }

  for (int ecr = 0; ecr < nvar * nitems; ecr++)
    dbout->setLocatorByAttribute(iattout_start + ecr, _locatorType, ecr + 1);
}

void NamingConvention::_setNames(Db *dbout,
                                 int iattout_start,
                                 VectorString names,
                                 const String& suffix,
                                 int nitems) const
{
  int ecr = 0;
  int nvar = (names.empty()) ? 1 : names.size();

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    String local = (names.empty()) ? String() : names[ivar];

    for (int item = 0; item < nitems; item++)
    {
      String locnum = (nitems <= 1) ? String() : std::to_string(item+1);
      String name = concatenateStrings(_delim,_radix, local, suffix, locnum);
      dbout->setName(iattout_start + ecr, name);
      ecr++;
    }
  }
}
