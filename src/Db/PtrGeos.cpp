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
#include "Db/PtrGeos.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"
#include "Basic/Utilities.hpp"

#include "geoslib_enum.h"

#include <string.h>
#include <sstream>

typedef struct
{
  char SREF[LOCAL_SIZE]; /* Name of the Locator */
  int  IREF; /* Unicity of the locator */
  char COMMENT[STRING_LENGTH]; /* Meaning */
} Def_Locator;

// TODO : DEF_LOCATOR static table refactoring. Sync with ELoc
static Def_Locator DEF_LOCATOR[] = { { "x",       0, "Coordinate" },
                                     { "z",       0, "Variable" },
                                     { "v",       0, "Variance of measurement error" },
                                     { "f",       0, "External Drift" },
                                     { "g",       0, "Gradient component" },
                                     { "lower",   0, "Lower bound of an inequality" },
                                     { "upper",   0, "Upper bound of an inequality" },
                                     { "p",       0, "Proportion" },
                                     { "w",       1, "Weight" },
                                     { "code",    1, "Code" },
                                     { "sel",     1, "Selection" },
                                     { "dom",     1, "Domain selection" },
                                     { "dblk",    0, "Block Extension" },
                                     { "adir",    1, "Dip direction angle" },
                                     { "adip",    1, "Dip angle" },
                                     { "size",    1, "Object height" },
                                     { "bu",      1, "Fault UP termination" },
                                     { "bd",      1, "Fault DOWN termination" },
                                     { "time",    0, "Time variable" },
                                     { "layer",   1, "Layer rank" },
                                     { "nostat",  0, "Non-stationary parameter" },
                                     { "tangent", 0, "Tangent" },
                                     { "ncsimu",  0, "Non-conditional simulation" },
                                     { "facies",  0, "Facies simulated" },
                                     { "gausfac", 0, "Gaussian value for Facies" },
                                     { "date",    1, "Date" },
                                     { "rklow",   0, "Disc. rank for Lower bound" },
                                     { "rkup",    0, "Disc. rank for Upper bound" },
                                     { "sum",     0, "Constraints on the sum" }
  };

void PtrGeos::clear()
{
  _r.clear();
}

void PtrGeos::erase(int locatorIndex)
{
  _r.erase(_r.begin() + locatorIndex);
}

int PtrGeos::findUIDInLocator(int iuid) const
{
  for (int locatorIndex = 0; locatorIndex < getLocatorNumber(); locatorIndex++)
    if (getLocatorByIndex(locatorIndex) == iuid) return (locatorIndex);
  return -1;
}

String PtrGeos::dumpLocator(int rank, const ELoc& locatorType) const
{
  std::stringstream sstr;

  int i = locatorType.getValue();
  sstr << rank+1 << " - Locator: " << DEF_LOCATOR[i].SREF << std::endl;
  sstr << "- Attributes = ";
  for (int locatorIndex = 0; locatorIndex < getLocatorNumber(); locatorIndex++)
    sstr << _r[locatorIndex] << " ";
  sstr << std::endl;

  return sstr.str();
}

bool PtrGeos::isLocatorIndexValid(int locatorIndex) const
{
  if (locatorIndex < 0 || locatorIndex >= getLocatorNumber())
  {
    mesArg("Locator Index", locatorIndex, getLocatorNumber());
    return false;
  }
  return true;
}

/**
 * Return the name of Locator
 * @param locatorType    Type of the Locator (can be negative for 'Rank')
 * @param locatorIndex   Rank within the locator starting from 1 (can be <0 for the keyword only)
 * @return
 */
String getLocatorName(const ELoc& locatorType, int locatorIndex)
{
  std::stringstream sstr;
  if (locatorType == ELoc::UNKNOWN)
  {
    sstr << STRING_NA;
  }
  else if (! isLocatorTypeValid(locatorType))
  {
    sstr << STRING_NA;
  }
  else
  {
    int i = locatorType.getValue();
    if (DEF_LOCATOR[i].IREF == 1)
      sstr << DEF_LOCATOR[i].SREF;
    else if (locatorIndex < 0)
      sstr << DEF_LOCATOR[i].SREF;
    else
      sstr << DEF_LOCATOR[i].SREF << locatorIndex+1;
  }
  return sstr.str();
}

/**
 * Check if the Locator type is valid or not
 * Note that the locator type is returned as -1 for non identified locator (such as rank)
 * @param locatorType The locator type to be identified
 * @param unknownValid True if ELoc::UNKNOWN is considered as valid
 * @return
 */
bool isLocatorTypeValid(const ELoc& locatorType, bool unknownValid)
{
  if (unknownValid) return true;
  if (locatorType == ELoc::UNKNOWN)
  {
    messerr("Locator Type must not be UNKNOWN");
    return false;
  }
  return true;
}

int getLocatorTypeFromName(const String& name_type)
{
  auto it = ELoc::getIterator();
  while (it.hasNext())
  {
    if (*it != ELoc::UNKNOWN)
    {
      int i = it.getValue();
      unsigned int lng = static_cast<unsigned int> (strlen(DEF_LOCATOR[i].SREF));
      if (name_type.compare(0,lng,DEF_LOCATOR[i].SREF) == 0) return i;
    }
    it.toNext();
  }
  return -1;
}

/**
 * Given a locator string, extract its characteristics
 * @param string     Locator string
 * @param ret_locatorType Resulting Locator type
 * @param ret_item   Resulting Locator rank (starting from 0)
 * @param ret_mult   Resulting Locator multiplicity (1: unique; 0: multiple)
 * @return Error code
 */
int locatorIdentify(String string, ELoc* ret_locatorType, int* ret_item, int* ret_mult)
{
  *ret_locatorType   = ELoc::UNKNOWN;
  *ret_item     = -1;
  *ret_mult     =  1;
  int  inum  = -1;
  int  found = -1;
  bool mult  =  0;

  // Transform the input argument into lower case for comparison
  String string_loc = string;
  toLower(string_loc);

  auto it = ELoc::getIterator();
  while (it.hasNext() && found < 0)
  {
    if (*it != ELoc::UNKNOWN)
    {
      int i = it.getValue();
      unsigned int lng = static_cast<unsigned int> (strlen(DEF_LOCATOR[i].SREF));
      if (string_loc.compare(0,lng,DEF_LOCATOR[i].SREF) == 0) found = i;
    }
    it.toNext();
  }
  if (found < 0)
  {
    // The locator has not been matched. It is returned as UNKNOWN
    *ret_locatorType = ELoc::UNKNOWN;
    *ret_item   = 0;
    *ret_mult   = 0;
    return 0;
  }

  // Decode the remaining characteristics
  unsigned int lng = static_cast<unsigned int> (strlen(DEF_LOCATOR[found].SREF));
  if (string_loc.size() > lng) inum = atoi(&string_loc[lng]);
  mult = (DEF_LOCATOR[found].IREF == 0);
  if (! mult && inum > 1)
  {
    // The locator has an index larger than 1 but the Locator should be Unique. Error
    string = STRING_NA;
    return 1;
  }

  /* Returning arguments */

  *ret_locatorType = ELoc::fromValue(found);
  *ret_item   = MAX(inum-1, 0);
  *ret_mult   = mult;
  return 0;
}

void printLocatorList()
{
  mestitle(0, "List of the available locators");
  auto it = ELoc::getIterator();
  while (it.hasNext())
  {
    if (*it != ELoc::UNKNOWN)
    {
      int i = it.getValue();
      if (DEF_LOCATOR[i].IREF)
        message(" %10s %s\n", DEF_LOCATOR[i].SREF, DEF_LOCATOR[i].COMMENT);
      else
        message(" %7s(*) %s\n", DEF_LOCATOR[i].SREF, DEF_LOCATOR[i].COMMENT);
    }
    it.toNext();
  }
  message("(*) These keywords must be followed by a number\n");
  return;
}

VectorString getLocatorNames()
{
  VectorString strings;
  auto it = ELoc::getIterator();
  while (it.hasNext())
  {
    if (*it != ELoc::UNKNOWN)
    {
      int i = it.getValue();
      strings.push_back(DEF_LOCATOR[i].SREF);
    }
    it.toNext();
  }
  return strings;
}

VectorInt getLocatorMultiples()
{
  VectorInt mult;
  auto it = ELoc::getIterator();
  while (it.hasNext())
  {
    if (*it != ELoc::UNKNOWN)
    {
      int i = it.getValue();
      mult.push_back(DEF_LOCATOR[i].IREF);
    }
    it.toNext();
  }
  return mult;
}
