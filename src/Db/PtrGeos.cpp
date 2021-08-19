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
#include "Db/PtrGeos.hpp"
#include "Basic/AStringable.hpp"
#include "geoslib_enum.h"
#include <string.h>

typedef struct
{
  char SREF[LOCAL_SIZE]; /* Name of the Locator */
  int  IREF; /* Unicity of the locator */
  char COMMENT[STRING_LENGTH]; /* Meaning */
} Def_Locator;

static Def_Locator DEF_LOCATOR[] = { { "x",       0, "Coordinate" },
                                     { "z",       0, "Variable" },
                                     { "v",       0, "Variance of measurement error" },
                                     { "f",       0, "External Drift" },
                                     { "g",       0, "Gradient component" },
                                     { "lower",   0,"Lower bound of an inequality" },
                                     { "upper",   0,"Upper bound of an inequality" },
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
                                     { "rkup",    0, "Disc. rank for Upper bound" }
  };

void PtrGeos::clear()
{
  _r.clear();
}

void PtrGeos::erase(int locatorIndex)
{
  _r.erase(_r.begin() + locatorIndex);
}

void PtrGeos::print(int rank, ENUM_LOCS locatorType) const
{
  message("%d - Locator: %s\n", rank + 1, DEF_LOCATOR[locatorType].SREF);
  message("- Attributes = ");
  for (int locatorIndex = 0; locatorIndex < getLocatorNumber(); locatorIndex++)
    message("%2d ", _r[locatorIndex]);
  message("\n");
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
String getLocatorName(ENUM_LOCS locatorType, int locatorIndex)
{
  std::stringstream sstr;
  if (locatorType < 0)
  {
    sstr << "NA";
  }
  else if (! isLocatorTypeValid(locatorType))
  {
    sstr << "NA";
  }
  else
  {
    if (DEF_LOCATOR[locatorType].IREF == 1)
      sstr << DEF_LOCATOR[locatorType].SREF;
    else if (locatorIndex < 0)
      sstr << DEF_LOCATOR[locatorType].SREF;
    else
      sstr << DEF_LOCATOR[locatorType].SREF << locatorIndex+1;
  }
  return sstr.str();
}

/**
 * Check if the Locator type is valid or not
 * Note that the locator type is returned as -1 for non identified locator (such as rank)
 * @param locatorType The locator type to be identified
 * @param unknownValid True if LOC_UNKNOWN is considered as valid
 * @return
 */
bool isLocatorTypeValid(ENUM_LOCS locatorType, bool unknownValid)
{
  int lower_bound = (unknownValid) ? -1 : 0;
  if (locatorType < lower_bound || locatorType >= MAXIMUM_LOC)
  {
    messerr("Error in the Locator Type (%d)",locatorType);
    messerr("It should be defined using the ENUM_LOCS");
    return false;
  }
  return true;
}

int getLocatorTypeFromName(const String& name_type)
{
  for (int i = 0; i < MAXIMUM_LOC; i++)
  {
    unsigned int lng = static_cast<unsigned int> (strlen(DEF_LOCATOR[i].SREF));
    if (name_type.compare(0,lng,DEF_LOCATOR[i].SREF) == 0) return i;
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
int locatorIdentify(String string, ENUM_LOCS *ret_locatorType, int *ret_item, int *ret_mult)
{
  *ret_locatorType   = LOC_UNKNOWN;
  *ret_item     = -1;
  *ret_mult     =  1;
  int  inum  = -1;
  int  found = -1;
  bool mult  =  0;
  for (int i = 0; i < MAXIMUM_LOC && found < 0; i++)
  {
    unsigned int lng = static_cast<unsigned int> (strlen(DEF_LOCATOR[i].SREF));
    if (string.compare(0,lng,DEF_LOCATOR[i].SREF) == 0) found = i;
  }
  if (found < 0) return 1;
  unsigned int lng = static_cast<unsigned int> (strlen(DEF_LOCATOR[found].SREF));
  if (string.size() > lng) inum = atoi(&string[lng]);
  mult = DEF_LOCATOR[found].IREF == 0;
  if (! mult && inum > 1)
  {
    string = "NA";
    return 1;
  }

  /* Returning arguments */

  *ret_locatorType = (ENUM_LOCS) found;
  *ret_item   = MAX(inum-1, 0);
  *ret_mult   = mult;
  return 0;
}

void printLocatorList()
{
  mestitle(0, "List of the available locators");
  for (int i = 0; i < MAXIMUM_LOC; i++)
  {
    if (DEF_LOCATOR[i].IREF)
      message(" %10s %s\n", DEF_LOCATOR[i].SREF, DEF_LOCATOR[i].COMMENT);
    else
      message(" %7s(*) %s\n", DEF_LOCATOR[i].SREF, DEF_LOCATOR[i].COMMENT);
  }
  message("(*) These keywords must be followed by a number\n");
  return;
}

VectorString getLocatorNames()
{
  VectorString strings;
  for (int i = 0; i < MAXIMUM_LOC; i++)
    strings.push_back(DEF_LOCATOR[i].SREF);
  return strings;
}

VectorInt getLocatorMultiples()
{
  VectorInt mult;
  for (int i = 0; i < MAXIMUM_LOC; i++)
    mult.push_back(DEF_LOCATOR[i].IREF);
  return mult;
}
