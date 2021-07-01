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
#pragma once

#include "geoslib_enum.h"
#include "Basic/Vector.hpp"

/**
 * Gives the Rank of the Attribute for:
 * - a given pointer type
 * - a given locatorIndex rank
 * The dimension of the internal vector is equal to the number of items for a
 * given pointer type
 */
class PtrGeos {
public:
  VectorInt _r;                  /* Rank of the attribute */

  bool isLocatorIndexValid(int locatorIndex) const;
  int  getLocatorByIndex(int locatorIndex) const { return _r[locatorIndex]; }
  void setLocatorByIndex(int locatorIndex, int value) { _r[locatorIndex] = value; }
  int  getLocatorNumber() const { return _r.size(); }
  void erase(int locatorIndex);
  void clear();
  void print(int rank, ENUM_LOCS locatorType) const;
  void resize(int count) { _r.resize(count,0); }
};

int    getLocatorTypeFromName(const String& name_type);
int    locatorIdentify(String string, ENUM_LOCS *locatorType, int *locatorIndex, int *mult);
bool   isLocatorTypeValid(ENUM_LOCS locatorType);
String getLocatorName(ENUM_LOCS locatorType, int locatorIndex=1);
void   printLocatorList();
VectorString getLocatorNames();
VectorInt    getLocatorMultiples();
