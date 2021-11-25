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

#include "gstlearn_export.hpp"
#include "Db/ELoc.hpp"
#include "Basic/Vector.hpp"

/**
 * Gives the Rank of the Attribute for:
 * - a given pointer type
 * - a given locatorIndex rank
 * The dimension of the internal vector is equal to the number of items for a
 * given pointer type
 */
class GSTLEARN_EXPORT PtrGeos {
public:
  VectorInt _r;    /* Rank of the attribute */

  bool isLocatorIndexValid(int locatorIndex) const;
  int  getLocatorByIndex(int locatorIndex) const { return _r[locatorIndex]; }
  void setLocatorByIndex(int locatorIndex, int value) { _r[locatorIndex] = value; }
  int  getLocatorNumber() const { return static_cast<int>(_r.size()); }
  void erase(int locatorIndex);
  void clear();
  void print(int rank, const ELoc& locatorType) const;
  void resize(int count) { _r.resize(count,0); }
};

GSTLEARN_EXPORT int    getLocatorTypeFromName(const String& name_type);
GSTLEARN_EXPORT int    locatorIdentify(String string, ELoc* locatorType, int* locatorIndex, int *mult);
GSTLEARN_EXPORT bool   isLocatorTypeValid(const ELoc& locatorType, bool unknownValid = false);
GSTLEARN_EXPORT String getLocatorName(const ELoc& locatorType, int locatorIndex=1);
GSTLEARN_EXPORT void   printLocatorList();
GSTLEARN_EXPORT VectorString getLocatorNames();
GSTLEARN_EXPORT VectorInt    getLocatorMultiples();
