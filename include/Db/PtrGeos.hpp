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
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Enum/ELoc.hpp"

namespace gstlrn
{

/**
 * Gives the Rank of the Attribute for:
 * - a given pointer type
 * - a given locatorIndex rank
 *
 * The dimension of the internal vector is equal to the number of items for a
 * given pointer type
 */
class GSTLEARN_EXPORT PtrGeos {
public:
  bool isLocatorIndexValid(Id locatorIndex) const;
  Id  getLocatorByIndex(Id locatorIndex) const { return _r[locatorIndex]; }
  void setLocatorByIndex(Id locatorIndex, Id value) { _r[locatorIndex] = value; }
  Id  getNLoc() const { return static_cast<Id>(_r.size()); }
  bool hasLocator() const { return ! _r.empty(); }
  Id  findUIDInLocator(Id iuid) const;
  void erase(Id locatorIndex);
  void clear();
  void resize(Id count) { _r.resize(count,0); }
  String dumpLocator(Id rank, const ELoc& locatorType) const;

private:
  VectorInt _r;    /* Rank of the attribute */
};

GSTLEARN_EXPORT Id getLocatorTypeFromName(const String& name_type);
GSTLEARN_EXPORT Id locatorIdentify(String string,
                                    ELoc* ret_locatorType,
                                    Id* ret_locatorIndex,
                                    Id* ret_mult);
GSTLEARN_EXPORT bool isLocatorTypeValid(const ELoc& locatorType,
                                        bool unknownValid = false);
GSTLEARN_EXPORT String getLocatorName(const ELoc& locatorType, Id locatorIndex=1);
GSTLEARN_EXPORT void   printLocatorList();
GSTLEARN_EXPORT VectorString getLocatorNames();
GSTLEARN_EXPORT VectorInt    getLocatorMultiples();
}
