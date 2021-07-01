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
#include "Basic/AStringable.hpp"
#include "Basic/Vector.hpp"
#include "geoslib_enum.h"

class Db;

class NamingConvention
{
private:
  String   _radix;
  String   _delim;
  ENUM_LOCS _locatorType;
  bool     _flagVarIter;
  bool     _flagClean;

public:
  NamingConvention(String radix = String(),
                   ENUM_LOCS locatorType = LOC_Z,
                   String delim = ".",
                   bool flagvariter = true,
                   bool flagclean = true);
  NamingConvention(const NamingConvention &m);
  NamingConvention& operator=(const NamingConvention &m);
  virtual ~NamingConvention();

  void setNamesAndLocators(Db* dbout,
                           int iattout_start,
                           const String& suffix = String(),
                           int nitems = 1,
                           bool flagLocator = false) const;
  void setNamesAndLocators(const VectorString& names,
                           Db* dbout,
                           int iattout_start,
                           const String& suffix = String(),
                           int nitems = 1,
                           bool flagLocator = false) const;
  void setNamesAndLocators(const String& namin,
                           Db* dbout,
                           int iattout_start,
                           const String& suffix = String(),
                           int nitems = 1,
                           bool flagLocator = false) const;
  void setNamesAndLocators(Db *dbin,
                           ENUM_LOCS locatorInType,
                           int nvar,
                           Db* dbout,
                           int iattout_start,
                           const String& suffix = String(),
                           int nitems = 1,
                           bool flagLocator = false) const;
  void setNamesAndLocators(Db *dbin,
                           const VectorInt& iatts,
                           Db* dbout,
                           int iattout_start,
                           const String& suffix = String(),
                           int nitems = 1,
                           bool flagLocator = false) const;
  void setNamesAndLocators(Db *dbin,
                           int iatt,
                           Db* dbout,
                           int iattout_start,
                           const String& suffix = String(),
                           int nitems = 1,
                           bool flagLocator = false) const;
  void setLocators(Db *dbout, int iattout_start, int nvar, int nitems = 1) const;

  void setDelim(const String& delim)    { _delim = delim; }
  void setLocatorType(ENUM_LOCS locatorType)  { _locatorType = locatorType; }
  void setRadix(const String& radix)    { _radix = radix; }
  void setFlagClean(bool flagClean)     { _flagClean = flagClean; }
  void setFlagVarIter(bool flagVarIter) { _flagVarIter = flagVarIter; }

private:
  void _setNames(Db *dbout,
                 int iattout_start,
                 VectorString names,
                 const String& suffix,
                 int nitems) const;
};
