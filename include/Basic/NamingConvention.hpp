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
#include "geoslib_define.h"

#include "Db/ELoc.hpp"

class Db;

class GSTLEARN_EXPORT NamingConvention
{
public:
  NamingConvention(String prefix = String(),
                   bool flag_varname = true,
                   bool flag_qualifier = true,
                   const ELoc& locatorOutType = ELoc::Z,
                   String delim = ".",
                   bool cleanSameLocator = false);
  NamingConvention(const NamingConvention &m);
  NamingConvention& operator=(const NamingConvention &m);
  virtual ~NamingConvention();

  void setNamesAndLocators(Db* dbout,
                           int iattout_start,
                           const String& qualifier = String(),
                           int nitems = 1,
                           bool flagSetLocator = true,
                           int locatorShift = 0) const;
  void setNamesAndLocators(const VectorString& names,
                           Db* dbout,
                           int iattout_start,
                           const String& qualifier = String(),
                           int nitems = 1,
                           bool flagSetLocator = true,
                           int locatorShift = 0) const;
  void setNamesAndLocators(Db* dbout,
                           int iattout_start,
                           const VectorString& names,
                           bool flagSetLocator = true,
                           int locatorShift = 0) const;
  void setNamesAndLocators(const String& namin,
                           Db* dbout,
                           int iattout_start,
                           const String& qualifier = String(),
                           int nitems = 1,
                           bool flagSetLocator = true,
                           int locatorShift = 0) const;
  void setNamesAndLocators(const Db *dbin,
                           const ELoc& locatorInType,
                           int nvar,
                           Db* dbout,
                           int iattout_start,
                           const String& qualifier = String(),
                           int nitems = 1,
                           bool flagSetLocator = true,
                           int locatorShift = 0) const;
  void setNamesAndLocators(const Db *dbin,
                           const VectorInt& iatts,
                           Db* dbout,
                           int iattout_start,
                           const String& qualifier = String(),
                           int nitems = 1,
                           bool flagSetLocator = true,
                           int locatorShift = 0) const;
  void setNamesAndLocators(const Db *dbin,
                           int iatt,
                           Db* dbout,
                           int iattout_start,
                           const String& qualifier = String(),
                           int nitems = 1,
                           bool flagSetLocator = true,
                           int locatorShift = 0) const;

  void setDelim(const String& delim)    { _delim = delim; }
  void setLocatorOutType(const ELoc& l) { _locatorOutType = l; }
  void setPrefix(const String& prefix)    { _prefix = prefix; }
  void setFlagClean(bool cleanSameLocator) { _cleanSameLocator = cleanSameLocator; }
  void setLocators(Db *dbout,
                   int iattout_start,
                   int nvar,
                   int nitems = 1,
                   int locatorShift = 0) const;

  bool isFlagQualifier() const { return _flagQualifier; }
  bool isFlagVarname() const { return _flagVarname; }

private:
  void _setNames(Db *dbout,
                 int iattout_start,
                 const VectorString& names,
                 const String& qualifier,
                 int nitems) const;

private:
  String _prefix;
  String _delim;
  bool   _flagVarname;
  bool   _flagQualifier;
  ELoc   _locatorOutType;
  bool   _cleanSameLocator;
};
