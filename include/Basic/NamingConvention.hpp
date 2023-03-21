/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Enum/ELoc.hpp"

class Db;

class GSTLEARN_EXPORT NamingConvention
{
public:
  NamingConvention(String prefix = String(),
                   bool flag_varname = true,
                   bool flag_qualifier = true,
                   bool flag_locator = true,
                   const ELoc& locatorOutType = ELoc::fromKey("Z"),
                   String delim = ".",
                   bool cleanSameLocator = true);
  NamingConvention(const NamingConvention &m);
  NamingConvention& operator=(const NamingConvention &m);
  virtual ~NamingConvention();

  static NamingConvention* create(String prefix = String(),
                                  bool flag_varname = true,
                                  bool flag_qualifier = true,
                                  bool flag_locator = true,
                                  const ELoc &locatorOutType = ELoc::fromKey("Z"),
                                  String delim = ".",
                                  bool cleanSameLocator = true);

  VectorString createNames(const VectorString &names,
                           const String &qualifier = String(),
                           int nitems = 1) const;

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
  bool   _flagLocator;
  ELoc   _locatorOutType;
  bool   _cleanSameLocator;
};
