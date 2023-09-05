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

class Db;

/**
 * Naming Convention facility
 * This class describes the way variables created within the current procedure
 * will be named afterwards and will possibly be assigned a locator.
 *
 * The generic name is generated as follows:
 *      'prefix'.'varname'.'qualifier'.'rank'
 *
 * - prefix: string provided in the constructor of this class
 * - varname: name of the (input) variable on which the procedure is performed
 * - qualifier: type of element stored in the variable
 * - rank: rank of the output variable (if several variables of the same type are generated)
 *
 * The choice of the 'prefix' is done by the user when launching the procedure
 * the other parameters are usually defined within the procedure.
 *
 * For example, when running 'kriging' function with several variables defined
 * in the input Db - say "Pb" abd "Zn" (they are assigned a Z-locator),
 * using the following command:
 *    kriging( ... namconv = NamingConvention("MyPrefix") )
 *
 * Then the 'kriging' procedure generates variables such as:
 *
 * - MyPrefix.Pb.estim (estimation of Pb by CoKriging)
 * - MyPrefix.Zn.estim (estimation of Zn by CoKriging)
 * - MyPrefix.Pb.stdev (St. Dev. of estimation error of Pb by CoKriging)
 * - MyPrefix.Zn.stdev (St. Dev. of estimation error of Zn by CoKriging)
 *
 * Ultimately, the newly created variables are assigned a locator.
 */
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
  VectorString _createNames(const VectorString &names,
                           const String &qualifier = String(),
                           int nitems = 1) const;

private:
  String _prefix; //!< String used as 'prefix'
  String _delim; //!< Character used as the 'delimitor' between different parts of the names
  bool   _flagVarname; //!< When TRUE, add the 'variable name'
  bool   _flagQualifier; //!< When TRUE, add the 'qualifier'
  bool   _flagLocator; //!< When TRUE, assign a locator to the newly created variables
  ELoc   _locatorOutType; //!< Type of locator assigned (if 'flagLocator' is TRUE)
  bool   _cleanSameLocator; //!< Clean variables with the same locator beforehand
};
