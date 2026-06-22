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

#include "Basic/AStringable.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Enum/ELoc.hpp"

namespace gstlrn
{
  class Db;

  /**
   * Naming Convention facility
   * This class describes the way variables created within the current procedure
   * will be named afterwards and will possibly be assigned a locator.
   *
   * The generic name is generated as follows:
   *      'prefix'.'varname'.'qualifier'|'rank'
   *
   * - prefix: string provided in the constructor of this class
   * - varname: name of the (input) variable on which the procedure is performed
   * - qualifier: type of element stored in the variable
   * - rank: rank of the output variable (if several simulations are generated)
   *
   * The choice of the 'prefix' is done by the user when launching the procedure
   * the other parameters are usually defined within the calling procedure.
   *
   * For example, when running 'kriging' function with several variables defined
   * in the input Db - say "Pb" abd "Zn" (they are assigned a Z-locator),
   * using the following command:
   *    kriging( ... namconv = NamingConvention("MyPrefix") )
   *
   * Then the kriging procedure generates variables such as:
   * - MyPrefix.Pb.estim (estimation of Pb by CoKriging)
   * - MyPrefix.Zn.estim (estimation of Zn by CoKriging)
   * - MyPrefix.Pb.stdev (St. Dev. of estimation error of Pb by CoKriging)
   * - MyPrefix.Zn.stdev (St. Dev. of estimation error of Zn by CoKriging)
   *
   * Then the non-conditional simulation procedure generates variables such as:
   * - MyPrefix.S1 (for first simulation)
   * - MyPrefix.S2 (for second simulation)
   * ...
   *
   *  Then the conditional simulation procedure generates variables such as:
   * - MyPrefix.var.S1 (for first simulation)
   * - MyPrefix.var.S2 (for second simulation)
   * ...
   *
   * For multivariate simulations, the new setNamesAndLocatorsForSimulations method
   * provides consistent naming with explicit Variable and Simulation indicators:
   *
   * Non-conditional multivariate simulations (e.g., 2 variables, 2 simulations):
   * - MyPrefix.V1.S1, MyPrefix.V1.S2, MyPrefix.V2.S1, MyPrefix.V2.S2
   *
   * Conditional multivariate simulations (e.g., variables Fe and Al, 2 simulations):
   * - MyPrefix.Fe.S1, MyPrefix.Fe.S2, MyPrefix.Al.S1, MyPrefix.Al.S2
   *
   * Ultimately, the newly created variables are assigned a locator.
   *
   * Note: the related method getNameEncoded provides a static way to retrieve
   * the variable name based on the same convention (see comments).
   */
  class GSTLEARN_EXPORT NamingConvention: public AStringable
  {
  public:
    NamingConvention(
      const String& prefix = "",
      bool flag_varname = true,
      bool flag_qualifier = true,
      bool flag_locator = true,
      const ELoc& locatorOutType = ELoc::fromKey("Z"),
      const String& delim = ".",
      bool cleanSameLocator = true);
    NamingConvention(const NamingConvention& m);
    NamingConvention& operator=(const NamingConvention& m);
    virtual ~NamingConvention();

    /// AStringable Interface
    String toString(const AStringFormat* strfmt = nullptr) const override;

    static NamingConvention* create(
      const String& prefix = "",
      bool flag_varname = true,
      bool flag_qualifier = true,
      bool flag_locator = true,
      const ELoc& locatorOutType = ELoc::fromKey("Z"),
      const String& delim = ".",
      bool cleanSameLocator = true);

    void setNamesAndLocators(
      Db* dbout,
      Id iattout_start,
      const String& qualifier = "",
      Id nitems = 1,
      bool flagSetLocator = true,
      Id locatorShift = 0) const;
    void setNamesAndLocators(
      const VectorString& names,
      Db* dbout,
      Id iattout_start,
      const String& qualifier = "",
      Id nitems = 1,
      bool flagSetLocator = true,
      Id locatorShift = 0) const;
    void setNamesAndLocators(
      Db* dbout,
      Id iattout_start,
      const VectorString& names,
      bool flagSetLocator = true,
      Id locatorShift = 0) const;
    void setNamesAndLocators(
      const String& namin,
      Db* dbout,
      Id iattout_start,
      const String& qualifier = "",
      Id nitems = 1,
      bool flagSetLocator = true,
      Id locatorShift = 0) const;
    void setNamesAndLocators(
      const Db* dbin,
      const VectorString& names,
      const ELoc& locatorInType,
      Id nvar,
      Db* dbout,
      Id iattout_start,
      const String& qualifier = "",
      Id nitems = 1,
      bool flagSetLocator = true,
      Id locatorShift = 0) const;
    void setNamesAndLocators(
      const Db* dbin,
      const VectorInt& iatts,
      Db* dbout,
      Id iattout_start,
      const String& qualifier = "",
      Id nitems = 1,
      bool flagSetLocator = true,
      Id locatorShift = 0) const;
    void setNamesAndLocators(
      const Db* dbin,
      Id iatt,
      Db* dbout,
      Id iattout_start,
      const String& qualifier = "",
      Id nitems = 1,
      bool flagSetLocator = true,
      Id locatorShift = 0) const;
    void setNamesAndLocatorsForSimulations(
      const Db* dbin,
      const VectorString& names,
      const ELoc& locatorInType,
      Id nvar,
      Db* dbout,
      Id iattout_start,
      Id nbsimu,
      bool flagSimuFirst = true,
      bool flagSetLocator = true,
      Id locatorShift = 0) const;

    void setDelim(const String& delim) { _delim = delim; }

    void setLocatorOutType(const ELoc& l) { _locatorOutType = l; }

    void setPrefix(const String& prefix) { _prefix = prefix; }

    void setFlagClean(bool cleanSameLocator)
    {
      _cleanSameLocator = cleanSameLocator;
    }

    void setLocators(
      Db* dbout,
      Id iattout_start,
      Id nvar,
      Id nitems = 1,
      Id locatorShift = 0) const;

    bool isFlagQualifier() const { return _flagQualifier; }

    bool isFlagVarname() const { return _flagVarname; }

    String getPrefix() const { return _prefix; }

    String getDelim() const { return _delim; }

    static String getNameEncoded(
      const String& prefix,
      const Db* db = nullptr,
      Id ivar = 0,
      Id nvar = 0,
      Id isimu = 0,
      Id nbsimu = 0,
      const String& extension = "",
      const String& delim = ".");

    static void Naming_Old_Style(bool status);

  private:
    void _setNames(
      Db* dbout,
      Id iattout_start,
      const VectorString& names,
      Id nvar,
      const String& qualifier,
      Id nitems) const;
    VectorString _createNames(
      const VectorString& names,
      Id nvar,
      const String& qualifier = "",
      Id nitems = 1) const;
    VectorString _createSimulationNames(
      const VectorString& names,
      Id nvar,
      Id nbsimu,
      bool flagSimuFirst) const;
    static Id _getNameCount(const VectorString& names, Id nvar);

  private:
    String _prefix; //!< String used as 'prefix'
    String
      _delim; //!< Character used as the 'delimitor' between different parts of the names
    bool _flagVarname; //!< When TRUE, add the 'variable name'
    bool _flagQualifier; //!< When TRUE, add the 'qualifier'
    bool
      _flagLocator; //!< When TRUE, assign a locator to the newly created variables
    ELoc
      _locatorOutType; //!< Type of locator assigned (if 'flagLocator' is TRUE)
    bool
      _cleanSameLocator; //!< Clean variables with the same locator beforehand
  };

  // typedef NamingConvention NC;
  class NC: public NamingConvention
  {
  };

} // namespace gstlrn
