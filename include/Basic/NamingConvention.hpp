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
   * @brief Naming Convention facility.
   *
   * This class describes the way variables created within the current procedure
   * will be named afterwards and will possibly be assigned a locator.
   *
   * The generic name is generated as follows:
   *      'prefix'.'varname'.'qualifier'
   *
   * - prefix: string provided in the constructor of this class
   * - varname: name of the (input) variable on which the procedure is performed
   * - qualifier: type of element stored in the variable
   *
   * The choice of the 'prefix' is done by the user when launching the procedure;
   * the other parameters are usually defined within the calling procedure.
   *
   * For example, when running 'kriging' function with several variables defined
   * in the input Db - say "Pb" and "Zn" (they are assigned a Z-locator),
   * using the following command:
   *    kriging( ... namconv = NamingConvention("MyPrefix") )
   *
   * Then the kriging procedure generates variables such as:
   * - MyPrefix.Pb.estim (estimation of Pb by CoKriging)
   * - MyPrefix.Zn.estim (estimation of Zn by CoKriging)
   * - MyPrefix.Pb.stdev (St. Dev. of estimation error of Pb by CoKriging)
   * - MyPrefix.Zn.stdev (St. Dev. of estimation error of Zn by CoKriging)
   *
   * ...
   *
   * Ultimately, the newly created variables are assigned a locator.
   *
   * Note: the related method getNameEncoded provides a static way to retrieve
   * the variable name based on the same convention (see comments).
   */
  class GSTLEARN_EXPORT NamingConvention: public AStringable
  {
  public:
    /**
     * @brief Constructor.
     *
     * @param prefix Prefix used for naming the output variables.
     * @param flag_varname If true, the variable name is included in the
     * generated name.
     * @param flag_qualifier If true, the qualifier is included in the
     * generated name.
     * @param flag_locator If true, a locator is assigned to the output
     * variables.
     * @param locatorOutType Type of locator assigned to the output variables.
     * @param delim Delimiter used to separate the components of the generated
     * name.
     * @param cleanSameLocator If true, variables with the same locator are
     * cleaned beforehand.
     */
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

    void setOutput(
      const VectorString& names,
      Id nvar,
      Db* dbout,
      Id iattout_start,
      const String& qualifier = "",
      Id nitems = 1,
      bool flagSetLocator = true,
      Id locatorShift = 0) const;

    void setLocatorOutType(const ELoc& l) { _locatorOutType = l; }

    void setLocators(
      Db* dbout,
      Id iattout_start,
      Id nvar,
      Id nitems = 1,
      Id locatorShift = 0) const;

    static void Naming_Old_Style(bool status);

    static String getNameEncoded(
      const String& prefix,
      const Db* db = nullptr,
      Id ivar = 0,
      Id nvar = 0,
      Id isimu = 0,
      Id nbsimu = 0,
      const String& extension = "",
      const String& delim = ".");

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

    static Id _getNameCount(const VectorString& names, Id nvar);

  private:
    String _prefix; //!< String used as 'prefix'
    String _delim; //!< Character used as 'delimiter'
    bool _flagVarname; //!< When TRUE, add the 'variable name'
    bool _flagQualifier; //!< When TRUE, add the 'qualifier'
    bool _flagLocator; //!< When TRUE, assign a locator to the new variables
    ELoc _locatorOutType; //!< Type of locator assigned ('flagLocator' is TRUE)
    bool _cleanSameLocator; //!< Clean variables with same locator beforehand
  };

  // typedef NamingConvention NC;
  class NC: public NamingConvention
  {
  };

} // namespace gstlrn
