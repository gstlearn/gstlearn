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
#include "Basic/NamingConvention.hpp"
#include "Basic/String.hpp"
#include "Db/Db.hpp"

#include <string>

namespace gstlrn
{
  // Default value for the Style used for Variable encoding
  bool Old_Style = false;

  /**
   * @brief Activate or deactivate the old naming convention.
   *
   * This method is provided for compatibility with the historical naming
   * convention.
   *
   * @param status If true, activate the old naming convention.
   */
  void NamingConvention::Naming_Old_Style(bool status)
  {
    Old_Style = status;
  }

  NamingConvention::NamingConvention(
    const String& prefix,
    bool flag_varname,
    bool flag_qualifier,
    bool flag_locator,
    const ELoc& locatorOutType,
    const String& delim,
    bool cleanSameLocator)
    : AStringable()
    , _prefix(prefix)
    , _delim(delim)
    , _flagVarname(flag_varname)
    , _flagQualifier(flag_qualifier)
    , _flagLocator(flag_locator)
    , _locatorOutType(locatorOutType)
    , _cleanSameLocator(cleanSameLocator)
  {
  }

  NamingConvention::NamingConvention(const NamingConvention& m)
    : AStringable(m)
    , _prefix(m._prefix)
    , _delim(m._delim)
    , _flagVarname(m._flagVarname)
    , _flagQualifier(m._flagQualifier)
    , _flagLocator(m._flagLocator)
    , _locatorOutType(m._locatorOutType)
    , _cleanSameLocator(m._cleanSameLocator)
  {
  }

  NamingConvention& NamingConvention::operator=(const NamingConvention& m)
  {
    if (this != &m)
    {
      AStringable::operator=(m);
      _prefix = m._prefix;
      _locatorOutType = m._locatorOutType;
      _flagVarname = m._flagVarname;
      _flagQualifier = m._flagQualifier;
      _flagLocator = m._flagLocator;
      _delim = m._delim;
      _cleanSameLocator = m._cleanSameLocator;
    }
    return *this;
  }

  NamingConvention::~NamingConvention() {}

  /**
   * @brief Create a NamingConvention object.
   *
   * This static method is a convenience function for creating a
   * NamingConvention object with the specified naming options.
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
   * @return Pointer to the newly created NamingConvention object.
   */
  NamingConvention* NamingConvention::create(
    const String& prefix,
    bool flag_varname,
    bool flag_qualifier,
    bool flag_locator,
    const ELoc& locatorOutType,
    const String& delim,
    bool cleanSameLocator)
  {
    return new NamingConvention(
      prefix, flag_varname, flag_qualifier, flag_locator, locatorOutType, delim,
      cleanSameLocator);
  }

  /**
   * @brief Generate names for output variables.
   *
   * The generated names are assigned to the output variables starting at
   * the specified attribute index. Depending on the options of the
   * NamingConvention object, the input variable name and qualifier can be
   * included in the generated names.
   *
   * The output variables can also be assigned the configured locator.
   *
   * @param names Names of the input variables.
   * @param nvar Number of variables.
   * @param dbout Output Db containing the variables to be named.
   * @param iattout_start Index of the first output variable.
   * @param qualifier Qualifier describing the output variables.
   * @param nitems Number of items generated for each variable.
   * @param flagSetLocator If true, assign the configured locator to the
   * output variables.
   * @param locatorShift Shift applied when assigning the locator.
   */
  void NamingConvention::setOutput(
    const VectorString& names,
    Id nvar,
    Db* dbout,
    Id iattout_start,
    const String& qualifier,
    Id nitems,
    bool flagSetLocator,
    Id locatorShift) const
  {
    // If (starting) UID is negative, simply skip the naming process
    // (this is a common case when no output variable is created)
    if (iattout_start < 0) return;

    // The corresponding column may have been deleted in the meantime,
    // so we check that the UID is still valid
    if (!dbout->isUIDValid(iattout_start)) return;

    auto nameloc = names;
    if (nameloc.empty())
    {
      // 'names' is not provided, 'nvar' prevails (if not defined, it is set to 1)
      if (nvar <= 0) nvar = 1;
    }
    else
    {
      // 'names' is provided
      auto namesize = static_cast<Id>(nameloc.size());
      if (nvar <= 0)
      {
        // If 'nvar' is not defined, argument 'names' prevails
        nvar = namesize;
      }
      else
      {
        // 'names' and 'nvar' are both defined: 'nvar' prevails
        if (namesize == 1 && nvar > 1)
        {
          // Particular case where 'nvar' > 1 but 'names' contains a single name: the name is expanded
          nameloc = generateMultipleNames(names[0], nvar);
        }
        else
        {
          // 'names' and 'nvar' are both defined: 'names' is reset to 'nvar' if needed
          nameloc.resize(nvar);
        }
      }
    }

    _setNames(dbout, iattout_start, nameloc, nvar, qualifier, nitems);

    if (flagSetLocator)
      setLocators(dbout, iattout_start, nvar, nitems, locatorShift);
  }

  /**
   * @brief Generate names for simulation output variables.
   *
   * This method generates names using both the variable and simulation
   * indices. The order of these two indices is controlled by
   * `flagSimuFirst`.
   *
   * For example, for two variables and two simulations, the generated names
   * can be:
   * - V1.S1, V1.S2, V2.S1, V2.S2 when `flagSimuFirst` is false;
   * - S1.V1, S1.V2, S2.V1, S2.V2 when `flagSimuFirst` is true.
   *
   * @param names Names of the input variables.
   * @param nvar Number of variables.
   * @param dbout Output Db containing the simulation variables.
   * @param iattout_start Index of the first output variable.
   * @param nbsimu Number of simulations.
   * @param flagSimuFirst If true, the simulation index is placed before the
   * variable index.
   * @param flagSetLocator If true, assign the configured locator to the
   * output variables.
   * @param locatorShift Shift applied when assigning the locator.
   */
  void NamingConvention::setOutputForSimulations(
    const VectorString& names,
    Id nvar,
    Db* dbout,
    Id iattout_start,
    Id nbsimu,
    bool flagSimuFirst,
    bool flagSetLocator,
    Id locatorShift) const
  {
    // If (starting) UID is negative, simply skip the naming process
    // (this is a common case when no output variable is created)
    if (iattout_start < 0) return;

    // The corresponding column may have been deleted in the meantime,
    // so we check that the UID is still valid
    if (!dbout->isUIDValid(iattout_start)) return;

    if (names.empty())
    {
      if (nvar <= 0) nvar = 1;
    }
    else
    {
      nvar = static_cast<Id>(names.size());
    }

    // Create simulation names
    VectorString outnames =
      _createSimulationNames(names, nvar, nbsimu, flagSimuFirst);

    // Set the names in the database
    Id ntotal = nvar * nbsimu;
    for (Id i = 0; i < ntotal; i++)
    {
      dbout->setNameByUID(iattout_start + i, outnames[i]);
    }

    if (flagSetLocator)
    {
      if (_flagLocator && _locatorOutType != ELoc::UNDEFINED)
      {
        // Erase already existing locators of the same Type
        if (_cleanSameLocator && locatorShift == 0)
          dbout->clearLocators(_locatorOutType);

        // Set the locator for all variables
        for (Id i = 0; i < ntotal; i++)
          dbout->setLocatorByUID(
            iattout_start + i, _locatorOutType, i + locatorShift);
      }
    }
  }

  /**
   * @brief Assign the configured locator to output variables.
   *
   * The locator is assigned to a set of output variables starting at
   * `iattout_start`.
   *
   * @param dbout Output Db containing the variables.
   * @param iattout_start Index of the first variable receiving the locator.
   * @param nvar Number of variables.
   * @param nitems Number of items associated with each variable.
   * @param locatorShift Shift applied when assigning the locator.
   */
  void NamingConvention::setLocators(
    Db* dbout,
    Id iattout_start,
    Id nvar,
    Id nitems,
    Id locatorShift) const
  {
    if (!_flagLocator || _locatorOutType == ELoc::UNDEFINED) return;

    // Erase already existing locators of the same Type
    // (this is only done if you are not precisely adding higher order version for given locator)
    if (_cleanSameLocator && locatorShift == 0)
      dbout->clearLocators(_locatorOutType);

    // Set the locator for all variables
    for (Id ecr = 0; ecr < nvar * nitems; ecr++)
      dbout->setLocatorByUID(
        iattout_start + ecr, _locatorOutType, ecr + locatorShift);
  }

  Id NamingConvention::_getNameCount(const VectorString& names, Id nvar)
  {
    if (nvar <= 0)
    {
      // Argument 'nvar' is not defined yet
      if (names.empty()) return 1;
      return static_cast<Id>(names.size());
    }

    // Argument 'nvar' is provided: is it consistent with 'names'
    if (names.empty()) return nvar;

    // Both 'nvar' and 'names' are provided. For safety reasons,
    // the number of variables is the minimum between the two
    return MIN(nvar, static_cast<Id>(names.size()));
  }

  void NamingConvention::_setNames(
    Db* dbout,
    Id iattout_start,
    const VectorString& names,
    Id nvar,
    const String& qualifier,
    Id nitems) const
  {
    auto nloc = _getNameCount(names, nvar);
    VectorString outnames = _createNames(names, nloc, qualifier, nitems);
    correctNamesForDuplicates(outnames, dbout->getAllNames());

    Id ecr = 0;
    for (Id ivar = 0; ivar < nloc; ivar++)
    {
      for (Id item = 0; item < nitems; item++)
      {
        dbout->setNameByUID(iattout_start + ecr, outnames[ecr]);
        ecr++;
      }
    }
  }

  VectorString NamingConvention::_createNames(
    const VectorString& names,
    Id nvar,
    const String& qualifier,
    Id nitems) const
  {
    VectorString outnames;

    for (Id ivar = 0; ivar < nvar; ivar++)
    {
      // Variable 'local' defined for each variable, is:
      // - extracted from the array 'names' (if defined)
      // - generated as the rank of the variable (if several)
      String loc_varname;
      String loc_number;
      if (_flagVarname)
      {
        if (static_cast<Id>(names.size()) == nvar) loc_varname = names[ivar];
        if (loc_varname.empty() && nvar > 1)
        {
          if (Old_Style)
            loc_varname = std::to_string(ivar + 1);
          else
            loc_varname = concatenateString("V", ivar + 1, "");
        }
      }
      else
      {
        // Build the rank from the variable number (possibly overwritten by item number)
        if (nvar > 1)
        {
          if (Old_Style)
            loc_number = std::to_string(ivar + 1);
          else
            loc_number = concatenateString("V", ivar + 1, "");
        }
      }

      for (Id item = 0; item < nitems; item++)
      {
        String loc_qualifier;

        if (_flagQualifier)
        {
          loc_qualifier = qualifier;
          if (nitems > 1)
          {
            if (Old_Style)
              loc_number = std::to_string(item + 1);
            else
              loc_number = concatenateString("S", item + 1, "");
          }
        }

        // Compose the variable name
        String name = concatenateStrings(
          _delim, _prefix, loc_varname, loc_qualifier, loc_number);

        if (name.empty()) name = "Dummy";
        outnames.push_back(name);
      }
    }
    return outnames;
  }

  VectorString NamingConvention::_createSimulationNames(
    const VectorString& names,
    Id nvar,
    Id nbsimu,
    bool flagSimuFirst) const
  {
    if (nvar <= 0 || nbsimu <= 0) return {};

    VectorString outnames;

    // Determine variable names
    VectorString varnames;
    bool isConditional = !names.empty();

    if (isConditional)
    {
      // Conditional simulation: use provided variable names
      varnames = names;
      if (static_cast<Id>(varnames.size()) != nvar && nvar > 0)
        varnames.resize(nvar);
    }
    else
    {
      // Non-conditional simulation: use V1, V2, ... format
      for (Id ivar = 0; ivar < nvar; ivar++)
      {
        String varname = 'V' + std::to_string(ivar + 1);
        varnames.push_back(varname);
      }
    }

    // Create names based on storage order
    if (flagSimuFirst)
    {
      // Simulation varies first: V1.S1, V1.S2, ..., V1.Sn, V2.S1, V2.S2, ...
      for (Id ivar = 0; ivar < nvar; ivar++)
      {
        for (Id isimu = 0; isimu < nbsimu; isimu++)
        {
          String simuname = 'S' + std::to_string(isimu + 1);
          String name =
            concatenateStrings(_delim, _prefix, varnames[ivar], simuname);
          if (name.empty()) name = "Dummy";
          outnames.push_back(name);
        }
      }
    }
    else
    {
      // Variable varies first: V1.S1, V2.S1, ..., Vn.S1, V1.S2, V2.S2, ...
      for (Id isimu = 0; isimu < nbsimu; isimu++)
      {
        for (Id ivar = 0; ivar < nvar; ivar++)
        {
          String simuname = 'S' + std::to_string(isimu + 1);
          String name =
            concatenateStrings(_delim, _prefix, varnames[ivar], simuname);
          if (name.empty()) name = "Dummy";
          outnames.push_back(name);
        }
      }
    }

    return outnames;
  }

  String NamingConvention::toString(const AStringFormat* /*strfmt*/) const
  {
    std::stringstream sstr;

    sstr << toStrTitle(0, "Naming Convention");
    sstr << "- Prefix  = " << _prefix << std::endl;
    sstr << "- Delimitor = '" << _delim << "'" << std::endl;
    sstr << "- Add the Variable Name = " << _flagVarname << std::endl;
    sstr << "- Add the Qualifier     = " << _flagQualifier << std::endl;
    sstr << "- Assign a Locator to output variables = " << _flagLocator
         << std::endl;
    sstr << "- Type of assigned locator = " << _locatorOutType.getDescr()
         << std::endl;
    sstr << "- Clean any other similar locator = " << _cleanSameLocator
         << std::endl;

    return sstr.str();
  }

  /**
   * @brief Generate a variable name according to the naming convention.
   *
   * This static method provides a way to generate a variable name without
   * creating a NamingConvention object.
   *
   * Depending on the arguments, the generated name can include the input
   * variable name, variable rank, simulation rank and an extension.
   *
   * @param prefix Prefix used in the generated name.
   * @param db Db containing the input variable names. May be nullptr when
   * the variable name is not required.
   * @param ivar Index of the input variable.
   * @param nvar Number of variables.
   * @param isimu Index of the simulation.
   * @param nbsimu Number of simulations.
   * @param extension Additional extension appended to the generated name.
   * @param delim Delimiter used to separate the different components.
   * @return Generated variable name.
   */
  String NamingConvention::getNameEncoded(
    const String& prefix,
    const Db* db,
    Id ivar,
    Id nvar,
    Id isimu,
    Id nbsimu,
    const String& extension,
    const String& delim)
  {
    String loc_varname;
    if (db != nullptr)
    {
      if (db->getNLoc(ELoc::Z) > 0)
        loc_varname = db->getNameByLocator(ELoc::Z, ivar);
    }
    else
    {
      if (ivar < 0)
      {
        loc_varname = "*";
      }
      else
      {
        if (nvar > 1 && ivar > 0)
          loc_varname = concatenateString("V", ivar, "");
      }
    }

    String loc_qualifier;
    if (!extension.empty())
    {
      loc_qualifier = extension;
    }
    else
    {
      if (isimu < 0)
      {
        loc_qualifier = "*";
      }
      else
      {
        if (nbsimu > 1 && isimu > 0)
          loc_qualifier = concatenateString("S", isimu, "");
      }
    }

    String name = concatenateStrings(delim, prefix, loc_varname, loc_qualifier);
    if (name.empty()) name = "Dummy";

    return name;
  }

} // namespace gstlrn
