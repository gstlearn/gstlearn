/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2025) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Basic/VectorNumT.hpp"
#include "Enum/ERole.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include <initializer_list>

namespace gstlrn
{
  /**
   * @brief Thie class describes a set of Columns in a DbData. They can be described by:
   * - their column indexes
   * - their names (unique in the Data Base)
   * - their Role (name and rank)
   * Note that the Column Name(s) can take the regexp syntax into account,
   * allowing to select a set of Columns by their name(s) adding the '*' special symbol.
   *
   * Two special entries have been provided:
   * - using a single Name, as this name can be inflated to many names using the '*' special symbol
   * - using a single Column Index: this entry is weird as it only corresponds to a single column.
   * Nevertheless, it is provided for consistency with the other entries, and to allow a single column selection.
   *
   * The aim of this class is to return internally a list of ColID(s), where one Colid is asggined
   * to each one of the targeted Columns.
   */
  class GSTLEARN_EXPORT MultID
  {
  public:
    MultID(const String& name);
    MultID(const VectorString& names);
    MultID(const VectorInt& icols);
    MultID(const ERole& role);
    MultID(Id icol);

#ifndef SWIG
    template<int N>
    MultID(const char (&name)[N])
      : MultID(String(name))
    {
    }

    MultID(std::initializer_list<const char*> names);
    MultID(std::initializer_list<String> names);
    MultID(std::initializer_list<Id> icols);
#endif

    const VectorString& getNames() const { return _names; }

    const VectorInt& getICols() const { return _icols; }

    const ERole& getRole() const { return _role; }

    String getDescr() const;

    void setNames(const VectorString& names) { _names = names; }

    void setICols(const VectorInt& icols) { _icols = icols; }

    void setERole(const ERole& role) { _role = role; }

    bool hasRoleDefined() const;

    /**
     * Factory helpers mainly intended for language bindings (SWIG).
     */
#ifndef SWIG
    static MultID create(const String& name);
    static MultID create(const VectorString& names);
    static MultID create(const VectorInt& icols);
    static MultID create(const ERole& role);
    static MultID create(const MultID& MultID);
#endif

  private:
    VectorInt _icols{};
    VectorString _names{};
    ERole _role{ERole::UNDEFINED};
  };

} // namespace gstlrn
