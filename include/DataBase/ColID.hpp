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

#include "gstlearn_export.hpp"

#include "DataBase/RoleID.hpp"
#include "geoslib_define.h"

namespace gstlrn
{
  class GSTLEARN_EXPORT ColID
  {
  public:
    ColID(const std::string& name, Id version = 0);
    ColID(const std::string_view name, Id version = 0);
    ColID(Id icol, Id version = 0);
    ColID(const RoleID& roleID, Id version = 0);
    ColID(const ERole& role, Id index = 0, Id version = 0);

#ifndef SWIG
    // catch string literals and dispatch them to the string_view constructor
    template<int N>
    ColID(const char (&name)[N], Id version = 0)
      : ColID{std::string_view{name}, version}
    {
    }
#endif

    const String& getName() const { return _name; }

    Id getICol() const { return _icol; }

    const RoleID& getRoleID() const { return _roleID; }

    const ERole& getRole() const { return _roleID.getRole(); }

    Id getIndex() const { return _roleID.getIndex(); }

    Id getVersion() const { return _version; }

    String getDescr() const;

    void setName(const String& name) { _name = name; }

    void setICol(Id icol) { _icol = icol; }

    void setRoleID(const RoleID& roleID) { _roleID = roleID; }

    void setVersion(Id version) { _version = version; }

    bool hasRoleDefined() const;

    /**
     * Factory helpers mainly intended for language bindings (SWIG).
     */
#ifndef SWIG
    static ColID create(std::string_view name, Id version = 0);

    static ColID create(Id icol, Id version = 0);

    static ColID create(const RoleID& roleID, Id version = 0);

    static ColID create(const ERole& role, Id index = 0, Id version = 0);

    static ColID create(const ColID& colid, Id version = 0);
#endif

  private:
    Id _icol{-1};
    String _name{};
    RoleID _roleID{RoleID(ERole::UNDEFINED, 0)};
    Id _version{0};
  };

} // namespace gstlrn
