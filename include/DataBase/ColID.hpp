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

#include "Basic/AStringable.hpp"
#include "DataBase/RoleID.hpp"
#include "geoslib_define.h"

namespace gstlrn
{
  /**
   * \brief
   * Class Descriptor
   *
   * This structure enables the user to refer to a given Column of the Data Base using
   * a variety of identification methods:
   * - by its column index
   * - by its name (unique in the Data Base) or
   * - by its Role (name and rank)
   * (the Column identification also handles the version number)
   *
   * Obviously, at the end of the Column recognition, the returned value must be the index of the Column
   * in the Data Base, which is unique and unambiguous.
   *
   * Note that using the Column index is not recommended as it may change over time,
   * and the user should prefer the use of the Column name or RoleID.
   * Nevertheless, it remains the quickest (direct) solution; all other recognition methods
   * require a search throughout all the Columns of the Data Base.
   *
   * In order to increase the security of this Column designation, the concept of the UID is introduced.
   *
   * When a new Column is added to the Data Base, it is assigned a unique identifier (UID) which is never modified
   * during the life of the Data Base. When a Column is deleted, it does with its UID.
   * This UID is stored in the ColID structure when this ColID results from a new Column creation.
   *
   * This UID (when stored in the ColID structure)is used to check that the Column designated by a ColID
   * is indeed the same Column as the one originally created. If the Column designation (using Column Index,
   * Column Name or RoleID) corresponds to a different Column than the one originally created, the UID prevails:
   * the Column actually identified is the one matching the UID rather than any other idenification modes.
   */
  class GSTLEARN_EXPORT ColID: public AStringable
  {
  public:
    ColID(const std::string& name, Id version = 0);
    ColID(const std::string_view name, Id version = 0);
    ColID(Id icol, Id version = 0);
    ColID(const RoleID& roleID, Id version = 0);
    ColID(const ERole& role, Id version = 0);
    // Important notes:
    // - Constructor from Enum does not allow setting the 'index' directly.
    // - Constructor with 3 arguments is discarded in order to avoid
    // ambiguity between "index" and "version".

#ifndef SWIG
    // catch string literals and dispatch them to the string_view constructor
    template<int N>
    ColID(const char (&name)[N], Id version = 0)
      : ColID{std::string_view{name}, version}
    {
    }
#endif

    /// AStringable Interface
    String toString(const AStringFormat* strfmt = nullptr) const override;

    const String& getName() const { return _name; }

    Id getICol() const { return _icol; }

    const RoleID& getRoleID() const { return _roleID; }

    const ERole& getRole() const { return _roleID.getRole(); }

    Id getIndex() const { return _roleID.getIndex(); }

    Id getVersion() const { return _version; }

    void setName(const String& name) { _name = name; }

    void setICol(Id icol) { _icol = icol; }

    void setRoleID(const RoleID& roleID) { _roleID = roleID; }

    void setVersion(Id version) { _version = version; }

    bool hasRoleDefined() const;

    Id getUniqueIndex() const { return _uid; }

    void setUniqueIndex(Id uid) { _uid = uid; }

    /**
     * Factory helpers mainly intended for language bindings (SWIG).
     */
#ifndef SWIG
    static ColID create(std::string_view name, Id version = 0);
    static ColID create(Id icol, Id version = 0);
    static ColID create(const RoleID& roleID, Id version = 0);
    static ColID create(const ColID& colid, Id version = 0);
#endif

  private:
    Id _icol{-1};
    String _name{};
    RoleID _roleID{RoleID(ERole::UNDEFINED, 0)};
    Id _version{0};
    Id _uid{-1};
  };

} // namespace gstlrn
