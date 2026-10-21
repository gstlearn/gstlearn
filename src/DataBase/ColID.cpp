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
#include "DataBase/ColID.hpp"

namespace gstlrn
{
  // Kept for SWIG, which currently resolves this overload more reliably than
  // the std::string_view constructor.
  ColID::ColID(const std::string& name, Id version)
    : ColID(std::string_view(name), version)
  {
  }

  /**
   * @brief Construct a new ColID from a Name stored as char*
   *
   * @param name
   * @param version
   */
  ColID::ColID(const std::string_view name, Id version)
    : _icol(-1)
    , _name(name)
    , _roleID()
    , _version(version)
  {
  }

  /**
   * @brief Construct a new ColID from a Column Index (0 based)
   *
   * @param icol Column Index
   * @param version
   */
  ColID::ColID(Id icol, Id version)
    : _icol(icol)
    , _name()
    , _roleID()
    , _version(version)
  {
  }

  ColID::ColID(const RoleID& roleID, Id version)
    : _icol(-1)
    , _name()
    , _roleID(roleID)
    , _version(version)
  {
  }

#ifndef SWIG
  ColID::ColID(const ERole& role, Id index, Id version)
    : _icol(-1)
    , _name()
    , _roleID(role, index)
    , _version(version)
  {
  }
#endif

  bool ColID::hasRoleDefined() const
  {
    return _roleID.getRole().isDifferent(ERole::UNDEFINED);
  }

#ifndef SWIG
  ColID ColID::create(std::string_view name, Id version)
  {
    return ColID(name, version);
  }

  ColID ColID::create(Id icol, Id version)
  {
    return ColID(icol, version);
  }

  ColID ColID::create(const RoleID& roleID, Id version)
  {
    return ColID(roleID, version);
  }

  ColID ColID::create(const ERole& role, Id index, Id version)
  {
    return ColID(role, index, version);
  }

  ColID ColID::create(const ColID& colid, Id version)
  {
    ColID res(colid);
    if (version != 0) res.setVersion(version);
    return res;
  }
#endif

  String ColID::toString(const AStringFormat* strfmt) const
  {
    DECLARE_UNUSED(strfmt);
    std::stringstream sstr;

    sstr << "Name: " << _name << "\n";
    sstr << "Column Index: " << _icol << "\n";
    sstr << "Role: " << _roleID.getName() << "\n";
    sstr << "Version: " << _version << "\n";

    return sstr.str();
  }

} // namespace gstlrn
