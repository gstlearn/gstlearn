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

  ColID::ColID(const ERole& role, Id index, Id version)
    : _icol(-1)
    , _name()
    , _roleID(role, index)
    , _version(version)
  {
  }

  /**
   * @brief Constructs a String which gives all possible information on a Column Identifier
   *
   * @return String
   */
  String ColID::getDescr() const
  {
    String descr;

    if (!_name.empty()) descr += "Name: " + _name + ", ";

    if (_icol >= 0) descr += "Index: " + std::to_string(_icol) + ", ";

    if (hasRoleDefined()) descr += "Role: " + _roleID.getName() + ", ";

    descr += "Version: " + std::to_string(_version);

    return descr;
  }

  bool ColID::hasRoleDefined() const
  {
    return _roleID.getRole().isDifferent(ERole::UNDEFINED);
  }

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

} // namespace gstlrn
