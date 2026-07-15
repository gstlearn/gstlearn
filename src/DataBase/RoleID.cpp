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
#include "DataBase/RoleID.hpp"
#include <algorithm>
#include <cctype>

namespace gstlrn
{
  namespace
  {
    // Helper accepting any string type (String, std::string_view, const char*)
    String toLowerCopy(std::string_view str)
    {
      String res(str);
      std::transform(
        res.begin(), res.end(), res.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
      return res;
    }
  } // namespace

  RoleID::RoleID(const ERole& role, Id index)
    : _role(role)
    , _index(index)
  {
  }

  bool RoleID::isUnique() const
  {
    auto it = ERoleAttr.find(_role.getKey());
    return (it != ERoleAttr.end()) ? it->second.isUnique : false;
  }

  String RoleID::getName() const
  {
    if (!isDefined()) return STRING_NA;

    String name = toLowerCopy(_role.getKey());

    if (!isUnique()) name += std::to_string(_index + 1);

    return name;
  }

  bool RoleID::isEqual(const RoleID& roleID, bool checkIndex) const
  {
    if (_role.isDifferent(roleID.getRole())) return false;
    return !checkIndex || (_index == roleID.getIndex());
  }

  std::optional<RoleID> RoleID::createFromName(String string)
  {
    if (string.empty()) return RoleID();

    toLower(string);

    // Separate alphabetic characters (role) from numeric characters (index)
    String roleName;
    String indexString;

    for (char c: string)
    {
      if (std::isdigit(static_cast<unsigned char>(c)))
        indexString += c;
      else
        roleName += c;
    }

    // Search for the corresponding role in the enumeration
    ERole role = ERole::UNDEFINED;
    auto it = ERole::getIterator();
    while (it.hasNext())
    {
      auto current = *it;
      if (current != ERole::UNDEFINED
          && toLowerCopy(current.getKey()) == roleName)
      {
        role = current;
        break;
      }
      it.toNext();
    }

    if (role == ERole::UNDEFINED) return RoleID();

    // Compute the internal index (1-based to 0-based)
    Id index = 0;
    if (!indexString.empty())
      index = std::max(std::atoi(indexString.c_str()) - 1, 0);

    // Check the role uniqueness constraint
    RoleID result(role, index);
    if (result.isUnique() && index > 0) return std::nullopt;

    return result;
  }

  RoleID roleIDIdentify(const String& name)
  {
    auto res = RoleID::createFromName(name);
    return res.value_or(RoleID());
  }

  void RoleID::removeRole()
  {
    _role = ERole::UNDEFINED;
    _index = 0;
  }

} // namespace gstlrn
