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
#include "DataBase/MultID.hpp"

namespace gstlrn
{
  MultID::MultID(const String& name)
    : _icols()
    , _names()
    , _role(ERole::UNDEFINED)
  {
    _names.push_back(name);
  }

  MultID::MultID(const VectorInt& icols)
    : _icols(icols)
    , _names()
    , _role(ERole::UNDEFINED)
  {
  }

  MultID::MultID(const VectorString& names)
    : _icols()
    , _names(names)
    , _role(ERole::UNDEFINED)
  {
  }

  MultID::MultID(const ERole& role)
    : _icols()
    , _names()
    , _role(role)
  {
  }

  MultID::MultID(Id icol)
    : _icols()
    , _names()
    , _role(ERole::UNDEFINED)
  {
    _icols.push_back(icol);
  }

  /**
   * @brief Constructs a String which gives all possible information on a Column Identifier
   *
   * @return String
   */
  String MultID::getDescr() const
  {
    String descr;

    if (!_names.empty())
    {
      descr = "Names: ";
      for (const auto& name: _names)
      {
        descr += name + " ";
      }
      return descr;
    }

    if (!_icols.empty())
    {
      descr += "Column Indices: ";
      for (const auto& icol: _icols)
      {
        descr += std::to_string(icol) + " ";
      }
      return descr;
    }

    if (hasRoleDefined())
    {
      descr = "Role: ";
      descr += _role.getKey();
      descr += ", ";
    }
    return descr;
  }

  bool MultID::hasRoleDefined() const
  {
    return _role.isDifferent(ERole::UNDEFINED);
  }

  MultID MultID::create(const String& name)
  {
    return MultID(name);
  }

  MultID MultID::create(const VectorInt& icols)
  {
    return MultID(icols);
  }

  MultID MultID::create(const ERole& role)
  {
    return MultID(role);
  }

  MultID MultID::create(const MultID& multid)
  {
    return MultID(multid);
  }

} // namespace gstlrn
