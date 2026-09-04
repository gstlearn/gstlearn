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

// Next line is needed to overload getNA and isNA for Category
#include "Basic/Utilities.hpp" // IWYU pragma: keep
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
  /**
   * @brief Category identifier and label class.
   *
   * Stores a category identifier (Id) and its associated textual label.
   */
  class GSTLEARN_EXPORT Category
  {
  public:
    Category() = default;

    Category(Id id, String name)
      : _id(id)
      , _name(std::move(name))
    {
    }

    Id getId() const { return _id; }

    const String& getName() const { return _name; }

    void setId(Id id) { _id = id; }

    void setName(const String& name) { _name = name; }

    bool operator==(const Category& other) const
    {
      return _id == other._id && _name == other._name;
    }

    Id getIdentifier() const { return _id; }

    bool operator!=(const Category& other) const { return !(*this == other); }

  private:
    Id _id{-1};
    String _name{};
  };

  template<>
  inline Category getNA<Category>()
  {
    return Category{-1, ""};
  }

  template<>
  inline bool isNA(const Category& v)
  {
    return v.getId() < 0;
  }

} // namespace gstlrn
