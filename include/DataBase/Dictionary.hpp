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

#include "geoslib_define.h"

#include <map>
#include <optional>
#include <string_view>
#include <utility>

namespace gstlrn
{
  /**
   * @brief Dictionary associating category identifiers with their labels.
   *
   * A dictionary stores a set of categories identified by an integer key.
   * Each category is associated with a textual label.
   *
   * The category is represented by the @c Category type, which is a pair
   * containing the category identifier and its label.
   */
  class Dictionary
  // class GSTLEARN_EXPORT Dictionary
  {
  public:
    /**
     * @brief Category identifier and label.
     *
     * The first element is the category identifier and the second element
     * is its associated label.
     */
    using Category = std::pair<Id, std::string_view>;

    Dictionary() = default;

    /**
     * @brief Constructs a dictionary from an existing map.
     *
     * The map is moved into the dictionary.
     *
     * @param data Map associating category identifiers with their labels.
     */
    Dictionary(std::map<Id, String>&& data)
      : _data{std::move(data)}
    {
    }

    bool addCategory(const Id key, const String& val);

    bool hasCategory(const Category& cat) const;

#ifndef SWIG
    /**
     * @brief Retrieves a category by its identifier.
     *
     * @param key Category identifier.
     *
     * @return The corresponding category, or an empty optional if the
     *         identifier is not present in the dictionary.
     */
    std::optional<Category> operator[](const Id key) const
    {
      const auto it = this->_data.find(key);
      if (it == this->_data.end()) return {};
      return *it;
    }
#endif

  private:
    std::map<Id, String> _data;
  };

} // namespace gstlrn
